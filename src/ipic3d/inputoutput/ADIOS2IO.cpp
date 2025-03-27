#ifdef USE_ADIOS2

#include "ADIOS2IO.hpp"
#include "iPic3D.h"
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Particles3D.h"
#include "Collective.h"

#include "mpi.h"
#include "adios2.h"

#include "ipic3d_cali.h"

#include "MPIdata.h"

#include <chrono>

#define CALI_MARK_ADIOS_OPEN_BEGIN CALI_MARK_BEGIN("adios_open_file")
#define CALI_MARK_ADIOS_OPEN_END CALI_MARK_END("adios_open_file")
#define CALI_MARK_ADIOS_CLOSE_BEGIN CALI_MARK_BEGIN("adios_close_file")
#define CALI_MARK_ADIOS_CLOSE_END CALI_MARK_END("adios_close_file")

namespace ADIOS2IO {

using namespace std;


void ADIOS2Manager::initOutputFiles(string fieldTag, string particleTag, int sample, iPic3D::c_Solver& KCode) {

    CALI_CXX_MARK_SCOPE("adios_mgr_init_output_files");
    if (open) {
        closeOutputFiles();
    }

    this->cartisianRank = KCode.vct->getCartesian_rank();
    this->saveDirName = KCode.col->getSaveDirName();
    this->restartDirName = KCode.col->getRestartDirName();
    this->restartTag = KCode.col->getRestartOutputCycle()? "proc_topology+E+B+rhos+Js+pressure+position+velocity+q+ID"s : ""s;

    this->fieldTag = fieldTag;
    this->particleTag = particleTag;
    this->sample = sample;

    this->col = KCode.col;
    this->vct = KCode.vct;
    this->grid = KCode.grid;
    this->EMf = KCode.EMf;

    this->part = KCode.outputPart;
    this->ns = KCode.col->getNs();
    this->testpart = KCode.testpart;
    this->nstestpart = KCode.col->getNsTestPart();


    // ADIOS2
    this->adios = adios2::ADIOS(MPIdata::get_PicGlobalComm());

    // open files
    if (!fieldTag.empty()) { throw runtime_error("Field output is not supported yet"); 
        CALI_CXX_MARK_SCOPE("adios_mgr_init_field_output_file");
        this->ioField = adios.DeclareIO("FieldOutput");
        this->ioField.SetEngine("BP5");
        auto filePath = saveDirName + "/field_" + to_string(cartisianRank) + ".bp";
        CALI_MARK_ADIOS_OPEN_BEGIN;
        engineField = ioField.Open(filePath, adios2::Mode::Write);
        CALI_MARK_ADIOS_OPEN_END;

    }

    if (!particleTag.empty()) {
        CALI_CXX_MARK_SCOPE("adios_mgr_init_particle_output_file");
        this->ioParticle = adios.DeclareIO("ParticleOutput");
        this->ioParticle.SetEngine("BP5");
        auto filePath = saveDirName + "/particle_" + to_string(cartisianRank) + ".bp";
        CALI_MARK_ADIOS_OPEN_BEGIN;
        engineParticle = ioParticle.Open(filePath, col->getRestart_status() == 0 ? adios2::Mode::Write : adios2::Mode::Append, MPI_COMM_SELF);
        CALI_MARK_ADIOS_OPEN_END;

        // parse the tag and prepae the map
        particleTag.erase(remove(particleTag.begin(), particleTag.end(), ' '), particleTag.end());
        vector<string> tags;
        stringstream ss(particleTag);
        string tag;
        while (getline(ss, tag, '+')) {
            tags.push_back(tag);
        }

        // find the function in the map and register it to vector
        for (auto tag : tags) {
            if (outputTagOptions.find(tag) != outputTagOptions.end()) {
                particleOptions.push_back(outputTagOptions[tag]);
            } else {
                throw runtime_error("Particle output tag is not supported: " + tag);
            }
        }


    }


    if (!restartTag.empty()) { 
        
        CALI_CXX_MARK_SCOPE("adios_mgr_init_restart_output_file");
        this->ioRestart = adios.DeclareIO("RestartOutput");
        this->ioRestart.SetEngine("BP5");
        auto filePath = restartDirName + "/restart_" + to_string(cartisianRank) + ".bp";
        CALI_MARK_ADIOS_OPEN_BEGIN;
        engineRestart = ioRestart.Open(filePath, col->getRestart_status() == 0 ? adios2::Mode::Write : adios2::Mode::Append, MPI_COMM_SELF);
        CALI_MARK_ADIOS_OPEN_END;

        // parse the tag and prepae the map
        restartTag.erase(remove(restartTag.begin(), restartTag.end(), ' '), restartTag.end());
        vector<string> tags;
        stringstream ss(restartTag);
        string tag;
        while (getline(ss, tag, '+')) {
            tags.push_back(tag);
        }
        // find the function in the map and register it to vector
        for (auto tag : tags) {
            if (outputTagOptions.find(tag) != outputTagOptions.end()) {
                restartOptions.push_back(outputTagOptions[tag]);
            } else {
                throw runtime_error("Restart output tag is not supported: " + tag);
            }
        }

    }


    open = true;

}



void ADIOS2Manager::appendFieldOutput(int cycle) {
    throw runtime_error("Field output is not supported yet");
}

void ADIOS2Manager::appendParticleOutput(int cycle) {
    if (particleOptions.empty()) return;

    CALI_MARK_BEGIN("adios_mgr_append_particle_output");
    engineParticle.BeginStep();

    auto cycleVar = _variableHelper<int>(ioParticle, "cycle");
    engineParticle.Put<int>(cycleVar, cycle);

    // auto timeVar = _variableHelper<int>(ioParticle, "IOTimeMS");
    // auto start = chrono::high_resolution_clock::now();

    for (auto option : particleOptions) {
        option(ioParticle, engineParticle);
    }

    // engineParticle.PerformPuts(); // do the heavy job here

    // auto stop = chrono::high_resolution_clock::now();
    // auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    // engineParticle.Put<int>(timeVar, duration.count());

    engineParticle.EndStep();
    CALI_MARK_END("adios_mgr_append_particle_output");

}


void ADIOS2Manager::appendRestartOutput(int cycle) {

    CALI_MARK_BEGIN("adios_mgr_append_restart_output");
    engineRestart.BeginStep();

    auto cycleVar = _variableHelper<int>(ioRestart, "cycle");
    engineRestart.Put<int>(cycleVar, cycle);
    for (auto option : restartOptions) {
        option(ioRestart, engineRestart);
    }

    engineRestart.EndStep();
    CALI_MARK_END("adios_mgr_append_restart_output");

}


void ADIOS2Manager::appendOutput(int cycle) {
    if (!open) throw runtime_error("Output files are not open");


    appendParticleOutput(cycle);

    outputCount++;
    lastCycle = cycle;

}


void ADIOS2Manager::closeOutputFiles() {
    CALI_CXX_MARK_SCOPE("adios_mgr_close_output_files");

    if (!open) return;

    if (!fieldTag.empty()) {
        CALI_MARK_ADIOS_CLOSE_BEGIN;
        engineField.Close();
        CALI_MARK_ADIOS_CLOSE_END;
    }

    if (!particleTag.empty()) {
        CALI_MARK_ADIOS_CLOSE_BEGIN;
        engineParticle.Close();
        CALI_MARK_ADIOS_CLOSE_END;
    }

    if (!restartTag.empty()) {
        CALI_MARK_ADIOS_CLOSE_BEGIN;
        engineRestart.Close();
        CALI_MARK_ADIOS_CLOSE_END;
    }

    open = false;
}







} // end namespace ADIOS2IO

#endif