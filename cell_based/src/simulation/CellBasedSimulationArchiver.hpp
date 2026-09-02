/*

Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CELLBASEDSIMULATIONARCHIVER_HPP_
#define CELLBASEDSIMULATIONARCHIVER_HPP_

// Must be included before any other serialization headers
#include "CheckpointArchiveTypes.hpp"

#include <string>

// The headers below are needed by the templated implementation, which is in this header
#include <fstream>
#include <sstream>

#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"
#include "FileFinder.hpp"
#include "SimulationTime.hpp"

/**
 * CellBasedSimulationArchiver handles the checkpointing (saving and loading)
 * of all the various AbstractCellBasedSimulation objects. It has no explicit constructor
 * (just uses a default one) and no member variables.
 *
 * The template options for INPUT_ARCHIVE_TYPE and OUTPUT_ARCHIVE_TYPE can (in theory) be set independently but should be chosen accordingly.
 * To help with this, below, we define type aliases to set both archive types for text and binary I/O respectively.
 * Default settings are chosen such that existing implementations depending on text format as default are kept intact.
 */
template <unsigned ELEMENT_DIM, class SIM, unsigned SPACE_DIM = ELEMENT_DIM, class INPUT_ARCHIVE_TYPE = boost::archive::text_iarchive, class OUTPUT_ARCHIVE_TYPE = boost::archive::text_oarchive>
class CellBasedSimulationArchiver
{
public:
    /**
     * Loads a saved cell-based simulation to run further.
     *
     * @return the unarchived simulation object
     * @param rArchiveDirectory  the name of the simulation to load
     *   (specified originally by simulation.SetOutputDirectory("wherever"); )
     * @param rTimeStamp  the time at which to load the simulation (this must
     *   be one of the times at which simulation.Save() was called)
     */
    static SIM* Load(const std::string& rArchiveDirectory, const double& rTimeStamp);

    /**
     * Saves the whole cell-based simulation for restarting later.
     *
     * Puts it in the archive folder under the simulation's OutputDirectory,
     * in the file "cell_population_sim_at_time_<SIMULATION TIME>.arch".
     * The mesh is written to files in the same folder.
     *
     * First archives simulation time (and other singletons, if used)
     * then the simulation itself.
     *
     * @param pSim pointer to the simulation
     */
    static void Save(SIM* pSim);
};

template <unsigned ELEMENT_DIM, class SIM, unsigned SPACE_DIM, class INPUT_ARCHIVE_TYPE, class OUTPUT_ARCHIVE_TYPE>
SIM* CellBasedSimulationArchiver<ELEMENT_DIM, SIM, SPACE_DIM, INPUT_ARCHIVE_TYPE, OUTPUT_ARCHIVE_TYPE>::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
{
    /**
     * Find the right archive (and mesh) to load.  The files are contained within
     * the 'archive' folder in rArchiveDirectory, with the archive itself called
     * 'cell_population_sim_at_time_`rTimeStamp`.arch'.  The path to this file is returned.
     *
     * The path to the mesh is stored in ArchiveLocationInfo for use by the
     * CellPopulation de-serialization routines.
     */
    std::ostringstream time_stamp;
    time_stamp << rTimeStamp;
    std::string archive_filename = "cell_population_sim_at_time_" + time_stamp.str() + ".arch";
    std::string mesh_filename = "mesh_" + time_stamp.str();
    FileFinder archive_dir(rArchiveDirectory + "/archive/", RelativeTo::ChasteTestOutput);
    ArchiveLocationInfo::SetMeshPathname(archive_dir, mesh_filename);

    // Create an input archive
    ArchiveOpener<INPUT_ARCHIVE_TYPE, std::ifstream> arch_opener(archive_dir, archive_filename);
    INPUT_ARCHIVE_TYPE* p_arch = arch_opener.GetCommonArchive();

    // Load the simulation
    SIM* p_sim;
    (*p_arch) >> p_sim;
    return p_sim;
}

template <unsigned ELEMENT_DIM, class SIM, unsigned SPACE_DIM, class INPUT_ARCHIVE_TYPE, class OUTPUT_ARCHIVE_TYPE>
void CellBasedSimulationArchiver<ELEMENT_DIM, SIM, SPACE_DIM, INPUT_ARCHIVE_TYPE, OUTPUT_ARCHIVE_TYPE>::Save(SIM* pSim)
{
    // Get the simulation time as a string
    const SimulationTime* p_sim_time = SimulationTime::Instance();
    assert(p_sim_time->IsStartTimeSetUp());
    std::ostringstream time_stamp;
    time_stamp << p_sim_time->GetTime();

    // Set up folder and filename of archive
    FileFinder archive_dir(pSim->GetOutputDirectory() + "/archive/", RelativeTo::ChasteTestOutput);
    std::string archive_filename = "cell_population_sim_at_time_" + time_stamp.str() + ".arch";
    ArchiveLocationInfo::SetMeshFilename(std::string("mesh_") + time_stamp.str());

    // Create output archive
    ArchiveOpener<OUTPUT_ARCHIVE_TYPE, std::ofstream> arch_opener(archive_dir, archive_filename);
    OUTPUT_ARCHIVE_TYPE* p_arch = arch_opener.GetCommonArchive();

    // Archive the simulation (const-ness would be a pain here)
    (*p_arch) & pSim;
}

/*
 * Specialization of the CellBasedSimulationArchiver for text based formats to reduce requirement to set both archive type template parameters
 */
template <unsigned ELEMENT_DIM, class SIM, unsigned SPACE_DIM = ELEMENT_DIM>
using CellBasedSimulationArchiverText = CellBasedSimulationArchiver<ELEMENT_DIM, SIM, SPACE_DIM, boost::archive::text_iarchive, boost::archive::text_oarchive>;

/*
 * Specialization of the CellBasedSimulationArchiver for binary based formats to reduce requirement to set both archive type template parameters
 */
template <unsigned ELEMENT_DIM, class SIM, unsigned SPACE_DIM = ELEMENT_DIM>
using CellBasedSimulationArchiverBinary = CellBasedSimulationArchiver<ELEMENT_DIM, SIM, SPACE_DIM, boost::archive::binary_iarchive, boost::archive::binary_oarchive>;

#endif /*CELLBASEDSIMULATIONARCHIVER_HPP_*/
