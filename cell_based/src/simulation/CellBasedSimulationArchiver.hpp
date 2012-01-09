/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef CELLBASEDSIMULATIONARCHIVER_HPP_
#define CELLBASEDSIMULATIONARCHIVER_HPP_

// Must be included before any other serialization headers
#include "CheckpointArchiveTypes.hpp"

#include <string>
#include <iostream>

#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"
#include "FileFinder.hpp"

/**
 * CellBasedSimulationArchiver handles the checkpointing (saving and loading)
 * of all the various AbstractCellBasedSimulation objects. It has no explicit constructor
 * (just uses a default one) and no member variables.
 */
template<unsigned DIM, class SIM>
class CellBasedSimulationArchiver
{
public:

    /**
     * Loads a saved cell-based simulation to run further.
     *
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


template<unsigned DIM, class SIM>
SIM* CellBasedSimulationArchiver<DIM, SIM>::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
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
    ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_filename);
    boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

    // Load the simulation
    SIM* p_sim;
    (*p_arch) >> p_sim;
    return p_sim;
}

template<unsigned DIM, class SIM>
void CellBasedSimulationArchiver<DIM, SIM>::Save(SIM* pSim)
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
    ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_filename);
    boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

    // Archive the simulation (const-ness would be a pain here)
    (*p_arch) & pSim;
}

#endif /*CELLBASEDSIMULATIONARCHIVER_HPP_*/
