/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef CARDIACSIMULATIONARCHIVER_HPP_
#define CARDIACSIMULATIONARCHIVER_HPP_

#include <string>

#include "FileFinder.hpp"


/**
 * CardiacSimulationArchiver is a helper class for checkpointing of cardiac simulations.
 *
 * The class is templated over the class defining the simulation (e.g. MonodomainProblem).
 * It would be more natural to template over dimensions as for other classes, and just deal with pointers
 * to AbstractCardiacProblem.  However, this breaks archive backwards compatibility...
 */
template<class PROBLEM_CLASS>
class CardiacSimulationArchiver
{
public:
    /**
     * Archives a simulation in the directory specified.
     *
     * @note Must be called collectively, i.e. by all processes.
     *
     * @param rSimulationToArchive object defining the simulation to archive
     * @param rDirectory directory where the multiple files defining the checkpoint will be stored
     *     (relative to CHASTE_TEST_OUTPUT)
     * @param clearDirectory whether the directory needs to be cleared or not.
     */
    static void Save(PROBLEM_CLASS& rSimulationToArchive, const std::string& rDirectory, bool clearDirectory=true);


    /**
     * Unarchives a simulation from the directory specified.
     *
     * Does a migrate if necessary (this is actually just a wrapper around the
     * Migrate method now).
     *
     * @note Must be called collectively, i.e. by all processes.
     *
     * @param rDirectory directory where the multiple files defining the checkpoint are located
     *     (relative to CHASTE_TEST_OUTPUT)
     */
    static PROBLEM_CLASS* Load(const std::string& rDirectory);

    /**
     * Unarchives a simulation from the directory specified.
     *
     * Does a migrate if necessary (this is actually just a wrapper around the
     * Migrate method now).
     *
     * @note Must be called collectively, i.e. by all processes.
     *
     * @param rDirectory directory where the multiple files defining the checkpoint are located
     */
    static PROBLEM_CLASS* Load(const FileFinder& rDirectory);

    /**
     * Load a simulation, saved by any number of processes, into the correct
     * number of processes for those currently launched.
     *
     * @note Must be called collectively, i.e. by all processes.
     *
     * Uses the DistributedVectorFactory saved in the process 0 archive to work out
     * how many secondary archive files to read, and loads the cells and boundary
     * conditions from these too.
     *
     * Uses a dumb partition to work out how to distribute the mesh and cells over
     * the processes.  If we are loading on the same number of processes as the
     * simulation was saved on, it uses exactly the same distribution as before.
     *
     * @param rDirectory directory where the multiple files defining the checkpoint are located
     */
    static PROBLEM_CLASS* Migrate(const FileFinder& rDirectory);
};

#endif /*CARDIACSIMULATIONARCHIVER_HPP_*/
