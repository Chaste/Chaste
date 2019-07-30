/*

Copyright (c) 2005-2019, University of Oxford.
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
     * @return a pointer to the unarchived cardiac problem class
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
     * @return a pointer to the unarchived cardiac problem class
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
     * @return a pointer to the migrated cardiac problem class
     */
    static PROBLEM_CLASS* Migrate(const FileFinder& rDirectory);
};

#endif /*CARDIACSIMULATIONARCHIVER_HPP_*/
