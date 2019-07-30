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

#ifndef TESTARCHIVEFORMAT_HPP_
#define TESTARCHIVEFORMAT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include <iomanip>
#include <boost/foreach.hpp>
#include "CryptSimulation2d.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "CellMutationStatesCountWriter.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * This class consists of a single crypt simulation archiving test.
 */
class TestArchiveFormat : public CxxTest::TestSuite
{
public:

    /**
     * This test is required because Test2DCryptRepresentativeSimulation loads
     * an archive stored in crypt/test/data. When the archiving of
     * OffLatticeSimulation and associate classes is updated, the stored archive
     * needs to be updated. This test checks that the archive can be loaded,
     * and will seg fault if not. It does nothing more, so it runs quickly
     * and can be in the continuous test pack.
     *
     * IF THIS TEST FAILS:
     * - You have probably changed an archiving method somewhere
     * - You need to remake crypt/test/data/<test below>/archive/
     * - To do this re-run TestGenerateSteadyStateCrypt.hpp
     * - Archives produced can then be copied to crypt/test/data/<test below>/archive/
     *
     * Note that when updating the archive, you can run TestGenerateSteadyStateCrypt.hpp with build=GccOpt to speed up the test.
     *
     * Note: from Chaste release 3.3 onward the earliest version of Boost supported is 1.40.
     *
     * NB: Produce archives with
     *  scons build=GccOpt_hostconfig,boost=1-40,use-cvode=0 test_suite=crypt/test/simulation/TestGenerateSteadyStateCrypt.hpp
     *  cp /tmp/$USER/testoutput/SteadyStateCrypt/archive/?*_150.* crypt/test/data/SteadyStateCrypt/archive/
     *
     * OR to produce archives in CMake:
     *  cmake -DBOOST_ROOT=/path/to/boost1.40 -DChaste_USE_CVODE=OFF /path/to/Chaste
     *  make TestGenerateSteadyStateCrypt_simulation_Runner
     *  ctest -R TestGenerateSteadyStateCrypt
     *  cp /path/to/Chaste/testoutput/SteadyStateCrypt/archive/?*_150.* /path/to/Chaste/crypt/test/data/SteadyStateCrypt/archive/
     *
     */
    void TestLoadArchive()
    {
        // Set start time
        SimulationTime::Instance()->SetStartTime(0.0);

        {
            // Horrendous hack to get this test to pass when using static libraries!
            Cylindrical2dMesh* p_mesh = new Cylindrical2dMesh(1.0);
            // Each of the following classes gives an 'unregistered class' error if not instantiated here...
            MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh);
            StochasticWntCellCycleModel ccm;
            SloughingCellKiller<2> killer(&crypt, 1.0);
            GeneralisedLinearSpringForce<2> force;
            crypt.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        }

        // Directory in which the stored results were archived
        std::string test_to_profile = "SteadyStateCrypt";

        // Simulation time at which the stored results were archived
        double t = 150;

        // The archive must be copied from crypt/test/data/<test_to_profile>
        FileFinder test_data_directory("crypt/test/data/" + test_to_profile + "/archive", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(test_data_directory.IsDir());

        // to the testoutput/archive directory to continue running the simulation
        OutputFileHandler archive_handler(test_to_profile + "/archive");

        // Following is done in two lines to avoid a bug in Intel compiler v12.0
        std::vector<FileFinder> temp_files = test_data_directory.FindMatches("*");
        BOOST_FOREACH(FileFinder temp_file, temp_files)
        {
            archive_handler.CopyFileTo(temp_file);
        }
#ifndef CHASTE_CVODE
        std::cout << "Warning: CVODE is off.  If this configuration of the test suite fails, but none of the others, then the archive was built with CVODE on.  If should be built with CVODE off." << std::endl;
#endif //CHASTE_CVODE

        // Load and run crypt simulation
        CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load(test_to_profile, t);
        p_simulator->SetEndTime(t + 1);

        ///\todo we should also test the rest of the simulation setup

        // Tidy up
        delete p_simulator;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTARCHIVEFORMAT_HPP_*/
