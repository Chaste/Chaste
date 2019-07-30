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

#ifndef TESTREPRESENTATIVESIMULATION_HPP_
#define TESTREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/foreach.hpp>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "SloughingCellKiller.hpp"

// Need to include all (or at least some) of these files!

#include "GeneralisedLinearSpringForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"

// Needed for NodesOnlyMesh
#include "PetscSetupAndFinalize.hpp"

/**
 * This class consists of a single test, in which a 2D model
 * of a colorectal crypt with representative parameter values
 * is loaded from an archive and simulated for a further period
 * of time.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class Test2DCryptRepresentativeSimulation : public CxxTest::TestSuite
{
public:

    void TestRepresentativeSimulationForProfiling()
    {
        EXIT_IF_PARALLEL;

        // Set start time
        SimulationTime::Instance()->SetStartTime(0.0);

        // Directory in which the stored results were archived
        std::string test_to_load = "SteadyStateCrypt";

        // Simulation time at which the stored results were archived
        double t = 150;

        // Directory in which to store profiling results
        std::string test_to_profile = "CryptProfiling";

        // How long to run the loaded crypt simulation for (in hours)
        double run_for = 10;

        // The archive must be copied from crypt/test/data/<test_to_profile>
        FileFinder test_data_directory("crypt/test/data/" + test_to_load +"/archive",RelativeTo::ChasteSourceRoot);
        TS_ASSERT(test_data_directory.IsDir());

        // to the testoutput/archive directory to continue running the simulation
        OutputFileHandler archive_handler(test_to_profile + "/archive");

        // Following is done in two lines to avoid a bug in Intel compiler v12.0
        std::vector<FileFinder> temp_files = test_data_directory.FindMatches("*");
        BOOST_FOREACH(FileFinder temp_file, temp_files)
        {
            archive_handler.CopyFileTo(temp_file);
        }

        // Load and run crypt simulation
        CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load(test_to_profile,t);
        p_simulator->SetEndTime(t+run_for); // start time + duration
        p_simulator->Solve();

        // Tidy up
        delete p_simulator;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTREPRESENTATIVESIMULATION_HPP_*/
