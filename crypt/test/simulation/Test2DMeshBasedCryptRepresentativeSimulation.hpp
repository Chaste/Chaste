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

#ifndef TESTREPRESENTATIVESIMULATION_HPP_
#define TESTREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "SloughingCellKiller.hpp"

// Need to include all (or at least some) of these files!

#include "GeneralisedLinearSpringForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"

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

    void TestRepresentativeSimulationForProfiling() throw (Exception)
    {
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

        // Create a new clean directory
        OutputFileHandler file_handler(test_to_profile, true);

        // The archive must be copied from crypt/test/data/<test_to_profile>
        // to the testoutput directory to continue running the simulation
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string test_data_directory = "crypt/test/data/" + test_to_load +"/";
        std::string command = "cp -Rf --remove-destination " + test_data_directory +"* "+ test_output_directory +"/" + test_to_profile + "/";

        // Test that the above command was implemented successfully
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

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
