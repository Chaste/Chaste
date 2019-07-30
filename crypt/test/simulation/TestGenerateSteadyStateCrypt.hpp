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

#ifndef TESTGENERATESTEADYSTATECRYPT_HPP_
#define TESTGENERATESTEADYSTATECRYPT_HPP_

#include <cxxtest/TestSuite.h>

#include <cstring>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "Version.hpp"
#include "FakePetscSetup.hpp"

class TestGenerateSteadyStateCrypt : public CxxTest::TestSuite
{
public:

    /*
     * This test can be used to generate steady state crypts for use
     * as the starting points of other simulations.
     *
     * You need to specify :
     * the kind of cell-cycle model to use on line 64,
     * WntConcentration on line 69,
     * change any model parameters around line 90,
     * and give the simulator options around line 95.
     */
    void TestGenerateSteadyStateCryptArchives()
    {
        std::string output_directory = "SteadyStateCrypt";

        double end_of_simulation = 150.0; // hours
        double time_of_each_run = 10.0; // for each run - the more saves and loads the better for testing this

        // Create mesh
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        double crypt_length = cells_up*(sqrt(3.0)/2.0)*crypt_width/cells_across;

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices, false, 30.0); // Last parameter adjusts Ghost spring stiffness in line with the linear_force later on

        // Set cell population to output cell types
        crypt.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetWntConcentrationParameter(1.0/3.0);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory(output_directory);

        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&simulator.rGetCellPopulation(), crypt_length));
        simulator.AddCellKiller(p_killer);

        // UNUSUAL SET UP HERE /////////////////////////////////////
        p_linear_force->SetMeinekeSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)

        simulator.UseJiggledBottomCells();
        // END OF UNUSUAL SET UP! //////////////////////////////////

        simulator.Solve();
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        /*
        * HOW_TO_TAG Cell Based/Simulation
        * Save and load ('checkpoint') a cell-based simulation to file.
        */
        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
        {
            CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("SteadyStateCrypt", t);
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();
            CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;
        }

        /*
         * This test is extremely unstable, and outcomes vary between compiler vendors, compiler versions, and
         * compiler optimizations.
         *
         * Previously, this test had the following expected_cell_count numbers:
         *      474u (GccOpt with CVODE on)
         *      450u (IntelProduction with CVODE on - only for Intel 2013, not 2015 or 2017)
         *      445u (GccOpt with CVODE off)
         *
         * More recently, it has been determined that expected_cell_count is:
         *      430u (Intel 2017 Release (O2 and O3)
         *      457u (Intel 2017 with O0)
         *
         * Given the fragility of this test, and that the numbers have been periodically updated when being built on
         * with different software versions, this test now simply checks that the resulting number of cells is in a
         * "sensible" range.  Further, we implicitly assume the test "passes" by not throwing at any point.
         *
         * \todo: Can we make this test more useful or robust?
         */
        {
            // Testing here that the number of cells at the end doesn't change
            CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("SteadyStateCrypt", end_of_simulation);

            // Define some "sensible" bounds on the number
            unsigned expected_cell_count_ub = 480u;
            unsigned expected_cell_count_lb = 425u;

            TS_ASSERT_LESS_THAN(p_simulator->rGetCellPopulation().GetNumRealCells(), expected_cell_count_ub);
            TS_ASSERT_LESS_THAN(expected_cell_count_lb, p_simulator->rGetCellPopulation().GetNumRealCells());

            delete p_simulator;
        }

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTGENERATESTEADYSTATECRYPT_HPP_*/
