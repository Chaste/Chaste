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

#ifndef TESTGENERATESTEADYSTATECRYPT_HPP_
#define TESTGENERATESTEADYSTATECRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "SmartPointers.hpp"

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
    void TestGenerateSteadyStateCryptArchives() throw (Exception)
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
        crypt.SetOutputCellMutationStates(true);

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

        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
        {
            CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load("SteadyStateCrypt",t);
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();
            CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;
        }

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTGENERATESTEADYSTATECRYPT_HPP_*/
