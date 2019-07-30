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

#ifndef TESTCRYPTSIMULATION2DNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "CryptSimulation2d.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CryptCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WntConcentration.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellBasedEventHandler.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SmartPointers.hpp"

// Cell population writers
#include "CellProliferativeTypesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCryptSimulation2dNightly : public AbstractCellBasedWithTimingsTestSuite
{
private:
    void setUp()
    {
        CellBasedEventHandler::Disable(); // these tests fail with event-handling on
        AbstractCellBasedWithTimingsTestSuite::setUp();
    }

public:
///////// NON-PERIODIC TESTS - these test the spring system and cell birth etc. /////////

    /**
     * Provides a reasonable test for the ghost node system...
     */
    void Test2DHoneycombMeshNotPeriodic()
    {
        // Create mesh
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        TS_ASSERT_EQUALS(cells.size(), location_indices.size());
        unsigned original_number_of_ghosts = p_mesh->GetNumNodes() - cells.size();
        TS_ASSERT_EQUALS(original_number_of_ghosts, 84u);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(12.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length, true, crypt_width));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DHoneycombMesh", false);

        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = p_mesh->GetNumNodes();

        /// \future This test is really fragile - gets different answers with Intel and Gcc compilers

        // 61 <= number_of_cells <= 63
        TS_ASSERT_LESS_THAN_EQUALS(number_of_cells, 63u);
        TS_ASSERT_LESS_THAN_EQUALS(61u, number_of_cells);
        // 145 <= number_of_nodes <= 147
        TS_ASSERT_LESS_THAN_EQUALS(number_of_nodes, 147u);
        TS_ASSERT_LESS_THAN_EQUALS(145u, number_of_nodes);

        std::set<unsigned> ghost_indices = crypt.GetGhostNodeIndices();
        TS_ASSERT_EQUALS(number_of_cells + ghost_indices.size(), number_of_nodes);

        std::vector<unsigned> cell_type_count = crypt.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 4u);

        TS_ASSERT_EQUALS(cell_type_count[0], 6u);   // Stem

        TS_ASSERT_LESS_THAN_EQUALS(cell_type_count[1], 23u); // 21 <= Transit <= 23
        TS_ASSERT_LESS_THAN_EQUALS(21u, cell_type_count[1]);

        TS_ASSERT_LESS_THAN_EQUALS(cell_type_count[2], 35u); // 33 <= Differentiated <= 35
        TS_ASSERT_LESS_THAN_EQUALS(33u, cell_type_count[2]);

        TS_ASSERT_EQUALS(cell_type_count[3], 0u);   // Default
    }

    void TestMonolayer()
    {
        // Create mesh
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        //double crypt_length = num_cells_depth - 1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, -1.0);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.AddPopulationWriter<VoronoiDataWriter>();
        crypt.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // Set the first cell to be logged
        crypt.Begin()->SetLogged();

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Monolayer");
        simulator.SetEndTime(1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = p_mesh->GetNumNodes();
        TS_ASSERT_EQUALS(number_of_cells, 69u);
        TS_ASSERT_EQUALS(number_of_nodes, 153u);

        std::set<unsigned> ghost_indices = crypt.GetGhostNodeIndices();
        TS_ASSERT_EQUALS(number_of_cells + ghost_indices.size(), number_of_nodes);

        std::vector<unsigned> cell_type_count = crypt.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 4u);
        TS_ASSERT_EQUALS(cell_type_count[0], 0u);   // Stem
        TS_ASSERT_EQUALS(cell_type_count[1], 33u);  // Transit
        TS_ASSERT_EQUALS(cell_type_count[2], 36u);  // Differentiated
        TS_ASSERT_EQUALS(cell_type_count[3], 0u);   // Default

        TS_ASSERT_EQUALS(crypt.GetVoronoiTessellation()->GetNumElements(), number_of_nodes);
        TS_ASSERT_EQUALS(crypt.GetVoronoiTessellation()->GetNumNodes(), 280u);
    }

    /**
     * Starting with a small mesh with one stem cell and the rest
     * differentiated, check that the number of cells at the end
     * of the simulation is as expected.
     */
    void Test2DCorrectCellNumbers()
    {
        // Create mesh
        unsigned num_cells_width = 7;

        unsigned num_cells_depth = 5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        double crypt_width = num_cells_width - 1.0;
        double crypt_length = num_cells_depth - 1.0;

        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = location_indices.size();
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            unsigned generation = 4;
            if (location_indices[i] == 27) // middle of bottom row of cells
            {
                generation = 0;
            }
            p_model->SetGeneration(generation);

            // Check the stem cell cycle time is still 24 hrs, otherwise this test might not pass
            //TS_ASSERT_DELTA(p_model->GetStemCellG1Duration(), 14, 1e-12); //These lines may trip up the Intel compiler with heavy optimization - don't know why?
            //TS_ASSERT_DELTA(p_model->GetTransitCellG1Duration(), 2, 1e-12);
            TS_ASSERT_DELTA(p_model->GetSG2MDuration(), 10, 1e-12);

            CellPtr p_cell(new Cell(p_state, p_model));

            p_cell->SetCellProliferativeType(p_diff_type);
            if (location_indices[i] == 27) // middle of bottom row of cells
            {
                p_cell->SetCellProliferativeType(p_stem_type);
            }

            p_cell->SetBirthTime(-1.0);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DSpringsCorrectCellNumbers");
        simulator.SetEndTime(40); // hours

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length, true, crypt_width));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Now count the number of each type of cell
        unsigned num_stem = 0;
        unsigned num_transit = 0;
        unsigned num_differentiated = 0;

        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            if (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
            {
                num_stem++;
            }
            else if (cell_iter->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
            {
                num_transit++;
            }
            else if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
            {
                num_differentiated++;
            }
            else
            {
                // shouldn't get here
                TS_ASSERT(false);
            }
        }

        TS_ASSERT_EQUALS(num_stem, 1u);
        TS_ASSERT_EQUALS(num_transit, 2u);

        TS_ASSERT_LESS_THAN(num_differentiated, 25u);
        TS_ASSERT_LESS_THAN(15u, num_differentiated);
    }

///////// PERIODIC TESTS - These test the system as a whole /////////

    void Test2DPeriodicNightly()
    {
        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        double crypt_length = cells_up*(sqrt(3.0)/2);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicNightly");
        simulator.SetEndTime(12.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_cells, 87u);
        TS_ASSERT_EQUALS(number_of_nodes, 135u);
    }

    void TestCrypt2DPeriodicWntNightly()
    {
        CellBasedEventHandler::Enable();

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        double crypt_length = cells_up*(sqrt(3.0)/2);

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");
        simulator.SetEndTime(24.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)

        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
#ifdef CHASTE_CVODE
        // divisions occur marginally earlier with CVODE
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 99u);
        TS_ASSERT_EQUALS(number_of_nodes, 147u);
#else
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 99u);
        TS_ASSERT_EQUALS(number_of_nodes, 147u);
#endif //CHASTE_CVODE

        // Tidy up
        WntConcentration<2>::Destroy();

        /*
         * HOW_TO_TAG Cell Based/Simulation
         * Time various aspects of a cell-based simulation using `CellBasedEventHandler`.

           Do not forget
            #include "CellBasedEventHandler.hpp"
           and to call
            CellBasedEventHandler::Enable();
           at the top of the test.

           If you are running multiple simulations then either
           place a CellBasedEventHandler::Reset(); before or after each solve
           or add the reset to the tearDown() method as above
         */

        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();
    }

    /**
     * This test is dontTest-ed out and not run every night as it
     * doesn't really test anything. It does show how to set up a
     * mutant simulation. Mutant viscosities are tested elsewhere
     * directly.
     */
    void dontRunTestWithMutantCellsUsingDifferentViscosities()
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        double crypt_length = cells_up*(sqrt(3.0)/2);

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            double dist_from_3_6 = sqrt((x-3)*(x-3)+(y-6)*(y-6));

            MAKE_PTR(WildTypeCellMutationState, p_healthy);
            MAKE_PTR(ApcTwoHitCellMutationState, p_apc2);
            if (dist_from_3_6 < 1.1)
            {
                cells[i]->SetMutationState(p_apc2);
            }
            else
            {
                cells[i]->SetMutationState(p_healthy);
            }
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicMutant");
        simulator.SetEndTime(12.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();

        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());

        unsigned number_of_cells = 0;
        unsigned number_of_mutant_cells = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            number_of_cells++;
            if (cell_iter->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
            {
                number_of_mutant_cells++;
            }
        }

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    void TestRandomDeathWithPeriodicMesh()
    {
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DRandomDeathPeriodic");
        simulator.SetEndTime(4.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&crypt,  0.700619609));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // There should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 1u);
    }
};

#endif /*TESTCRYPTSIMULATION2DNIGHTLY_HPP_*/
