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

#ifndef TESTCRYPTSIMULATION1D_HPP_
#define TESTCRYPTSIMULATION1D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CryptSimulation1d.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SloughingCellKiller.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "NumericFileComparison.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCryptSimulation1d : public AbstractCellBasedTestSuite
{
public:

    /**
     * In this test, there is no birth as only differentiated cells are used.
     * There is also no death because the we have not passed a cell killer
     * into the simulation. We just perturb one of the nodes and check that
     * each spring relaxes to its rest length.
     */
    void Test1dCryptWithNoBirthOrDeath()
    {
        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(21);

        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0] = 10.5;
        mesh.SetNode(10, shifted_point);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0 - node_index;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1dWithNoBirthAndNoDeath");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        // Note: converges very slowly so large tolerance of 0.1
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            Node<1>* p_node = mesh.GetNode(index);

            TS_ASSERT_EQUALS(p_node->IsDeleted(), false);

            double coord = p_node->rGetLocation()[0];
            TS_ASSERT_DELTA(coord, index, 0.05);
        }

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt1dWithNoBirthAndNoDeath", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        NumericFileComparison( node_results_file, "crypt/test/data/Crypt1dWithNoBirthAndNoDeath/results.viznodes").CompareFiles(0.0002 /*abstol*/);

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        FileComparison( cell_type_results_file, "crypt/test/data/Crypt1dWithNoBirthAndNoDeath/results.vizcelltypes").CompareFiles();

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        FileComparison( elem_results_file, "crypt/test/data/Crypt1dWithNoBirthAndNoDeath/results.vizelements").CompareFiles();
    }

    /**
     * In this test, we pass a sloughing cell killer into the simulation, and
     * check that a cell starting at the end of the crypt is sloughed off.
     */
    void Test1dCryptWithDeathButNoBirth()
    {
        double crypt_length = 22.0;

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(24);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0 - node_index;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1dWithDeathButNoBirth");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Add sloughing cell killer to simulation
        MAKE_PTR_ARGS(SloughingCellKiller<1>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 25u);

        // Run simulation
        simulator.Solve();

        // Since the default crypt length is 22, two cells should have been sloughed
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 23u);

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt1dWithDeathButNoBirth", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        FileComparison( node_results_file, "crypt/test/data/Crypt1dWithDeathButNoBirth/results.viznodes").CompareFiles();

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        FileComparison( cell_type_results_file, "crypt/test/data/Crypt1dWithDeathButNoBirth/results.vizcelltypes").CompareFiles();

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        FileComparison( elem_results_file, "crypt/test/data/Crypt1dWithDeathButNoBirth/results.vizelements").CompareFiles();
    }

    /**
     * In this test, we allow cells to proliferate.
     */
    void Test1dCryptWithBirthButNoDeath()
    {
        // Get pointers to singleton objects
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Set up cells by iterating through the nodes
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();

            double birth_time = 0;
            if (i == 0)
            {
                birth_time = -p_rand_gen->ranf()*(p_model->GetStemCellG1Duration() + p_model->GetSG2MDuration()); // hours - doesn't matter for stem cell
            }
            else if (i < 15)
            {
                birth_time = -p_rand_gen->ranf()*(p_model->GetTransitCellG1Duration() + p_model->GetSG2MDuration()); // hours
            }

            // Do this here due to fixing #1972 and not wanting to change saved results
            p_model->SetStemCellG1Duration(1.0);
            p_model->SetTransitCellG1Duration(1.0);

            // Set the cell-cycle model's generation
            unsigned generation = 4;
            if (i == 0)
            {
                generation = 0;
            }
            else if (i < 15)
            {
                generation = 1 + (i - 1) / 5;
            }
            p_model->SetGeneration(generation);

            CellPtr p_cell(new Cell(p_healthy_state, p_model));

            // Set the cell proliferative type
            if (i == 0)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
            }
            else if (i < 15)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff_type);
            }

            p_cell->InitialiseCellCycleModel();

            // Set the cell's birth time
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1dWithCells");
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        // There should be 34 cells at the end of the simulation
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 37u);

        // Work out where the previous test wrote its files
        OutputFileHandler handler("Crypt1dWithCells", false);

        std::string node_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        FileComparison( node_results_file, "crypt/test/data/Crypt1dWithCells/results.viznodes").CompareFiles();

        std::string cell_type_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        FileComparison( cell_type_results_file, "crypt/test/data/Crypt1dWithCells/results.vizcelltypes").CompareFiles();

        std::string elem_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        FileComparison( elem_results_file, "crypt/test/data/Crypt1dWithCells/results.vizelements").CompareFiles();
    }

    /**
     * In this test, we check that the daughters of a cell that has just divided
     * are put in the correct positions.
     */
    void TestCalculateCellDivisionVector()
    {
        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(3);

        // Set up cells by iterating through the nodes
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double birth_time;
            if (i == 1)
            {
                birth_time = -23.0;
            }
            else
            {
                birth_time = -1.0;
            }

            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("TestCalculateCellDivisionVector");
        simulator.SetEndTime(1.001);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        // There should be 5 cells at the end of the simulation
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 5u);

        // The second cell should have just divided, so be a distance 0.15 away from its initial location
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(1)->rGetLocation()[0], 1.15, 1e-3);

        // The new cell should also be a distance 0.15 away from the second cell's initial location
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(4)->rGetLocation()[0], 0.85, 1e-3);

        // The other cells should still be at their initial locations
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(2)->rGetLocation()[0], 2.0, 1e-3);
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(3)->rGetLocation()[0], 3.0, 1e-3);
    }

    /**
     * In this test, we include cell birth and cell death.
     */
    void Test1dCryptWithBirthAndDeath()
    {
        double crypt_length = 22.0;

        // Get pointers to singleton objects
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // Set up cells by iterating through the nodes
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();

            CellPtr p_cell(new Cell(p_healthy_state, p_model));

            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                generation = 0;
                birth_time = -p_rand_gen->ranf()*(p_model->GetStemCellG1Duration()
                                                  + p_model->GetSG2MDuration()); // hours - doesn't matter for stem cell
            }
            else if (i < 15)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 1 + (i - 1) / 5;
                birth_time = -p_rand_gen->ranf()*(p_model->GetTransitCellG1Duration()
                                                    + p_model->GetSG2MDuration()); // hours
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff_type);
                generation = 4;
                birth_time = 0; // hours
            }

            // Do this here due to fixing #1972 and not wanting to change saved results
            p_model->SetStemCellG1Duration(1.0);
            p_model->SetTransitCellG1Duration(1.0);
            p_model->SetGeneration(generation);

            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1dWithCellsAndGrowth");
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Add sloughing cell killer to simulation
        MAKE_PTR_ARGS(SloughingCellKiller<1>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // There should be 34 cells at the end of the simulation
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 31u);

        // Work out where the previous test wrote its files
        OutputFileHandler output_file_handler("Crypt1dWithCellsAndGrowth", false);

        std::string node_results_file = output_file_handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        FileComparison( node_results_file , "crypt/test/data/Crypt1dWithCellsAndGrowth/results.viznodes").CompareFiles();

        std::string cell_type_results_file = output_file_handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        FileComparison( cell_type_results_file, "crypt/test/data/Crypt1dWithCellsAndGrowth/results.vizcelltypes").CompareFiles();

        std::string elem_results_file = output_file_handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizelements";
        FileComparison( elem_results_file, "crypt/test/data/Crypt1dWithCellsAndGrowth/results.vizelements").CompareFiles();
    }

    void Test1DChainWithTysonNovakCellsAndNoDeath()
    {
        // Get pointers to singleton objects
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // Set up cells by iterating through the nodes
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        unsigned num_cells_at_start = mesh.GetNumNodes();
        std::vector<CellPtr> cells;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();

            CellPtr p_cell(new Cell(p_healthy_state, p_model));

            double birth_time;
            if (i == 0)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                /* Note we can't set the age of the cell to be too large or the ODE solver may
                 * not converge in the first timestep.
                 * Numerical value chosen so results are unchanged with old revisions.
                 * See #2788#comment:56
                 */
                birth_time = -p_rand_gen->ranf()*0.15; //(p_model->GetAverageStemCellCycleTime());
            }
            else if (i < 15)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                /* Note we can't set the age of the cell to be too large or the ODE solver may
                 * not converge in the first timestep.
                 * Numerical value chosen so results are unchanged with old revisions.
                 * See #2788#comment:56
                 */
                birth_time = -p_rand_gen->ranf()*0.15; //(p_model->GetAverageTransitCellCycleTime());
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff_type);
                birth_time = 0;
            }
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("CryptWithTysonNovakCells");
        simulator.SetEndTime(1.35);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);

        // Set the MeinekeSpringGrowthDuration to be the shorter than the default MPhase Duration as Tyson Novak sells divide every 1.25 hours.
        p_linear_force->SetMeinekeSpringGrowthDuration(0.1);

        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        /*
         * Check that several cell divisions have occurred. The Tyson-Novak cell
         * cycle has period 1.25, so by the end of the simulation the number of
         * cells should have doubled.
         *
         * Note that this test will fail if it is run for much longer, because the
         * Tyson-Novak cell-cycle period is so short that the cells because too
         * squashed together.
         */
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), num_cells_at_start + 23u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 2*num_cells_at_start);
    }

    /**
     * Create a crypt containing a single stem cell and all other cells differentiated.
     * Check that there is the correct number of cells at the end of a simulation, and
     * that they are in the correct order.
     */
    void Test1dChainCorrectCellNumbers()
    {
        double crypt_length = 5.0;

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(5);

        // For coverage, shift a node to a negative position
        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0] = -0.5;
        mesh.SetNode(0, shifted_point);

        // Set up cells by iterating through the nodes
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        unsigned num_cells = mesh.GetNumNodes();
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_model));

            unsigned generation = 4;
            if (i == 0)
            {
                generation = 0;
            }
            p_model->SetGeneration(generation);

            if (i == 0)
            {
                p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
            }
            else
            {
                p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
            }

            // The stem cell cycle time must still be 24 h, otherwise this test may not pass
            //TS_ASSERT_DELTA(p_model->GetStemCellG1Duration(), 14.0, 1e-12);//These lines may trip up the Intel compiler with heavy optimization - don't know why?
            //TS_ASSERT_DELTA(p_model->GetTransitCellG1Duration(), 2.0, 1e-12);
            TS_ASSERT_DELTA(p_model->GetSG2MDuration(), 10.0, 1e-12);

            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(0);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);
        crypt.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1dTestCorrectCellNumbers");
        simulator.SetEndTime(40);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Add sloughing cell killer to simulation
        MAKE_PTR_ARGS(SloughingCellKiller<1>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        unsigned num_cells_at_end = simulator.rGetCellPopulation().GetNumRealCells();

        /*
         * At the end of the simulation, there should be 6 live cells and one sloughed cell.
         * Note that throughout the simulation 2 cells are sloughed off, but one place is reused.
         */
        TS_ASSERT_EQUALS(num_cells_at_end, 6u);

        // Check we have the correct number of each cell type
        std::vector<unsigned> cell_type_count = simulator.rGetCellPopulation().GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 4u);
        TS_ASSERT_EQUALS(cell_type_count[0], 1u);
        TS_ASSERT_EQUALS(cell_type_count[1], 2u);
        TS_ASSERT_EQUALS(cell_type_count[2], 3u);
        TS_ASSERT_EQUALS(cell_type_count[3], 0u);

        for (AbstractCellPopulation<1>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
             cell_iter != simulator.rGetCellPopulation().End();
             ++cell_iter)
        {
            c_vector<double, 1> cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);
            double x = cell_location[0];

            if (fabs(x) < 1e-2)
            {
                TS_ASSERT_EQUALS(cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);
            }
            else if ((fabs(x-1) < 1e-2) || (fabs(x-2) < 1e-2))
            {
                TS_ASSERT_EQUALS(cell_iter->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);
            }
        }
    }

    /**
     * Test with Wnt-dependent cells.
     */
    void TestWntCellsCannotMoveAcrossYEqualsZero()
    {
        double crypt_length = 5.0;

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(3);

        // Set up cells by iterating through the nodes
        // (don't use any stem cells as we want to test the jiggling)
        unsigned num_cells = mesh.GetNumNodes();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        for (unsigned i=0; i<num_cells; i++)
        {
            WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();

            // The stem cell cycle time must still be 24 h, otherwise this test may not pass
            //TS_ASSERT_DELTA(p_cell_cycle_model1->GetStemCellG1Duration(), 14.0, 1e-12);
            //TS_ASSERT_DELTA(p_cell_cycle_model1->GetTransitCellG1Duration(), 2.0, 1e-12);
            TS_ASSERT_DELTA(p_cell_cycle_model1->GetSG2MDuration(), 10.0, 1e-12);

            p_cell_cycle_model1->SetDimension(1);
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_cycle_model1));
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            p_cell->SetBirthTime(0.0);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);
        crypt.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        crypt.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        AbstractCellPopulation<1>::Iterator cell_iterator = crypt.Begin();
        cell_iterator->SetBirthTime(-1.0);   // Make cell-cycle models do minimum work
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);

        boost::shared_ptr<AbstractCellProperty> p_apc1(crypt.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(crypt.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_label(crypt.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_iterator->AddCellProperty(p_label);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(p_apc1);
        ++cell_iterator;
        cell_iterator->SetBirthTime(-1.0);
        cell_iterator->SetMutationState(p_bcat1);

        // Create an instance of a Wnt concentration
        WntConcentration<1>::Instance()->SetType(LINEAR);
        WntConcentration<1>::Instance()->SetCellPopulation(crypt);
        WntConcentration<1>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1DWntMatureCells");
        simulator.SetEndTime(0.01);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        // Check that nothing has moved below y=0
        for (AbstractCellPopulation<1>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_LESS_THAN(-1e-15, crypt.GetLocationOfCellCentre(*cell_iter)[0]);
        }

        std::vector<unsigned> cell_mutation_state_count = simulator.rGetCellPopulation().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_state_count.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[0], 2u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_state_count[2], 0u); // No APC two hit
        TS_ASSERT_EQUALS(cell_mutation_state_count[3], 1u);

        std::vector<unsigned> cell_type_count = simulator.rGetCellPopulation().GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_type_count.size(), 4u);
        TS_ASSERT_EQUALS(cell_type_count[0], 0u);
        TS_ASSERT_EQUALS(cell_type_count[1], 4u);
        TS_ASSERT_EQUALS(cell_type_count[2], 0u);
        TS_ASSERT_EQUALS(cell_type_count[3], 0u);

        // Tidy up
        WntConcentration<1>::Destroy();
    }

    /**
     * Test saving a CryptSimulation1d object.
     */
    void TestSave()
    {
        double crypt_length = 22.0;

        // Get pointers to singleton object
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // Set up cells by iterating through the nodes
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_cell(new Cell(p_healthy_state, p_model));

            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                generation = 0;
                birth_time = -p_rand_gen->ranf()*(p_model->GetStemCellG1Duration()
                                                  + p_model->GetSG2MDuration());
            }
            else if (i < 15)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 1 + (i - 1)/5;
                birth_time = -p_rand_gen->ranf()*(p_model->GetTransitCellG1Duration()
                                                    + p_model->GetSG2MDuration());
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff_type);
                generation = 4;
                birth_time = 0;
            }
            p_model->SetGeneration(generation);
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1DSaveAndLoad");

        // Our full end time is 0.25, here we run until 0.1 then load and run more below.
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<1>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<1, CryptSimulation1d>::Save(&simulator);
    }

    /**
     * Test loading a CryptSimulation1d object.
     */
    void TestLoad()
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        CryptSimulation1d* p_simulator1;

        p_simulator1 = CellBasedSimulationArchiver<1, CryptSimulation1d>::Load("Crypt1DSaveAndLoad", 0.1);

        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();

        // Save then reload and run from 0.2 to 0.25
        MutableMesh<1,1>& r_mesh1 = (static_cast<MeshBasedCellPopulation<1>*>(&(p_simulator1->rGetCellPopulation())))->rGetMesh();
        CellBasedSimulationArchiver<1, CryptSimulation1d>::Save(p_simulator1);

        CryptSimulation1d* p_simulator2 = CellBasedSimulationArchiver<1, CryptSimulation1d>::Load("Crypt1DSaveAndLoad", 0.2);
        MutableMesh<1,1>& r_mesh2 = (static_cast<MeshBasedCellPopulation<1>*>(&(p_simulator2->rGetCellPopulation())))->rGetMesh();

        TS_ASSERT_EQUALS(r_mesh1.GetNumAllNodes(), r_mesh2.GetNumAllNodes());
        TS_ASSERT_EQUALS(r_mesh1.GetNumNodes(), r_mesh2.GetNumNodes());
        TS_ASSERT_EQUALS(r_mesh1.GetNumBoundaryNodes(), r_mesh2.GetNumBoundaryNodes());

        for (unsigned i=0; i<r_mesh1.GetNumAllNodes(); i++)
        {
            Node<1>* p_node = r_mesh1.GetNode(i);
            Node<1>* p_node2 = r_mesh2.GetNode(i);
            TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
            TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
            TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
            TS_ASSERT_DELTA(p_node->rGetLocation()[0], p_node2->rGetLocation()[0], 1e-16);
        }

        TS_ASSERT_EQUALS(r_mesh1.GetNumElements(), r_mesh2.GetNumElements());
        TS_ASSERT_EQUALS(r_mesh1.GetNumAllElements(), r_mesh2.GetNumAllElements());
        TS_ASSERT_EQUALS(r_mesh1.GetNumBoundaryElements(), r_mesh2.GetNumBoundaryElements());
        TS_ASSERT_EQUALS(r_mesh1.GetNumAllBoundaryElements(), r_mesh2.GetNumAllBoundaryElements());

        MutableMesh<1,1>::ElementIterator it = r_mesh1.GetElementIteratorBegin();
        MutableMesh<1,1>::ElementIterator it2 = r_mesh2.GetElementIteratorBegin();
        for (;
             it != r_mesh1.GetElementIteratorEnd();
             ++it, ++it2)
        {
            TS_ASSERT_EQUALS(it->GetNumNodes(), it2->GetNumNodes());
            for (unsigned i=0; i<it->GetNumNodes(); i++)
            {
                TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(i), it2->GetNodeGlobalIndex(i));
            }
        }

        p_simulator2->SetEndTime(0.25);

        // Run simulation
        p_simulator2->Solve();

        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 4.000, 1e-4);

        std::vector<double> node_11_location = p_simulator2->GetNodeLocation(11);
        TS_ASSERT_DELTA(node_11_location[0], 11.000, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }

    /**
     * In this test, we include cell birth and cell death.
     */
    void TestCryptSimulation1DParameterOutput()
    {
        double crypt_length = 22.0;

        // Get pointers to singleton objects
        RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

        // Create a mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(22);

        // Set up cells by iterating through the nodes
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_model));

            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                generation = 0;
                birth_time = -p_rand_gen->ranf()*(p_model->GetStemCellG1Duration()
                                                  + p_model->GetSG2MDuration()); // hours - doesn't matter for stem cell
            }
            else if (i < 15)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 1 + (i - 1) / 5;
                birth_time = -p_rand_gen->ranf()*(p_model->GetTransitCellG1Duration()
                                                    + p_model->GetSG2MDuration()); // hours
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff_type);
                generation = 4;
                birth_time = 0; // hours
            }

            p_model->SetGeneration(generation);
            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up crypt simulation
        CryptSimulation1d simulator(crypt);
        simulator.SetOutputDirectory("Crypt1dWithCellsAndGrowth");
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Add sloughing cell killer to simulation
        MAKE_PTR_ARGS(SloughingCellKiller<1>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Test Output methods

        ///\todo uncomment see #1453
        //TS_ASSERT_EQUALS(simulator.GetIdentifier(), "CryptSimulation1d");

        std::string output_directory = "TestCryptSimulation1dOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);
        out_stream parameter_file = output_file_handler.OpenOutputFile("crypt_sim_1d_results.parameters");
        simulator.OutputSimulationParameters(parameter_file);
        parameter_file->close();

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison( results_dir + "crypt_sim_1d_results.parameters", "crypt/test/data/TestCryptSimulationOutputParameters/crypt_sim_1d_results.parameters").CompareFiles();

        ///\todo check output of simulator.OutputSimulationSetup();
    }
};

#endif /*TESTCRYPTSIMULATION1D_HPP_*/
