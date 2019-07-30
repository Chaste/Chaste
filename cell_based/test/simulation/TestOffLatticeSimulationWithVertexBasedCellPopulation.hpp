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

#ifndef TESTOFFLATTICESIMULATIONWITHVERTEXBASEDCELLPOPULATION_HPP_
#define TESTOFFLATTICESIMULATIONWITHVERTEXBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "WelikyOsterForce.hpp"
#include "AbstractCellKiller.hpp"
#include "TargetedCellKiller.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexMeshWriter.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "FakePetscSetup.hpp"

class TestOffLatticeSimulationWithVertexBasedCellPopulation : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestSingleCellRelaxationNagaiHonda()
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, true, cos(theta), sin(theta)));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements(), std::vector<unsigned>(), p_diff_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSingleCellRelaxationNagaiHonda");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Test relaxes to circle (can be more stringent with more nodes and more time)
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetVolumeOfElement(0), 1.0, 0.05);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetSurfaceAreaOfElement(0), 3.5449077, 0.1);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    void TestDoNotRestrictVertexMovement()
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, true, cos(theta), sin(theta)));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements(), std::vector<unsigned>(), p_diff_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        TS_ASSERT(cell_population.GetRestrictVertexMovementBoolean())
        cell_population.SetRestrictVertexMovementBoolean(false);
        TS_ASSERT(!cell_population.GetRestrictVertexMovementBoolean())

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexSimulationWithoutMovementRestriction");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Test relaxes to circle (can be more stringent with more nodes and more time)
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetVolumeOfElement(0), 1.0, 0.05);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetSurfaceAreaOfElement(0), 3.549077, 0.1);

        // Test that the vertex restriction warning does not occur
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
    }

    void TestSingleCellRelaxationWelikyOster()
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, true, cos(theta), sin(theta)));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements(), std::vector<unsigned>(), p_diff_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSingleCellRelaxationWelikyOster");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(WelikyOsterForce<2>, p_weliky_oster_force);
        simulator.AddForce(p_weliky_oster_force);

        // Run simulation
        simulator.Solve();

        // Test relaxes to circle (can be more stringent with more nodes and more time)
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetVolumeOfElement(0), 1.0, 0.1);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetSurfaceAreaOfElement(0), 3.5449077, 0.1);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        Warnings::QuietDestroy();
    }

    void TestSingleCellDividing()
    {
        // Create a simple 2D MutableVertexMesh with only one cell
        HoneycombVertexMeshGenerator generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cell.
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_transit_type);
        double birth_time = -20.0; // divides straight away
        p_cell->SetBirthTime(birth_time);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSingleCellDividing");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Check that cell divided successfully
        unsigned new_num_nodes = simulator.rGetCellPopulation().GetNumNodes();
        unsigned new_num_elements = (static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumElements();
        unsigned new_num_cells = simulator.rGetCellPopulation().GetNumRealCells();

        TS_ASSERT_EQUALS(new_num_nodes, 8u);
        TS_ASSERT_EQUALS(new_num_elements, 2u);
        TS_ASSERT_EQUALS(new_num_cells, 2u);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    void TestVertexMonolayerWithCellBirth()
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(5, 5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            if (i==12)
            {
                cells[i]->SetCellProliferativeType(p_stem_type);
                cells[i]->SetBirthTime(-23.95);
            }
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        unsigned old_num_nodes = cell_population.GetNumNodes();
        unsigned old_num_elements = cell_population.GetNumElements();
        unsigned old_num_cells = cell_population.GetNumRealCells();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);

        // Test that the member mDt has been initialised correctly
        TS_ASSERT_DELTA(simulator.GetDt(), 0.002, 1e-6);

        simulator.SetOutputDirectory("TestVertexMonolayerWithCellBirth");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Check that cell 12 divided successfully
        unsigned new_num_nodes = simulator.rGetCellPopulation().GetNumNodes();
        unsigned new_num_elements = (static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumElements();
        unsigned new_num_cells = simulator.rGetCellPopulation().GetNumRealCells();

        TS_ASSERT_EQUALS(new_num_nodes, old_num_nodes+2); // as division of element is longer than threshold so is divided
        TS_ASSERT_EQUALS(new_num_elements, old_num_elements+1);
        TS_ASSERT_EQUALS(new_num_cells, old_num_cells+1);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    // This test uses a larger timestep to run faster.
    void TestVertexMonolayerWithVoid()
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create a horseshoe-shaped mesh
        p_mesh->DeleteElementPriorToReMesh(0);
        p_mesh->DeleteElementPriorToReMesh(1);
        p_mesh->DeleteElementPriorToReMesh(4);
        p_mesh->DeleteElementPriorToReMesh(6);

        p_mesh->ReMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-(double)i -19.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexMonolayerWithVoid");
        simulator.SetEndTime(20.0);
        simulator.SetSamplingTimestepMultiple(1);
        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier#2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        ////////////////////////////////////////////
        /// Strange setup to speed up simulation ///
        ////////////////////////////////////////////
        p_mesh->SetCellRearrangementThreshold(0.05);
        simulator.SetDt(0.1);
        ////////////////////////////////////////////

        // Run simulation
        simulator.Solve();

        // Check that void has been removed and vertices are in the correct position
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), 40u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumElements(),15u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 15u);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    void TestVertexMonolayerWithCellDeath()
    {
        /*
         * We don't want apoptosing cells to be labelled as dead after a certain time in
         * vertex simulations, so set the apoptosis time to something large.
         */

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(4,4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        /*
         * Set a high threshold for the maximum cell area for a T2 swap,
         * so that this occurs as soon a any element becomes triangular.
         */
        p_mesh->SetT2Threshold(1.0);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        cells[6]->StartApoptosis(false);
        cells[14]->StartApoptosis(true);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        unsigned old_num_nodes = cell_population.GetNumNodes();
        unsigned old_num_elements = cell_population.GetNumElements();
        unsigned old_num_cells = cell_population.GetNumRealCells();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexMonolayerWithCellDeath");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add cell killers to simulation (note we must account for element index changes following each kill)
        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_cell0_killer, (&cell_population, 0)); // element on the SW corner
        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_cell2_killer, (&cell_population, 2)); // element on the S boundary
        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_cell9_killer, (&cell_population, 9)); // element on the interior

        simulator.AddCellKiller(p_cell0_killer);
        simulator.AddCellKiller(p_cell2_killer);
        simulator.AddCellKiller(p_cell9_killer);

        // Run simulation
        simulator.Solve();

        // Check that cells 6 and 14 have now been removed
        unsigned new_num_nodes = simulator.rGetCellPopulation().GetNumNodes();
        unsigned new_num_elements = (static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumElements();
        unsigned new_num_cells = simulator.rGetCellPopulation().GetNumRealCells();

        TS_ASSERT_EQUALS(new_num_nodes, old_num_nodes-5);    // Due to the cells on the boundary that get killed and the apoptotic cell that does a T2 swap
        TS_ASSERT_EQUALS(new_num_elements, old_num_elements-4);
        TS_ASSERT_EQUALS(new_num_cells, old_num_cells-4);
        TS_ASSERT_EQUALS(new_num_cells, new_num_elements);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "A Cell is removed without performing a T2 swap. This could leave a void in the mesh.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    /*
     * This test visualizing cells of 2 mutation types, wildtype and labelled type.
     * It asserts that neighbouring cells have the correct adhesion parameter for difference
     * pairs of nodes.
     */
    void TestVertexMonolayerWithTwoMutationTypes()
    {
        // Create a simple 2D MutableVertexMesh with four cells
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_transit_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-2.0);
        }

        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        cells[0]->AddCellProperty(p_label);
        cells[2]->AddCellProperty(p_label);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexMonolayerWithTwoMutationStates");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        ///\todo test against a saved simulation or something similar, i.e check the positions of some vertices.
        TS_ASSERT_EQUALS(p_mesh->GetNode(13)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().front()->HasCellProperty<CellLabel>(), true);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    /*
     * This illustrates how the force class NagaiHondaDifferentialAdhesionForce
     * may be configured and used in a vertex-based simulation. See also the
     * tutorial TestRunningDifferentialAdhesionSimulationsTutorial.hpp.
     */
    void TestVertexMonolayerWithDifferentialAdhesion()
    {
        // Create a simple 2D MutableVertexMesh with four cells
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-2.0);
        }

        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        cells[0]->AddCellProperty(p_label);
        cells[2]->AddCellProperty(p_label);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexMonolayerWithDifferentialAdhesion");
        simulator.SetEndTime(1.0);

        // Set up force law and pass it to the simulation
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(6.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(3.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12.0);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(40.0);
        simulator.AddForce(p_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        ///\todo test against a saved simulation or something similar, i.e check the positions of some vertices.
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().front()->HasCellProperty<CellLabel>(), true);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    void TestSingleCellRelaxationAndApoptosis()
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, true, cos(theta), sin(theta)));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        MutableVertexMesh<2,2> mesh(nodes, elements);
        mesh.SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements(), std::vector<unsigned>(), p_diff_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexSingleCellApoptosis");
        simulator.SetEndTime(2.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Test relaxes to circle (can be more stringent with more nodes and more time)
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetVolumeOfElement(0), 1.0, 1e-1);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetSurfaceAreaOfElement(0), 3.5449077, 1e-1);

        CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(0);
        p_cell->StartApoptosis(false);

        simulator.SetEndTime(2.25); // any longer and cell target area is zero but element can't be removed as its the only one.

        // Run simulation
        simulator.Solve();

        TS_ASSERT_DELTA(cell_population.rGetMesh().GetVolumeOfElement(0), 0.5098, 1e-4);
        TS_ASSERT_DELTA(cell_population.rGetMesh().GetSurfaceAreaOfElement(0), 2.5417, 1e-3);

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    /*
     * This test is to stress test the vertex simulations by creating a massive monolayer,
     * it would also be useful for benchmarking.
     *
     * \todo In order to work the mesh archiving must support boundary nodes (see #1076)
     */
    void noTestVertexStressTest()
    {
        double start_time = 0.0;
        double end_time = 100.0;
        std::string output_directory = "StressTestVertex";

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(end_time);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Archive now and then reload

        // Save simulation in steady state
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Now save and reload to find where it breaks!
        for (unsigned i=0; i<40; i++)
        {
            start_time = end_time;
            end_time = end_time + 10.0;

            OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,start_time);
            p_simulator->SetDt(0.002);
            p_simulator->SetSamplingTimestepMultiple(50);
            p_simulator->SetEndTime(end_time);
            p_simulator->Solve();

            CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(p_simulator);
            delete p_simulator;
        }
        Warnings::QuietDestroy();
    }

    // Test archiving of a OffLatticeSimulation that uses a VertexBasedCellPopulation.
    void TestArchiving()
    {
        // Set end time
        double end_time = 0.1;

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(6, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithVertexBasedCellPopulationSaveAndLoad");
        simulator.SetEndTime(end_time);

        TS_ASSERT_DELTA(simulator.GetDt(), 0.002, 1e-12);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run and save simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 36u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), 96u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumElements(), 36u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-9);
        CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(23);
        TS_ASSERT_DELTA(p_cell->GetAge(), 23.1, 1e-4);

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Load simulation
        OffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestOffLatticeSimulationWithVertexBasedCellPopulationSaveAndLoad", end_time);

        p_simulator->SetEndTime(0.2);

        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 36u);
        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumNodes(), 96u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->GetNumElements(), 36u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-9);
        CellPtr p_cell2 = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(23);
        TS_ASSERT_DELTA(p_cell2->GetAge(), 23.1, 1e-4);

        // Test we can access the force collection, and that there is one force in it
        TS_ASSERT_EQUALS(p_simulator->rGetForceCollection().size(), 1u);

        // Get a pointer to the Nagai Honda Force, and verify we can access its methods
        boost::shared_ptr<NagaiHondaForce<2> > p_force =
                boost::static_pointer_cast<NagaiHondaForce<2> >(p_simulator->rGetForceCollection()[0]);

        double param_value = p_force->GetNagaiHondaDeformationEnergyParameter(); // Get the current value
        p_force->SetNagaiHondaDeformationEnergyParameter(1.23);
        TS_ASSERT_DELTA(p_force->GetNagaiHondaDeformationEnergyParameter(), 1.23, 1e-6);
        p_force->SetNagaiHondaDeformationEnergyParameter(param_value); // Set back to the original value


        // Run simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // Tidy up
        delete p_simulator;

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHVERTEXBASEDCELLPOPULATION_HPP_*/
