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

#ifndef TESTDELTANOTCHMODIFIER_HPP_
#define TESTDELTANOTCHMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "DeltaNotchSrnModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "DeltaNotchTrackingModifier.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellIdWriter.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestDeltaNotchModifier : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestUpdateAtEndOfTimeStepNodeBased()
    {
        EXIT_IF_PARALLEL;

        // Create a small 2D NodeBasedCellPopulation
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // ASsociate each cell with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Initial condition for delta, notch, mean_delta
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
            p_cc_model->SetDimension(2);

            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Create and configure cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchNodeBasedUpdateAtEndOfTimeStep");
        simulator.SetEndTime(0.01);

        // No mechanics so cells don't move

        // Create a Delta-Notch tracking modifier and add it to the simulation
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Run simulation
        simulator.Solve();

        // Test that the node velocities file exists
        OutputFileHandler output_file_handler("TestDeltaNotchNodeBasedUpdateAtEndOfTimeStep", false);
        FileFinder generated = output_file_handler.FindFile("results_from_time_0/nodevelocities.dat");
        TS_ASSERT(generated.Exists());

        // Check levels in cell 0
        CellPtr cell0 = cell_population.rGetCells().front();
        double notch = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 1.0000, 1e-04);
    }

    void TestHeterogeneousDeltaNotchOnUntetheredTwoCellSystem()
    {
        EXIT_IF_PARALLEL;

        // Two cells close to each other
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.5, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Establish a CCM for the cells and randomise the birth times
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Cell #1:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
        p_cc_model->SetDimension(2);

        DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
        p_srn_model->SetInitialConditions(starter_conditions);
        CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(0.0);
        cells.push_back(p_cell);

        // Cell #2:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions_2;
        starter_conditions_2.push_back(0.19);
        starter_conditions_2.push_back(0.5);

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model_2 = new UniformCellCycleModel();
        p_cc_model_2->SetDimension(2);

        DeltaNotchSrnModel* p_srn_model_2 = new DeltaNotchSrnModel();
        p_srn_model_2->SetInitialConditions(starter_conditions_2);
        CellPtr p_cell_2(new Cell(p_state, p_cc_model_2, p_srn_model_2));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell_2->SetCellProliferativeType(p_diff_type);
        p_cell_2->SetBirthTime(0.0);
        cells.push_back(p_cell_2);

        // Create the cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Set up the simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchTwoCell_heterogee");
        simulator.SetEndTime(10.0);

        // No mechanics so cells don't move

        // Add Delta-Notch tracking modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Run the simulation
        simulator.Solve();

        // Acquire cell pointers
        CellPtr p_cell_0b = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1b = cell_population.GetCellUsingLocationIndex(1);

        // Check that the simulation converges on the expected values
        double notch_0b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_0b->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch_0b, 0.9640326, 1e-02);  //Default solution at t=10
        double delta_0b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_0b->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta_0b, 0.0122205, 1e-04);  //Default solution at t=10
        double notch_1b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_1b->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch_1b, 0.0261745, 1e-03);  //Default solution at t=10
        double delta_1b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_1b->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta_1b, 0.8151536, 1e-02);  //Default solution at t=10


        // Now move cell so the cells have no neighboursthen they both run to a homogeneous steady state (note here mean delta=0)
        c_vector<double,2> old_point;
        old_point = static_cast<NodesOnlyMesh<2>* >(&(simulator.rGetCellPopulation().rGetMesh()))->GetNode(1)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+2.0;
        new_point.rGetLocation()[1] = old_point[1];
        static_cast<NodesOnlyMesh<2>* >(&(simulator.rGetCellPopulation().rGetMesh()))->SetNode(1, new_point);

        // Run the simulation
        simulator.SetEndTime(20.0);
        simulator.Solve();

        // Check that the simulation converges on the expected values
        notch_0b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_0b->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch_0b, 0.0, 1e-03);  //Default solution at t=20 for mean delta=0
        delta_0b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_0b->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta_0b, 1.0, 1e-03);  //Default solution at t=20 for mean delta=0
        notch_1b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_1b->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch_1b, 0.0, 1e-03);  //Default solution at t=20 for mean delta=0
        delta_1b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_1b->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta_1b, 1.0, 1e-03);  //Default solution at t=20 for mean delta=0


        // Avoid memory leaks
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestHomogeneousDeltaNotchOnUntetheredTwoCellSystem()
    {
        EXIT_IF_PARALLEL;

        // Two cells close to each other
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.5, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Establish a CCM for the cells and randomise the birth times
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Cell #1:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.9);
        starter_conditions.push_back(0.5);

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
        p_cc_model->SetDimension(2);

        DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
        p_srn_model->SetInitialConditions(starter_conditions);
        CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(0.0);
        cells.push_back(p_cell);

        // Cell #2:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions_2;
        starter_conditions_2.push_back(0.9);
        starter_conditions_2.push_back(0.5);

        // Establish a Stochastic CCM and a DN SRN for each of the cells
        UniformCellCycleModel* p_cc_model_2 = new UniformCellCycleModel();
        p_cc_model_2->SetDimension(2);

        DeltaNotchSrnModel* p_srn_model_2 = new DeltaNotchSrnModel();
        p_srn_model_2->SetInitialConditions(starter_conditions_2);
        CellPtr p_cell_2(new Cell(p_state, p_cc_model_2, p_srn_model_2));

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        p_cell_2->SetCellProliferativeType(p_diff_type);
        p_cell_2->SetBirthTime(0.0);
        cells.push_back(p_cell_2);


        // Create the cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Set up the simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchTwoCell_homgee");
        simulator.SetEndTime(10.0);

        // Add Delta-Notch tracking modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Define the radius of interaction as we're dealing with a node-based simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(0.75);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Acquire cell pointers
        CellPtr p_cell_0b = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1b = cell_population.GetCellUsingLocationIndex(1);

        // Check that the simulation converges on the expected values
        double notch_0b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_0b->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch_0b, 0.3538417, 1e-04);  //Default solution at t=10
        double delta_0b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_0b->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta_0b, 0.0740040, 1e-04);  //Default solution at t=10
        double notch_1b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_1b->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch_1b, 0.3538417, 1e-04);  //Default solution at t=10
        double delta_1b = dynamic_cast<DeltaNotchSrnModel*>(p_cell_1b->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta_1b, 0.0740040, 1e-04);  //Default solution at t=10

        // Avoid memory leaks
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestUpdateAtEndOfTimeStepVertex()
    {
        EXIT_IF_PARALLEL;

        // Create a regular vertex mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();


        // Initial condition for delta, notch
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        // Create some cells, each with a cell-cycle model and srn that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Create and configure cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchVertex2D");
        simulator.SetEndTime(0.01);

        // Add Delta-Notch tracking modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create force law and add to simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check levels in cell 0
        CellPtr cell0 = cell_population.rGetCells().front();
        double notch = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 0.9921, 1e-04);
    }

    void TestUpdateAtEndOfTimeStepMeshBasedWithGhostes()
    {
        EXIT_IF_PARALLEL;

        // Create a 2D honeycomb mesh
        HoneycombMeshGenerator generator(2, 2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();//**Changed**//

        // Initial condition for delta, notch
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        // Create some cells, each with a cell-cycle model and srn that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
            p_cc_model->SetDimension(2);

            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create and configure cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchMeshBasedWithGhostNodesUpdateAtEndOfTimeStep");
        simulator.SetEndTime(0.01);

        // Add Delta-Notch tracking modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Set up force law and add to simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
//        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        simulator.Solve();

        // Check levels in cell 0
        CellPtr cell0 = cell_population.rGetCells().front();
        double notch = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 1.0000, 1e-04);
    }

    void TestUpdateAtEndOfTimeStepPottsBased()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Initial condition for delta, notch
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        // Create some cells, each with a cell-cycle model and srn that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
            p_cc_model->SetDimension(2);

            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();

        // Create and configure cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchPottsBasedUpdateAtEndOfTimeStep");
        simulator.SetDt(0.01);
        simulator.SetEndTime(0.01);

        // Add Delta-Notch tracking modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Run simulation
        simulator.Solve();

        // Check levels in cell 0
        CellPtr cell0 = cell_population.rGetCells().front();
        double notch = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 1.0000, 1e-04);
    }

    void TestUpdateAtEndOfTimeStepCaBased()
    {
        EXIT_IF_PARALLEL;

        // Create cell population
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices;
        location_indices.push_back(0);
        location_indices.push_back(1);
        location_indices.push_back(4);
        location_indices.push_back(5);
        location_indices.push_back(6);
        location_indices.push_back(10);
        location_indices.push_back(11);
        location_indices.push_back(13);
        location_indices.push_back(14);
        location_indices.push_back(18);
        location_indices.push_back(19);

        // Initial condition for delta, notch
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        // Create some cells, each with a cell-cycle model and srn that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();
            p_cc_model->SetDimension(2);

            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellIdWriter>();

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchCaBasedUpdateAtEndOfTimeStep");
        simulator.SetDt(0.01);
        simulator.SetEndTime(0.01);

        // Add Delta-Notch tracking modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Run simulation
        simulator.Solve();

        // Check levels in cell 0 (this should be the same as for the Potts test, considering the configuration)
        CellPtr cell0 = cell_population.rGetCells().front();
        double notch = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchSrnModel*>(cell0->GetSrnModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 1.0000, 1e-04);
    }

    void TestArchiving()
    {
        EXIT_IF_PARALLEL;

        // Create a regular 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Initial condition for delta, notch
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        // Create some cells, each with a cell-cycle model and srn that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -1.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Create and configure cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchTrackingModifierSaveAndLoad");
        double end_time = 0.01;
        simulator.SetEndTime(end_time);

        // Add Delta-Notch tracking modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // NagaiHondaForce requires a growth modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run and save simulation
        simulator.Solve();

        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), 16u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumElements(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell->GetAge(), 1.01, 1e-4);

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Load simulation
        OffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestDeltaNotchTrackingModifierSaveAndLoad", end_time);

        p_simulator->SetEndTime(0.2);

        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumNodes(), 16u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->GetNumElements(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell2 = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell2->GetAge(), 1.01, 1e-4);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // Tidy up
        delete p_simulator;

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
        Warnings::QuietDestroy();
    }

    void TestDeltaNotchModifierOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestDeltaNotchModifierOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "DeltaNotchTrackingModifier-2");

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("DeltaNotchTrackingModifier.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("DeltaNotchTrackingModifier.parameters");
            FileFinder reference("cell_based/test/data/TestSimulationModifierOutputParameters/DeltaNotchTrackingModifier.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTDELTANOTCHMODIFIER_HPP_*/
