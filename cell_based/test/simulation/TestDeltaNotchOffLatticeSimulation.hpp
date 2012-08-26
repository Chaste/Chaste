/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTDELTANOTCHOFFLATTICESIMULATION_HPP_
#define TESTDELTANOTCHOFFLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <ctime>
#include "NodesOnlyMesh.hpp"
#include "DeltaNotchCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "DeltaNotchOffLatticeSimulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestDeltaNotchOffLatticeSimulation : public AbstractCellBasedTestSuite
{
private:
    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }
public:

    void TestUpdateAtEndOfTimeStepNodeBased() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a small population
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create some cells, each with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Initial condition for delta, notch, mean_delta
        std::vector<double> initial_conditions;
        initial_conditions.push_back(1.0);
        initial_conditions.push_back(1.0);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetInitialConditions(initial_conditions);
            p_model->SetDimension(2);
            p_model->SetMaxTransitGenerations(UINT_MAX);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);
        cell_population.SetCellAncestorsToLocationIndices();

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchNodeBasedUpdateAtEndOfTimeStep");
        simulator.SetEndTime(0.01);
        simulator.SetOutputNodeVelocities(true);

        // Set up force law and add to simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Check levels in cell 0
        double notch = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 1.0000, 1e-04);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestHeterogeneousDeltaNotchOnUntetheredTwoCellSystem()
    {
        EXIT_IF_PARALLEL;

        // Two cells close to each other
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.5, 0.0));

        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(nodes);

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

        // Establish a DNCCM for each of the cells
        DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
        p_model->SetDimension(2);
        p_model->SetInitialConditions(starter_conditions);

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(0.0);
        cells.push_back(p_cell);

        // Cell #2:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions_2;
        starter_conditions_2.push_back(0.19);
        starter_conditions_2.push_back(0.5);

        // Establish a DNCCM for each of the cells
        DeltaNotchCellCycleModel* p_model_2= new DeltaNotchCellCycleModel();
        p_model_2->SetDimension(2);
        p_model_2->SetInitialConditions(starter_conditions_2);

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        CellPtr p_cell_2(new Cell(p_state, p_model_2));
        p_cell_2->SetCellProliferativeType(p_diff_type);
        p_cell_2->SetBirthTime(0.0);
        cells.push_back(p_cell_2);

        // Create the cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up the simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchTwoCell_heterogee");
        simulator.SetEndTime(10.0);
        simulator.SetOutputNodeVelocities(true);

        // Define the radius of interaction as we're dealing with a node-based simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(0.75);
        simulator.AddForce(p_linear_force);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Run the simulation
        simulator.Solve();

        // Acquire cell pointers
        CellPtr p_cell_0b = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1b = cell_population.GetCellUsingLocationIndex(1);

        // Check that the simulation converges on the expected values
        double notch_0b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_0b->GetCellCycleModel())->GetNotch();
        TS_ASSERT_DELTA(notch_0b, 0.9640326, 1e-02);  //Default solution at t=10
        double delta_0b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_0b->GetCellCycleModel())->GetDelta();
        TS_ASSERT_DELTA(delta_0b, 0.0122205, 1e-04);  //Default solution at t=10
        double notch_1b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_1b->GetCellCycleModel())->GetNotch();
        TS_ASSERT_DELTA(notch_1b, 0.0261745, 1e-03);  //Default solution at t=10
        double delta_1b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_1b->GetCellCycleModel())->GetDelta();
        TS_ASSERT_DELTA(delta_1b, 0.8151536, 1e-02);  //Default solution at t=10

        // Avoid memory leaks
        delete p_mesh;
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

        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(nodes);

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

        // Establish a DNCCM for each of the cells
        DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
        p_model->SetDimension(2);
        p_model->SetInitialConditions(starter_conditions);

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(0.0);
        cells.push_back(p_cell);

        // Cell #2:
        //-----------------------------------------------------------------------------------
        // Create a vector of initial conditions
        std::vector<double> starter_conditions_2;
        starter_conditions_2.push_back(0.9);
        starter_conditions_2.push_back(0.5);

        // Establish a DNCCM for each of the cells
        DeltaNotchCellCycleModel* p_model_2 = new DeltaNotchCellCycleModel();
        p_model_2->SetDimension(2);
        p_model_2->SetInitialConditions(starter_conditions_2);

        // Ensure that all cells have birth time of zero in order to avoid
        // problems during the ODE solve of the ReadyToDivide() call, prior to
        // entering the main simulation timeloop
        CellPtr p_cell_2(new Cell(p_state, p_model_2));
        p_cell_2->SetCellProliferativeType(p_diff_type);
        p_cell_2->SetBirthTime(0.0);
        cells.push_back(p_cell_2);

        // Create the cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up the simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchTwoCell_homgee");
        simulator.SetEndTime(10.0);
        simulator.SetOutputNodeVelocities(true);

        // Define the radius of interaction as we're dealing with a node-based simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(0.75);
        simulator.AddForce(p_linear_force);

        cell_population.SetMechanicsCutOffLength(1.5);
        simulator.Solve();

        // Acquire cell pointers
        CellPtr p_cell_0b = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1b = cell_population.GetCellUsingLocationIndex(1);

        // Check that the simulation converges on the expected values
        double notch_0b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_0b->GetCellCycleModel())->GetNotch();
        TS_ASSERT_DELTA(notch_0b, 0.3538417, 1e-04);  //Default solution at t=10
        double delta_0b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_0b->GetCellCycleModel())->GetDelta();
        TS_ASSERT_DELTA(delta_0b, 0.0740040, 1e-04);  //Default solution at t=10
        double notch_1b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_1b->GetCellCycleModel())->GetNotch();
        TS_ASSERT_DELTA(notch_1b, 0.3538417, 1e-04);  //Default solution at t=10
        double delta_1b = dynamic_cast<DeltaNotchCellCycleModel*>(p_cell_1b->GetCellCycleModel())->GetDelta();
        TS_ASSERT_DELTA(delta_1b, 0.0740040, 1e-04);  //Default solution at t=10

        // Avoid memory leaks
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestUpdateAtEndOfTimeStepVertex() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a regular vertex mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchVertex2D");
        simulator.SetEndTime(0.01);
        simulator.SetOutputNodeVelocities(true);

        // Create force law and add to simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check levels in cell 0
        double notch = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 0.9921, 1e-04);
    }

    void TestUpdateAtEndOfTimeStepMeshBased() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a 2D honeycomb mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchMeshBasedUpdateAtEndOfTimeStep");
        simulator.SetEndTime(0.01);
        simulator.SetOutputNodeVelocities(true);

        // Set up force law and add to simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check levels in cell 0
        double notch = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetNotch();
        TS_ASSERT_DELTA(notch, 0.9999, 1e-04);
        double delta = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetDelta();
        TS_ASSERT_DELTA(delta, 0.9901, 1e-04);
        double mean_delta = dynamic_cast<DeltaNotchCellCycleModel*>(cells[0]->GetCellCycleModel())->GetMeanNeighbouringDelta();
        TS_ASSERT_DELTA(mean_delta, 1.0000, 1e-04);
    }

    void TestArchiving() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a regular vertex mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -1.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchOffLatticeSimulationSaveAndLoad");
        double end_time=0.01;
        simulator.SetEndTime(end_time);
        simulator.SetOutputNodeVelocities(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // Run and save simulation
        simulator.Solve();

        CellBasedSimulationArchiver<2, DeltaNotchOffLatticeSimulation<2> >::Save(&simulator);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), 16u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumElements(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell->GetAge(), 1.01, 1e-4);

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Load simulation
        DeltaNotchOffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, DeltaNotchOffLatticeSimulation<2> >::Load("TestDeltaNotchOffLatticeSimulationSaveAndLoad", end_time);

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
};

#endif /*TESTDELTANOTCHOFFLATTICESIMULATION_HPP_*/
