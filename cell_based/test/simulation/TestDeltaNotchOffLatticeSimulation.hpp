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
#include "CellwiseData.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
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

    void TestPostSolveNodeBased() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a small population
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_mesh);

        // Create some cells, each with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2);
            p_model->SetMaxTransitGenerations(UINT_MAX);
            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);
        cell_population.SetCellAncestorsToLocationIndices();

        // Create and initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumNodes(), 3);
        p_data->SetCellPopulation(&cell_population);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchNodeBasedPostSolve");
        simulator.SetEndTime(0.01);

        // Set up force law and add to simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check levels in cell 0 ///\todo see #1995
//        TS_ASSERT_DELTA(p_data->GetValue(cells[0],0),0.9384,1e-4);
//        TS_ASSERT_DELTA(p_data->GetValue(cells[1],0),0.9990,1e-4);
//        TS_ASSERT_DELTA(p_data->GetValue(cells[2],0),0.9588,1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestPostSolveVertex() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a regular vertex mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2u);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumElements(), 3);
        p_data->SetCellPopulation(&cell_population);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchVertex2D");
        simulator.SetEndTime(0.01);

        // Create force law and add to simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
 ///\todo see #1995
//        TS_ASSERT_DELTA(p_data->GetValue(cells[0],0),0.9386,5e-4);
//                TS_ASSERT_DELTA(p_data->GetValue(cells[1],0),0.9990,5e-4);
//        TS_ASSERT_DELTA(p_data->GetValue(cells[2],0),0.9589,5e-4);

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestPostSolveMeshBased() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a 2D honeycomb mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumNodes(), 3);
        p_data->SetCellPopulation(&cell_population);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchMeshBasedPostSolve");
        simulator.SetEndTime(0.01);

        // Set up force law and add to simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check levels in cell 0 ///\todo see #1995
//        TS_ASSERT_DELTA(p_data->GetValue(cells[0],0),0.9384,1e-4);
//        TS_ASSERT_DELTA(p_data->GetValue(cells[1],0),0.9990,1e-4);
//        TS_ASSERT_DELTA(p_data->GetValue(cells[2],0),0.9588,1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();
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
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(2u);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -1.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(p_mesh->GetNumElements(), 3);
        p_data->SetCellPopulation(&cell_population);

        // Create and configure cell-based simulation
        DeltaNotchOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchOffLatticeSimulationSaveAndLoad");
        double end_time=0.01;
        simulator.SetEndTime(end_time);

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
