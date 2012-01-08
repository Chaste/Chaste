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

#ifndef TESTOFFLATTICESIMULATIONWITHVOLUMETRACKED_HPP_
#define TESTOFFLATTICESIMULATIONWITHVOLUMETRACKED_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"

#include "SmartPointers.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"

#include "OffLatticeSimulation.hpp"
#include "VolumeTrackedOffLatticeSimulation.hpp"

#include "ContactInhibitionCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"

#include "MutableMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "NodesOnlyMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "CellsGenerator.hpp"
#include "Warnings.hpp"

///\todo This test suite would be better named TestVolumeTrackedOffLatticeSimulation
class TestOffLatticeSimulationWithVolumeTracked : public AbstractCellBasedTestSuite
{
public:

    void TestNodeBasedSimulationWithContactInhibition()
    {
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetCellProliferativeType(STEM);
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.7);
            p_cycle_model->SetEquilibriumVolume(1.0);
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            cells.push_back(p_cell);
        }

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Create a contact inhibition simulator
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedSimulationWithVolumeTracked");
        simulator.SetEndTime(1.0);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellwiseData
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), p_data->GetValue(*cell_iter,0), 1e-4);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestMeshBasedSimulationWithContactInhibition()
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetCellProliferativeType(STEM);
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.7);
            p_cycle_model->SetEquilibriumVolume(1.0);
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            cells.push_back(p_cell);
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force law
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);

        // Initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Create a contact inhibition simulator
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMeshBasedSimulationWithVolumeTracked");
        simulator.SetEndTime(1.0);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellwiseData
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), p_data->GetValue(*cell_iter,0), 1e-4);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestVertexBasedSimulationWithContactInhibition()
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cell state
        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;

        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetCellProliferativeType(STEM);
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.7);
            p_cycle_model->SetEquilibriumVolume(1.0);
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);

        // Initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumElements(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Create a contact inhibition simulator
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedSimulationWithVolumeTracked");
        simulator.SetEndTime(1.0);
        simulator.AddForce(p_nagai_honda_force);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellwiseData
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), p_data->GetValue(*cell_iter,0), 1e-4);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestVolumeTrackedOffLatticeSimulationArchiving() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cell state
        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetCellProliferativeType(STEM);
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-1.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.7);
            p_cycle_model->SetEquilibriumVolume(1.0);
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force law
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);

        // Initialize CellwiseData
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Create a contact inhibition simulator
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVolumeTrackedOffLatticeSimulationSaveAndLoad");
        double end_time = 0.01;
        simulator.SetEndTime(end_time);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        CellBasedSimulationArchiver<2, VolumeTrackedOffLatticeSimulation<2> >::Save(&simulator);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumRealCells(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell->GetAge(), 1.01, 1e-4);

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Load simulation
        VolumeTrackedOffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, VolumeTrackedOffLatticeSimulation<2> >::Load("TestVolumeTrackedOffLatticeSimulationSaveAndLoad", end_time);

        p_simulator->SetEndTime(0.2);

        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->GetNumRealCells(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell2 = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell2->GetAge(), 1.01, 1e-4);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // Tidy up
        delete p_simulator;
        CellwiseData<2>::Destroy();

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::QuietDestroy();
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHVOLUMETRACKED_HPP_*/
