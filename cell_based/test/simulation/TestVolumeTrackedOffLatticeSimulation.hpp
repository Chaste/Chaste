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

#ifndef TESTVOLUMETRACKEDOFFLATTICESIMULATION_HPP_
#define TESTVOLUMETRACKEDOFFLATTICESIMULATION_HPP_

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

class TestVolumeTrackedOffLatticeSimulation : public AbstractCellBasedTestSuite
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
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellwiseData
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), p_data->GetValue(*cell_iter,0), 1e-4);
        }

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
        simulator.AddForce(p_force);
        simulator.SetEndTime(simulator.GetDt()/2.0);

         // Run simulation
         simulator.Solve();

         // Test that the volumes of the cells are correct in CellwiseData
         for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
              cell_iter != cell_population.End();
              ++cell_iter)
         {
             TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), p_data->GetValue(*cell_iter,0), 1e-4);
         }

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
        simulator.AddForce(p_nagai_honda_force);

        simulator.SetEndTime(simulator.GetDt()/2.0);

         // Run simulation
         simulator.Solve();

         // Test that the volumes of the cells are correct in CellwiseData
         for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
              cell_iter != cell_population.End();
              ++cell_iter)
         {
             TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), p_data->GetValue(*cell_iter,0), 1e-4);
         }

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

#endif /*TESTVOLUMETRACKEDOFFLATTICESIMULATION_HPP_*/
