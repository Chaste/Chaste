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

#ifndef TESTVOLUMETRACKINGMODIFIER_HPP_
#define TESTVOLUMETRACKINGMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "OffLatticeSimulation.hpp"
#include "VolumeTrackingModifier.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellsGenerator.hpp"
#include "Warnings.hpp"
#include "CellVolumesWriter.hpp"
#include "FileComparison.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestVolumeTrackingModifier : public AbstractCellBasedTestSuite
{
public:

    void TestNodeBasedSimulationWithContactInhibition()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D NodeBasedCellPopulation
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetEquilibriumVolume(0.7854); //pi *(0.5)^2
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedSimulationWithVolumeTracked");

        TS_ASSERT_EQUALS(simulator.GetDt(), 1.0/120.0); // Default value for off-lattice simulations
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test that the node velocities file exists
        OutputFileHandler output_file_handler("TestNodeBasedSimulationWithVolumeTracked", false);
        FileFinder generated = output_file_handler.FindFile("results_from_time_0/nodevelocities.dat");
        TS_ASSERT(generated.Exists());

        // Test that the volumes of the cells are correct in CellData at the first timestep
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), cell_iter->GetCellData()->GetItem("volume"), 1e-4);
        }

        simulator.SetEndTime(2.0);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellData at the end time
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), cell_iter->GetCellData()->GetItem("volume"), 1e-4);
        }

        // Check that the correct number of cells are labelled (i.e. experiencing contact inhibition)
        TS_ASSERT_EQUALS(cell_population.GetCellPropertyRegistry()->Get<CellLabel>()->GetCellCount(), 1u);
    }

    void TestMeshBasedSimulationWithContactInhibition()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        HoneycombMeshGenerator generator(3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetEquilibriumVolume(0.866); //sqrt(3.0)/2
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMeshBasedSimulationWithVolumeTracked");
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellData at the first timestep
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), (*cell_iter)->GetCellData()->GetItem("volume"), 1e-4);
        }

        simulator.SetEndTime(2.0);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellData at the end time
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), (*cell_iter)->GetCellData()->GetItem("volume"), 1e-4);
        }

        // Check that the correct number of cells are labelled (i.e. experiencing contact inhibition)
        TS_ASSERT_EQUALS(cell_population.GetCellPropertyRegistry()->Get<CellLabel>()->GetCellCount(),9u);
    }

    void TestMeshBasedSimulationWithGhostNodesAndContactInhibition()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel.

        // Create a simple 2D MeshBasedCellPopulationWithGhostNodes
        HoneycombMeshGenerator generator(3, 3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetEquilibriumVolume(0.866); //sqrt(3.0)/2
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells,location_indices);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMeshBasedSimulationWithGhostNodesAndVolumeTracked");
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellData at the first timestep
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), (*cell_iter)->GetCellData()->GetItem("volume"), 1e-4);
        }

        simulator.SetEndTime(2.0);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellData at the end time
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), (*cell_iter)->GetCellData()->GetItem("volume"), 1e-4);
        }

        // Check that the correct number of cells are labelled (i.e. experiencing contact inhibition)
        TS_ASSERT_EQUALS(cell_population.GetCellPropertyRegistry()->Get<CellLabel>()->GetCellCount(), 12u);
    }

    void TestVertexBasedSimulationWithContactInhibition()
    {
        EXIT_IF_PARALLEL;    // Output in cell-based simulations doesn't work in parallel ///\todo #2356

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;

        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.99); // Very high as cells are not that compressed, as low in number and want to check that contact inhibition works.
            p_cycle_model->SetEquilibriumVolume(1.0); // Target volume in NagaiHonda force
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedSimulationWithVolumeTracked");
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellData at the first timestep
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
          cell_iter != cell_population.End();
          ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), (*cell_iter)->GetCellData()->GetItem("volume"), 1e-4);
        }

        simulator.SetEndTime(2.0);

        // Run simulation
        simulator.Solve();

        // Test that the volumes of the cells are correct in CellData at the end time
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
          cell_iter != cell_population.End();
          ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), (*cell_iter)->GetCellData()->GetItem("volume"), 1e-4);
        }

        // Check that the correct number of cells are labelled (i.e. experiencing contact inhibition)
        TS_ASSERT_EQUALS(cell_population.GetCellPropertyRegistry()->Get<CellLabel>()->GetCellCount(),8u);
    }

    void TestVolumeTrackedOffLatticeSimulationArchiving()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D MeshBasedCellPopulation
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-1.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.7);
            p_cycle_model->SetEquilibriumVolume(1.0);
            p_cycle_model->SetStemCellG1Duration(0.1);
            p_cycle_model->SetTransitCellG1Duration(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->InitialiseCellCycleModel();

            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a contact inhibition simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVolumeTrackedOffLatticeSimulationSaveAndLoad");
        double end_time = 0.01;
        simulator.SetEndTime(end_time);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumRealCells(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell->GetAge(), 1.01, 1e-4);

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Load simulation
        OffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestVolumeTrackedOffLatticeSimulationSaveAndLoad", end_time);

        p_simulator->SetEndTime(0.2);

        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->GetNumRealCells(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell2 = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell2->GetAge(), 1.01, 1e-4);

        // Run simulation
        p_simulator->Solve();

        // Tidy up
        delete p_simulator;

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        Warnings::QuietDestroy();
    }

    void TestVolumeTrackingModifierOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestVolumeTrackingModifierOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "VolumeTrackingModifier-2");

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("VolumeTrackingModifier.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("VolumeTrackingModifier.parameters");
            FileFinder reference("cell_based/test/data/TestSimulationModifierOutputParameters/VolumeTrackingModifier.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTVOLUMETRACKINGMODIFIER_HPP_*/
