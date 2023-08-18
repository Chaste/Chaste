/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTDIVISIONBIASTRACKINGMODIFIER_HPP_
#define TESTDIVISIONBIASTRACKINGMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "BiasedBernoulliTrialCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DivisionBiasTrackingModifier.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "Warnings.hpp"
#include "FileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestDivisionBiasTrackingModifier : public AbstractCellBasedTestSuite
{
public:

    void TestNodeBasedSimulationWithDivisionBias()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D NodeBasedCellPopulation
        HoneycombMeshGenerator generator(5, 5, 0);
        boost::shared_ptr<TetrahedralMesh<2,2> > p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            BiasedBernoulliTrialCellCycleModel* p_cycle_model = new BiasedBernoulliTrialCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetMaxDivisionProbability(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedSimulationWithDivisionBias");

        TS_ASSERT_EQUALS(simulator.GetDt(), 1.0/120.0); // Default value for off-lattice simulations
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a division bias tracking modifier and pass it to the simulation
        c_vector<double, 2> bias_vector;
        bias_vector(0) = 0.0;
        bias_vector(1) = 1.0;
        MAKE_PTR_ARGS(DivisionBiasTrackingModifier<2>, p_modifier, (bias_vector));
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test that the cell data are correct at the first timestep
        double y_min = cell_population.rGetMesh().CalculateBoundingBox().rGetLowerCorner()[1];
        double y_max = cell_population.rGetMesh().CalculateBoundingBox().rGetUpperCorner()[1];
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            double estimated_bias = (cell_population.GetLocationOfCellCentre(*cell_iter)[1] - y_min)/(y_max - y_min);
            
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
            TS_ASSERT_DELTA(estimated_bias, bias, 1e-3);
        }
    
        simulator.SetEndTime(2.0);
        simulator.Solve();

        // Test that the cell data are correct at the end time
        y_min = cell_population.rGetMesh().CalculateBoundingBox().rGetLowerCorner()[1];
        y_max = cell_population.rGetMesh().CalculateBoundingBox().rGetUpperCorner()[1];
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
          cell_iter != cell_population.End();
          ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            double estimated_bias = (cell_population.GetLocationOfCellCentre(*cell_iter)[1] - y_min)/(y_max - y_min);
            
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
            TS_ASSERT_DELTA(estimated_bias, bias, 1e-3);
        }
    }

    void TestMeshBasedSimulationWithDivisionBias()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        HoneycombMeshGenerator generator(3, 3);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            BiasedBernoulliTrialCellCycleModel* p_cycle_model = new BiasedBernoulliTrialCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetMaxDivisionProbability(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMeshBasedSimulationWithDivisionBias");
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a division bias tracking modifier and pass it to the simulation
        c_vector<double, 2> bias_vector;
        bias_vector(0) = 1.0;
        bias_vector(1) = 0.0;
        MAKE_PTR_ARGS(DivisionBiasTrackingModifier<2>, p_modifier, (bias_vector));
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test that the cell data are correct at the first timestep
        double x_min = cell_population.rGetMesh().CalculateBoundingBox().rGetLowerCorner()[0];
        double x_max = cell_population.rGetMesh().CalculateBoundingBox().rGetUpperCorner()[0];
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            double estimated_bias = (cell_population.GetLocationOfCellCentre(*cell_iter)[0] - x_min)/(x_max - x_min);
            
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
            TS_ASSERT_DELTA(estimated_bias, bias, 1e-3);
        }
        
        simulator.SetEndTime(2.0);
        simulator.Solve();

        // Test that the cell data are correct at the end time
        x_min = cell_population.rGetMesh().CalculateBoundingBox().rGetLowerCorner()[0];
        x_max = cell_population.rGetMesh().CalculateBoundingBox().rGetUpperCorner()[0];
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            double estimated_bias = (cell_population.GetLocationOfCellCentre(*cell_iter)[0] - x_min)/(x_max - x_min);
            
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
            TS_ASSERT_DELTA(estimated_bias, bias, 1e-3);
        }
    }

    void TestMeshBasedSimulationWithGhostNodesAndDivisionBias()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel.

        // Create a simple 2D MeshBasedCellPopulationWithGhostNodes
        HoneycombMeshGenerator generator(3, 3, 3);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            BiasedBernoulliTrialCellCycleModel* p_cycle_model = new BiasedBernoulliTrialCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetMaxDivisionProbability(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells,location_indices);

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMeshBasedSimulationWithGhostNodesAndVolumeTracked");
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a division bias tracking modifier and pass it to the simulation
        c_vector<double, 2> bias_vector;
        bias_vector(0) = 0.0;
        bias_vector(1) = 1.0;
        MAKE_PTR_ARGS(DivisionBiasTrackingModifier<2>, p_modifier, (bias_vector));
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test that the cell data are in the correct interval in CellData at the first timestep
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
        }
        
        simulator.SetEndTime(2.0);
        simulator.Solve();

        // Test that the cell data are in the correct interval in CellData at the end time
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
        }
    }

    void TestVertexBasedSimulationWithDivisionBias()
    {
        EXIT_IF_PARALLEL;    // Output in cell-based simulations doesn't work in parallel ///\todo #2356

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(2, 2);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;

        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            BiasedBernoulliTrialCellCycleModel* p_cycle_model = new BiasedBernoulliTrialCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetMaxDivisionProbability(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedSimulationWithDivisionBias");
        simulator.SetEndTime(simulator.GetDt()/2.0);

        // Create a division bias tracking modifier and pass it to the simulation
        c_vector<double, 2> bias_vector;
        bias_vector(0) = 0.0;
        bias_vector(1) = 1.0;
        MAKE_PTR_ARGS(DivisionBiasTrackingModifier<2>, p_modifier, (bias_vector));
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);

        // Pass a target area modifier to the simulation
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(1.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        // Test that the cell data are correct at the first timestep
        double y_min = 100;
        double y_max = -100;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            if (y >y_max)
            {
                y_max = y;
            }
            if (y < y_min)
            {
                y_min = y;
            }
        }
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            double estimated_bias = (cell_population.GetLocationOfCellCentre(*cell_iter)[1] - y_min)/(y_max - y_min);
            
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
            TS_ASSERT_DELTA(estimated_bias, bias, 1e-3);
        }
    
        simulator.SetEndTime(2.0);
        simulator.Solve();

        // Test that the cell data are correct at the end time
        y_min = 100;
        y_max = -100;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
            if (y >y_max)
            {
                y_max = y;
            }
            if (y < y_min)
            {
                y_min = y;
            }
        }
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double bias = cell_iter->GetCellData()->GetItem("bias");
            double estimated_bias = (cell_population.GetLocationOfCellCentre(*cell_iter)[1] - y_min)/(y_max - y_min);
            
            TS_ASSERT_LESS_THAN_EQUALS(0.0, bias);
            TS_ASSERT_LESS_THAN_EQUALS(bias, 1.0);
            TS_ASSERT_DELTA(estimated_bias, bias, 1e-3);
        }
    }

    void TestDivisionBiasOffLatticeSimulationArchiving()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D MeshBasedCellPopulation
        HoneycombMeshGenerator generator(2, 2, 0);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            BiasedBernoulliTrialCellCycleModel* p_cycle_model = new BiasedBernoulliTrialCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-10.0);
            p_cycle_model->SetMaxDivisionProbability(0.1);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDivisionBiasTrackedOffLatticeSimulationSaveAndLoad");
        double end_time = 0.01;
        simulator.SetEndTime(end_time);

        // Create a division bias tracking modifier and pass it to the simulation
        c_vector<double, 2> bias_vector;
        bias_vector(0) = 1.0;
        bias_vector(1) = 0.0;
        MAKE_PTR_ARGS(DivisionBiasTrackingModifier<2>, p_modifier, (bias_vector));
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
        TS_ASSERT_DELTA(p_cell->GetAge(), 10.01, 1e-4);

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Load simulation
        OffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestDivisionBiasTrackedOffLatticeSimulationSaveAndLoad", end_time);

        p_simulator->SetEndTime(0.2);

        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->GetNumRealCells(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell2 = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell2->GetAge(), 10.01, 1e-4);

        // Run simulation
        p_simulator->Solve();

        // Tidy up
        delete p_simulator;

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        Warnings::QuietDestroy();
    }

    void TestDivisionBiasTrackingModifierOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestDivisionBiasTrackingModifierOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        c_vector<double, 2> bias_vector;
        bias_vector(0) = 1.0;
        bias_vector(1) = 0.0;
        MAKE_PTR_ARGS(DivisionBiasTrackingModifier<2>, p_modifier, (bias_vector));
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "DivisionBiasTrackingModifier-2");

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("DivisionBiasTrackingModifier.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("DivisionBiasTrackingModifier.parameters");
            FileFinder reference("cell_based/test/data/TestSimulationModifierOutputParameters/DivisionBiasTrackingModifier.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTDIVISIONBIASTRACKINGMODIFIER_HPP_*/
