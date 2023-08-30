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

#ifndef TESTEXTRINSICPULLMODIFIER_HPP_
#define TESTEXTRINSICPULLMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "OffLatticeSimulation.hpp"
#include "ExtrinsicPullModifier.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FarhadifarForce.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "Warnings.hpp"
#include "FileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestExtrinsicPullModifier : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedSimulationWithExtrinsicPull()
    {
        EXIT_IF_PARALLEL;    // Output in cell-based simulations doesn't work in parallel ///\todo #2356

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(2, 2);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedSimulationWithExtrinsicPull");
        simulator.SetDt(0.1);
        simulator.SetEndTime(2.0*simulator.GetDt());

        // Create an extrinsic pull modifier
        MAKE_PTR(ExtrinsicPullModifier<2>, p_modifier);

        // Test get methods
        p_modifier->SetApplyExtrinsicPullToAllNodes(false);
        TS_ASSERT_EQUALS(p_modifier->GetApplyExtrinsicPullToAllNodes(), false);
        TS_ASSERT_DELTA(p_modifier->GetSpeed(), 1.0, 1e-12);

        // Pass the modifier to the simulation
        simulator.AddSimulationModifier(p_modifier);

        // Run simulation
        simulator.Solve();

        // Test only the right-most nodes (indices 10 and 13) have moved
        unsigned node_idx = 0;
        for (unsigned i=0; i<2; i++)
        {
            TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(node_idx)->rGetLocation()[0], i + 0.5, 1e-6);
            node_idx++;
        }
        for (unsigned j=1; j<5; j++)
        {
            for (unsigned i=0; i<=2; i++)
            {
                double x_coord = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
                if ((node_idx == 10) || (node_idx == 13))
                {
                    TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(node_idx)->rGetLocation()[0], 2.7, 1e-6);
                }
                else
                {
                    TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(node_idx)->rGetLocation()[0], x_coord, 1e-6);
                }
                node_idx++;
            }
        }
        for (unsigned i=1; i<2; i++)
        {
            TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(node_idx)->rGetLocation()[0], i, 1e-6);
            node_idx++;
        }
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(node_idx)->rGetLocation()[0], 2, 1e-6);

        // Test set/get methods
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2> > >::iterator iter = simulator.GetSimulationModifiers()->begin();
        TS_ASSERT(boost::static_pointer_cast<ExtrinsicPullModifier<2> >(*iter));        
        boost::static_pointer_cast<ExtrinsicPullModifier<2> >(*iter)->SetApplyExtrinsicPullToAllNodes(true);
        boost::static_pointer_cast<ExtrinsicPullModifier<2> >(*iter)->SetSpeed(2.0);

        simulator.SetEndTime(4.0*simulator.GetDt());

        // Run simulation
        simulator.Solve();

        // Test that some of the nodes have the correct locations
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()[0], 0.5, 1e-3);
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(1)->rGetLocation()[0], 1.7222, 1e-3);
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(5)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(6)->rGetLocation()[0], 1.1481, 1e-3);
    }

    void TestSimulationArchivingWithExtrinsicPull()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(2, 2);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model and srn that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            NoCellCycleModel* p_cc_model = new NoCellCycleModel();
            p_cc_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_cc_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(-1.0);
            cells.push_back(p_cell);
        }
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSimulationArchivingWithExtrinsicPullSaveAndLoad");
        double end_time = 0.01;
        simulator.SetEndTime(end_time);

        // Create an extrinsic pull modifier and pass it to the simulation
        MAKE_PTR(ExtrinsicPullModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Run simulation
        simulator.Solve();

        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNumRealCells(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell->GetAge(), 1.01, 1e-4);

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Load simulation
        OffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestSimulationArchivingWithExtrinsicPullSaveAndLoad", end_time);

        p_simulator->SetEndTime(0.2);

        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 4u);
        TS_ASSERT_EQUALS((static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation())))->GetNumRealCells(), 4u);

        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.01, 1e-9);
        CellPtr p_cell2 = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(3);
        TS_ASSERT_DELTA(p_cell2->GetAge(), 1.01, 1e-4);

        // Run simulation
        p_simulator->Solve();

        // Tidy up
        delete p_simulator;

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        Warnings::QuietDestroy();
    }

    void TestExtrinsicPullModifierOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestExtrinsicPullModifierOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        MAKE_PTR(ExtrinsicPullModifier<2>, p_modifier);
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "ExtrinsicPullModifier-2");

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("ExtrinsicPullModifier.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("ExtrinsicPullModifier.parameters");
            FileFinder reference("cell_based/test/data/TestExtrinsicPullModifierOutputParameters/ExtrinsicPullModifier.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTEXTRINSICPULLMODIFIER_HPP_*/
