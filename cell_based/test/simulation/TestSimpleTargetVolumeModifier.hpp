/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTSIMPLETARGETVOLUMEMODIFIER_HPP_
#define TESTSIMPLETARGETVOLUMEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

#include "SimpleTargetVolumeModifier.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellBasedEventHandler.hpp"
#include "FileComparison.hpp"

// For the end-to-end growth simulation test
#include "HoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralMonolayerVertexMeshForce.hpp"

// This test is only run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSimpleTargetVolumeModifier : public AbstractCellBasedTestSuite
{
public:

    void TestNonPhaseBasedCellCycleModelMethodsAndExceptions()
    {
        // First set up SimulationTime (this is usually handled by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);

        // A non-phase-based cell-cycle model
        UniformCellCycleModel* p_model = new UniformCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_type);

        MAKE_PTR(SimpleTargetVolumeModifier<3>, p_modifier);

        // Without a growth duration, a non-phase-based cell-cycle model must throw
        TS_ASSERT_THROWS_THIS(p_modifier->UpdateTargetVolumeOfCell(p_cell),
            "If SetGrowthDuration() has not been called, a subclass of AbstractPhaseBasedCellCycleModel must be used");

        p_modifier->SetReferenceTargetVolume(9.0);
        p_modifier->SetGrowthDuration(2.0);

        SimulationTime::Instance()->IncrementTimeOneStep();
        p_modifier->UpdateTargetVolumeOfCell(p_cell);

        // At age 1, the cell should be halfway between its initial target volume 9.0/2 = 4.5 and 9.0
        TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("target volume"), 6.75, 1e-6);

        // Coverage of the accessors
        TS_ASSERT_DELTA(p_modifier->GetReferenceTargetVolume(), 9.0, 1e-6);
        TS_ASSERT_DELTA(p_modifier->GetGrowthDuration(), 2.0, 1e-6);

        CellBasedEventHandler::Reset(); // Logging was started but not stopped due to the exception above
    }

    void TestUpdateTargetVolumeOfCellBranches()
    {
        // SimulationTime is needed for cell ages and apoptosis timing
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(4.0, 4);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        MAKE_PTR(SimpleTargetVolumeModifier<3>, p_modifier);
        p_modifier->SetReferenceTargetVolume(2.0);
        p_modifier->SetGrowthDuration(4.0);

        // (i) A proliferating cell of age 2 (born at t = -2), half-way through the growth duration
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(-2.0);
            p_cell->InitialiseCellCycleModel();

            p_modifier->UpdateTargetVolumeOfCell(p_cell);
            // target = ref * 0.5 * (1 + age/growth_duration) = 2 * 0.5 * (1 + 2/4) = 1.5
            TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("target volume"), 1.5, 1e-9);
        }

        // (ii) A cell well past its cycle, so ReadyToDivide(): daughters restart at half the reference
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(-1000.0);
            p_cell->InitialiseCellCycleModel();
            TS_ASSERT(p_cell->ReadyToDivide());

            p_modifier->UpdateTargetVolumeOfCell(p_cell);
            TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("target volume"), 0.5 * 2.0, 1e-9);
        }

        // (iii) An apoptotic cell: its target volume decreases over time towards zero
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(-1.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->StartApoptosis();

            p_modifier->UpdateTargetVolumeOfCell(p_cell);
            const double apoptotic_target_1 = p_cell->GetCellData()->GetItem("target volume");
            TS_ASSERT_LESS_THAN(apoptotic_target_1, 2.0);
            TS_ASSERT_LESS_THAN_EQUALS(0.0, apoptotic_target_1);

            SimulationTime::Instance()->IncrementTimeOneStep();
            p_modifier->UpdateTargetVolumeOfCell(p_cell);
            const double apoptotic_target_2 = p_cell->GetCellData()->GetItem("target volume");
            TS_ASSERT_LESS_THAN(apoptotic_target_2, apoptotic_target_1);
        }
    }

    void TestSimpleTargetVolumeModifierArchiving()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SimpleTargetVolumeModifier.arch";

        // Separate scope to write the archive
        {
            AbstractCellBasedSimulationModifier<3,3>* const p_modifier = new SimpleTargetVolumeModifier<3>();
            (static_cast<SimpleTargetVolumeModifier<3>*>(p_modifier))->SetReferenceTargetVolume(14.3);
            (static_cast<SimpleTargetVolumeModifier<3>*>(p_modifier))->SetGrowthDuration(9.2);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_modifier;
            delete p_modifier;
        }

        // Separate scope to read the archive
        {
            AbstractCellBasedSimulationModifier<3,3>* p_modifier2;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_modifier2;

            double reference_target_volume = (static_cast<SimpleTargetVolumeModifier<3>*>(p_modifier2))->GetReferenceTargetVolume();
            TS_ASSERT_DELTA(reference_target_volume, 14.3, 1e-9);

            double growth_duration = (static_cast<SimpleTargetVolumeModifier<3>*>(p_modifier2))->GetGrowthDuration();
            TS_ASSERT_DELTA(growth_duration, 9.2, 1e-9);

            delete p_modifier2;
        }
    }

    void TestTargetVolumeGrowthDrivesCellGrowth()
    {
        // #480 End-to-end check of the growth loop: a SimpleTargetVolumeModifier grows each cell's
        // target volume over time, and the monolayer volume force makes the cells' actual volumes
        // follow, i.e. the cells grow (as in Okuda et al. 2013, where actual volume is emergent from
        // relaxation towards the growing target). Non-dividing (NoCellCycleModel) cells are used so
        // the growth is isolated from division, and only the volume term of the force is active.
        EXIT_IF_PARALLEL;

        const double z_height = 1.0;
        const double target_area = 1.0;
        const unsigned num_cells_x = 3;
        const unsigned num_cells_y = 3;

        // Build a small 3D monolayer by extruding a 2D honeycomb mesh
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder("TestTargetVolumeGrowth");
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Record the initial mean cell volume (before any growth)
        double initial_mean_volume = 0.0;
        for (unsigned i = 0; i < cell_population.GetNumElements(); ++i)
        {
            initial_mean_volume += cell_population.rGetMesh().GetVolumeOfElement(i);
        }
        initial_mean_volume /= cell_population.GetNumElements();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestTargetVolumeGrowthDrivesCellGrowth");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(0.5);

        // Only the volume-compressibility term is active; its global target (the fallback) is the
        // initial volume, so the growth comes entirely from the per-cell target set by the modifier.
        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force);
        p_force->SetVolumeParameters(50, z_height * target_area);
        simulator.AddForce(p_force);

        // Grow each cell's target volume from V/2 up to 2.0 over a growth duration of 1.0
        MAKE_PTR(SimpleTargetVolumeModifier<3>, p_growth_modifier);
        p_growth_modifier->SetReferenceTargetVolume(2.0);
        p_growth_modifier->SetGrowthDuration(1.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();

        // No division should have occurred
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);

        // The per-cell target volumes should have grown above the initial value of 1.0 (towards 2.0)
        double final_mean_target = 0.0;
        for (VertexBasedCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            TS_ASSERT(cell_iter->GetCellData()->HasItem("target volume"));
            final_mean_target += cell_iter->GetCellData()->GetItem("target volume");
        }
        final_mean_target /= cell_population.GetNumRealCells();
        TS_ASSERT_LESS_THAN(1.0, final_mean_target);

        // The actual cell volumes should have followed the growing target
        double final_mean_volume = 0.0;
        for (unsigned i = 0; i < cell_population.GetNumElements(); ++i)
        {
            final_mean_volume += cell_population.rGetMesh().GetVolumeOfElement(i);
        }
        final_mean_volume /= cell_population.GetNumElements();
        TS_ASSERT_LESS_THAN(initial_mean_volume, final_mean_volume);
    }

    void TestSimpleTargetVolumeModifierOutputParameters()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler output_file_handler("TestSimpleTargetVolumeModifierOutputParameters", false);

        MAKE_PTR(SimpleTargetVolumeModifier<3>, p_modifier);
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "SimpleTargetVolumeModifier-3");

        p_modifier->SetReferenceTargetVolume(6.2);
        p_modifier->SetGrowthDuration(1.7);

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("SimpleTargetVolumeModifier.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        // The file should contain both this class's and its parent's parameters
        TS_ASSERT(output_file_handler.FindFile("SimpleTargetVolumeModifier.parameters").Exists());
    }
};

#endif /*TESTSIMPLETARGETVOLUMEMODIFIER_HPP_*/
