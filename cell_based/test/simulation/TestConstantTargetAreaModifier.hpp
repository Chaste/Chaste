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

#ifndef TESTCONSTANTTARGETAREAMODIFIER_HPP_
#define TESTCONSTANTTARGETAREAMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

#include "ConstantTargetAreaModifier.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "NagaiHondaForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "FileComparison.hpp"

// This test is only run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestConstantTargetAreaModifier : public AbstractCellBasedTestSuite
{
public:

    void TestTargetAreaOfDaughterCellsConstant()
    {
        // Create a simple 2D MutableVertexMesh with only one cell
        HoneycombVertexMeshGenerator generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_transit_type);
        double birth_time = -11.0; // The cell cycle duration is 12
        p_cell->SetBirthTime(birth_time);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestTargetAreaOfDaughterCellsSimple");
        simulator.SetEndTime(0.997);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // Create a ConstantTargetAreaModifier
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        p_growth_modifier->SetReferenceTargetArea(2.3);

        // Run simulation
        simulator.Solve();

        // We should only have one cell now
        unsigned num_cells_before_division = simulator.rGetCellPopulation().GetNumRealCells();
        TS_ASSERT_EQUALS(num_cells_before_division, 1u);

        // This is the cell from before; let's see what its target area is
        double target_area_before_division = p_cell->GetCellData()->GetItem("target area");
        TS_ASSERT_DELTA(target_area_before_division, 2.3, 1e-9);

        // We now adjust the end time and run the simulation a bit further
        simulator.SetEndTime(1.001);
        simulator.Solve();

        // We should now have two cells
        unsigned num_cells_at_division = simulator.rGetCellPopulation().GetNumRealCells();
        TS_ASSERT_EQUALS(num_cells_at_division, 2u);

        // Iterate over the cells, checking their target areas
        for (VertexBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double target_area_at_division = cell_iter->GetCellData()->GetItem("target area");
            TS_ASSERT_DELTA(target_area_at_division, 2.3, 1e-9);
        }
    }

    void TestConstantTargetAreaModifierArchiving()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ConstantTargetAreaModifier.arch";

        // Separate scope to write the archive
        {
            // Initialise a growth modifier and set a non-standard mature target area
            AbstractCellBasedSimulationModifier<2,2>* const p_modifier = new ConstantTargetAreaModifier<2>();
            (static_cast<ConstantTargetAreaModifier<2>*>(p_modifier))->SetReferenceTargetArea(14.3);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            output_arch << p_modifier;
            delete p_modifier;
        }

        // Separate scope to read the archive
        {
            AbstractCellBasedSimulationModifier<2,2>* p_modifier2;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_modifier2;

            // See whether we read out the correct member variables
            double reference_target_area = (static_cast<ConstantTargetAreaModifier<2>*>(p_modifier2))->GetReferenceTargetArea();
            TS_ASSERT_DELTA(reference_target_area, 14.3, 1e-9);

            delete p_modifier2;
        }
    }

    void TestConstantTargetAreaModifierOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestConstantTargetAreaModifierOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        MAKE_PTR(ConstantTargetAreaModifier<2>, p_modifier);
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "ConstantTargetAreaModifier-2");

        p_modifier->SetReferenceTargetArea(6.2);

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("ConstantTargetAreaModifier.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code
            FileFinder generated = output_file_handler.FindFile("ConstantTargetAreaModifier.parameters");
            FileFinder reference("cell_based/test/data/TestSimulationModifierOutputParameters/ConstantTargetAreaModifier.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTCONSTANTTARGETAREAMODIFIER_HPP_*/
