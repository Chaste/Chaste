/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTCABASEDDIVISIONRULE_HPP_
#define TESTCABASEDDIVISIONRULE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "CaBasedCellPopulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCaBasedDivisionRule.hpp"
#include "ExclusionCaBasedDivisionRule.hpp"
#include "PottsMeshGenerator.hpp"
#include "SmartPointers.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCaBasedDivisionRules : public AbstractCellBasedTestSuite
{
public:

    void TestAddCellwithExclusionBasedDivisionRule()
    {
        /**
         * In this test we basically test that the AbstractCaBasedDivisionRule is implemented and joined with the population
         * correctly. We make a new ExclusionCaBasedDivisionRule, divide a cell with it and check that the new cells
         * are in the correct locations.
         */

        // Make a simple Potts mesh
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create 6 cells in the bottom 2 rows
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<6; index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_temp_cell(new Cell(p_state, p_model));
        p_temp_cell->SetCellProliferativeType(p_stem_type);
        p_temp_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the exclusion division rule (note that this is the default

        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new ExclusionCaBasedDivisionRule<2>());
        cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and try to add new cell by dividing cell at site 0;
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();


        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT(!(p_division_rule->IsRoomToDivide(p_parent_cell,cell_population)));

        // Test adding the new cell in the population (note this calls CalculateDaughterNodeIndex)
        TS_ASSERT_THROWS_THIS(cell_population.AddCell(p_cell_0, p_parent_cell),
                              "Trying to divide when there is no room to divide, check your division rule");

        // Test adding it in a free space
        p_parent_cell = cell_population.GetCellUsingLocationIndex(4);
        TS_ASSERT(p_division_rule->IsRoomToDivide(p_parent_cell,cell_population));

        TS_ASSERT_EQUALS(p_division_rule->CalculateDaughterNodeIndex(p_cell_0,p_parent_cell,cell_population), 7u);

        // Test adding the new cell in the population (note this calls CalculateDaughterNodeIndex)
        cell_population.AddCell(p_cell_0, p_parent_cell);

        // Now check the cells are in the correct place
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 7u);
    }
};

#endif /*TESTCABASEDDIVISIONRULE_HPP_*/
