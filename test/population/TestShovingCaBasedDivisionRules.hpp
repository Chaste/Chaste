/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTSHOVINGCABASEDDIVISIONRULE_HPP_
#define TESTSHOVINGCABASEDDIVISIONRULE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "CaBasedCellPopulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCaBasedDivisionRule.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "PottsMeshGenerator.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * This test is to be appended to TestCaBasedDivisionRules when moving to the cell_based and crypt folders.
 */

class TestShovingCaBasedDivisionRules : public AbstractCellBasedTestSuite
{
public:

	void TestTest()
	{
		TS_ASSERT(true);
	}

    void TestAddCellwithShovingBasedDivisionRule()
    {
        /**
         * In this test we basically test that the AbstractCaBasedDivisionRule is implemented and joined with the population
         * correctly. We make a new ShovingCaBasedDivisionRule, divide a cell with it and check that the new cells
         * are in the correct locations.
         */

        /*
         * First we test where this is space around the cells.
         * This is the default setup.
         */
        {
            // Make a simple potts mesh
            PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
            PottsMesh<2>* p_mesh = generator.GetMesh();

            // Create 9 cells in the central nodes
            std::vector<unsigned> location_indices;
            for (unsigned row=1; row<4; row++)
            {
                location_indices.push_back(1+row*5);
                location_indices.push_back(2+row*5);
                location_indices.push_back(3+row*5);
            }

            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator;
            cells_generator.GenerateBasic(cells, location_indices.size());

            // Create cell population
            CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

            unsigned cell_locations[9] = {6, 7, 8, 11, 12, 13, 16, 17, 18};
            unsigned index = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),cell_locations[index])
                ++index;
            }

            // Make a new cell to add
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);

            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_new_cell(new Cell(p_state, p_model));
            p_new_cell->SetCellProliferativeType(p_stem_type);
            p_new_cell->SetBirthTime(-1);

            // Set the division rule for our population to be the shoving division rule
            boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new ShovingCaBasedDivisionRule<2>());
            cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

            // Get the division rule back from the population and try to add new cell by dividing cell at site 0;
            boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

            // Select central cell
            CellPtr p_cell_12 = cell_population.GetCellUsingLocationIndex(12);
            // This always returns true for this division rule
            TS_ASSERT((p_division_rule->IsRoomToDivide(p_cell_12,cell_population)));


            // Moves into node 13
            TS_ASSERT_EQUALS(p_division_rule->CalculateDaughterNodeIndex(p_new_cell,p_cell_12,cell_population), 13u);

            // Test adding the new cell in the population (note this calls CalculateDaughterNodeIndex)
            cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_cell_12);

            // Now check the cells are in the correct place
            TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 10u);

            // Note the cell on node 13 has been shoved to node 14 and the new cell is on node 13
            unsigned new_cell_locations[10] = {6, 7, 8, 11, 12, 14, 16, 17, 18, 13};
            index = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),new_cell_locations[index])
                ++index;
            }
        }
        /*
         * We now test where there is no room to divide without the cells being shoved to the edge of the mesh.
         */
        {
            // Make a simple potts mesh
            PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
            PottsMesh<2>* p_mesh = generator.GetMesh();

            // Create 25 cells, one for each node
            std::vector<unsigned> location_indices;
            for (unsigned index=0; index<25; index++)
            {
                location_indices.push_back(index);
            }

            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator;
            cells_generator.GenerateBasic(cells, location_indices.size());

            // Create cell population
            CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

            // Make a new cell to add
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);

            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_new_cell(new Cell(p_state, p_model));
            p_new_cell->SetCellProliferativeType(p_stem_type);
            p_new_cell->SetBirthTime(-1);

            // Set the division rule for our population to be the shoving division rule
            boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new ShovingCaBasedDivisionRule<2>());
            cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

            // Get the division rule back from the population and try to add new cell by dividing cell at site 0;
            boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

            // Select central cell
            CellPtr p_cell_12 = cell_population.GetCellUsingLocationIndex(12);
            // This always returns true for this division rule
            TS_ASSERT((p_division_rule->IsRoomToDivide(p_cell_12,cell_population)));


            // Moves into node 13
            TS_ASSERT_THROWS_THIS(p_division_rule->CalculateDaughterNodeIndex(p_new_cell,p_cell_12,cell_population),
                    "Cells reaching the boundary of the domain. Make the potts mesh larger.");
        }
    }
};

#endif /*TESTSHOVINGCABASEDDIVISIONRULE_HPP_*/
