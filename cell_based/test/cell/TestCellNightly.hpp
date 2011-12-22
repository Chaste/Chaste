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

#ifndef TESTCELLNIGHTLY_HPP_
#define TESTCELLNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <iostream>

#include "Cell.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestCellNightly: public AbstractCellBasedTestSuite
{
public:

    void TestTysonNovakImmortalStemCell()
    {
        double end_time = 100.0; // A good load of divisions to make sure nothing mucks up..
        // one division = 1.26 hours.
        int time_steps = 1000;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);

        unsigned divisions = 0;
        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();

            if (p_stem_cell->ReadyToDivide())
            {
                p_stem_cell->Divide();
                TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
                divisions++;
            }
        }
        TS_ASSERT_DELTA(divisions, (unsigned)(end_time/1.26), 1);
    }

    void Test0DBucketWithTysonNovak()
    {
        double end_time = 7.0; // not very long because cell cycle time is only 1.2
        // (75 mins) because Tyson Novaks is for yeast
        int time_steps = 100;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();

        std::vector<CellPtr> cells;
        std::vector<CellPtr> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<double> times(time_steps);

        cells.push_back(p_stem_cell);
        std::vector<CellPtr>::iterator cell_iterator;

        unsigned i = 0;
        while (p_simulation_time->GetTime()< end_time)
        {
            // Produce the offspring of the cells
            p_simulation_time->IncrementTimeOneStep();
            times[i] = p_simulation_time->GetTime();
            cell_iterator = cells.begin();
            unsigned j = 0;
            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->ReadyToDivide())
                {
                    newly_born.push_back((*cell_iterator)->Divide());
                }
                ++cell_iterator;
                j++;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Update cell counts
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                switch ((*cell_iterator)->GetCellCycleModel()->GetCellProliferativeType())
                {
                    case STEM:
                        stem_cells[i]++;
                        break;
                    case TRANSIT:
                        transit_cells[i]++;
                        break;
                    default:
                        differentiated_cells[i]++;
                        break;
                }

                ++cell_iterator;
            }

            i++;
        }

        TS_ASSERT_EQUALS(stem_cells[time_steps-1], 1u);
        TS_ASSERT_EQUALS(transit_cells[time_steps-1], 31u);
        TS_ASSERT_EQUALS(differentiated_cells[time_steps-1], 0u);
    }
};

#endif /*TESTCELLNIGHTLY_HPP_*/
