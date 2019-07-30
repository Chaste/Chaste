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

#ifndef TESTCELLNIGHTLY_HPP_
#define TESTCELLNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <fstream>
#include <iostream>

#include "Cell.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

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
        MAKE_PTR(StemCellProliferativeType, p_type);
        TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
        p_stem_cell->SetCellProliferativeType(p_type);
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
        MAKE_PTR(StemCellProliferativeType, p_type);

        TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
        p_stem_cell->SetCellProliferativeType(p_type);
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
                if ((*cell_iterator)->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                {
                    stem_cells[i]++;
                }
                else if ((*cell_iterator)->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                {
                    transit_cells[i]++;
                }
                else
                {
                    differentiated_cells[i]++;
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
