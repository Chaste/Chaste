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

#ifndef TESTCELLNIGHTLYFORCRYPT_HPP_
#define TESTCELLNIGHTLYFORCRYPT_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <iostream>

#include "Cell.hpp"
#include "WntCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestCellNightlyForCrypt: public AbstractCellBasedTestSuite
{
public:
    /*
     * We are checking that the CellPtrs work with the Wnt cell-cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationAPCONEHIT()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_transit_type(new TransitCellProliferativeType);
        CellPtr p_wnt_cell(new Cell(p_apc1, p_cell_cycle_model1));
        p_wnt_cell->SetCellProliferativeType(p_transit_type);
        p_wnt_cell->InitialiseCellCycleModel();

        double s_g2_duration = p_cell_cycle_model1->GetSG2MDuration();

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= 4.804 + s_g2_duration)
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);

        CellPtr p_wnt_cell2 = p_wnt_cell->Divide();

        double time_of_birth = p_wnt_cell->GetBirthTime();
        double time_of_birth2 = p_wnt_cell2->GetBirthTime();

        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result1 = p_wnt_cell->ReadyToDivide();
            bool result2 = p_wnt_cell2->ReadyToDivide();

            if (time >= 4.804 + s_g2_duration + time_of_birth)
            {
                TS_ASSERT_EQUALS(result1, true);
                TS_ASSERT_EQUALS(result2, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result1, false);
                TS_ASSERT_EQUALS(result2, false);
            }
        }

        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the CellPtrs work with the Wnt cell-cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationBetaCat()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 0.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_transit_type(new TransitCellProliferativeType);
        CellPtr p_wnt_cell(new Cell(p_bcat1, p_cell_cycle_model1));
        p_wnt_cell->SetCellProliferativeType(p_transit_type);
        p_wnt_cell->InitialiseCellCycleModel();

        double s_g2_duration = p_cell_cycle_model1->GetSG2MDuration();

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= 7.82 + s_g2_duration)
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);

        CellPtr p_wnt_cell2 = p_wnt_cell->Divide();

        double time_of_birth = p_wnt_cell->GetBirthTime();
        double time_of_birth2 = p_wnt_cell2->GetBirthTime();

        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result1 = p_wnt_cell->ReadyToDivide();
            bool result2 = p_wnt_cell2->ReadyToDivide();

            if (time >= 7.82 + s_g2_duration + time_of_birth)
            {
                TS_ASSERT_EQUALS(result1, true);
                TS_ASSERT_EQUALS(result2, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result1, false);
                TS_ASSERT_EQUALS(result2, false);
            }
        }

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the CellPtrs work with the Wnt cell-cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModelAndMutationAPC2()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 200;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 0.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_transit_type(new TransitCellProliferativeType);
        CellPtr p_wnt_cell(new Cell(p_apc2, p_cell_cycle_model1));
        p_wnt_cell->SetCellProliferativeType(p_transit_type);
        p_wnt_cell->InitialiseCellCycleModel();

        boost::shared_ptr<AbstractCellMutationState> p_this_state = p_wnt_cell->GetMutationState();

        TS_ASSERT_EQUALS(p_this_state->IsType<ApcTwoHitCellMutationState>(), true);

        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);

        p_wnt_cell->SetMutationState(p_apc1);

        p_this_state = p_wnt_cell->GetMutationState();

        TS_ASSERT_EQUALS(p_this_state->IsType<ApcTwoHitCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_this_state->IsType<ApcOneHitCellMutationState>(), true);

        p_wnt_cell->SetMutationState(p_apc2);

        double s_g2_duration = p_cell_cycle_model1->GetSG2MDuration();

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= 3.9435 + s_g2_duration)
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);

        CellPtr p_wnt_cell2 = p_wnt_cell->Divide();

        double time_of_birth = p_wnt_cell->GetBirthTime();
        double time_of_birth2 = p_wnt_cell2->GetBirthTime();

        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result1 = p_wnt_cell->ReadyToDivide();
            bool result2 = p_wnt_cell2->ReadyToDivide();

            if (time >= 3.9435 + s_g2_duration + time_of_birth)
            {
                TS_ASSERT_EQUALS(result1, true);
                TS_ASSERT_EQUALS(result2, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result1, false);
                TS_ASSERT_EQUALS(result2, false);
            }
        }

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTCELLNIGHTLYFORCRYPT_HPP_*/
