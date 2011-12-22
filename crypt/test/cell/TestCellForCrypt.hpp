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

#ifndef TESTCELL_HPP_
#define TESTCELL_HPP_

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <iostream>

#include "Cell.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "StochasticWntCellCycleModel.hpp"


class TestCell: public AbstractCellBasedTestSuite
{
public:

    /*
     * ReadyToDivide() now calls UpdateCellProliferativeType() where appropriate.
     * (at the moment in Wnt-dependent cells).
     */
    void TestUpdateCellProliferativeTypes() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200, 20);
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
        p_stem_cell->InitialiseCellCycleModel();
        p_stem_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->GetCellProliferativeType(),STEM);

        p_stem_cell->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);

        p_stem_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);

        // Test a Wnt dependent cell
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(0.0);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        CellPtr p_wnt_cell(new Cell(p_healthy_state, p_cell_cycle_model1));

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);

        p_wnt_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),DIFFERENTIATED);

        p_wnt_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),DIFFERENTIATED);

        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(1.0);

        // Go forward through time
        for (unsigned i=0; i<20; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }

        p_wnt_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the CellPtrs work with the Wnt cell-cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModel() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        CellPtr p_wnt_cell(new Cell(p_healthy_state, p_cell_cycle_model1));
        p_wnt_cell->InitialiseCellCycleModel();

        double SG2MDuration = p_cell_cycle_model1->GetSG2MDuration();

#ifdef CHASTE_CVODE
        const double expected_g1_duration = 5.96441;
#else
        const double expected_g1_duration = 5.971;
#endif //CHASTE_CVODE

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= expected_g1_duration+SG2MDuration)
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

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

            if (time >= expected_g1_duration + SG2MDuration + time_of_birth)
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
     * We are checking that the CellPtrs work with the StochasticWnt cell-cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithStochasticWntCellCycleModel() throw(Exception)
    {
//        double first_random_num = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
//        double second_random_num = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
//        double third_random_num = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
//        RandomNumberGenerator::Instance()->Reseed(0);


        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 101;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        StochasticWntCellCycleModel* p_cell_model = new StochasticWntCellCycleModel();
        p_cell_model->SetDimension(2);
        p_cell_model->SetCellProliferativeType(TRANSIT);
        CellPtr p_wnt_cell(new Cell(p_healthy_state, p_cell_model));
        p_wnt_cell->InitialiseCellCycleModel();

        // These are the first three normal random with mean of usual G2 Duration (4hrs), s.d. 0.9 and this seed (0)
        double SG2MDuration1 = p_cell_model->GetSDuration() + 5.00104 + p_cell_model->GetMDuration();
        double SG2MDuration2 = p_cell_model->GetSDuration() + 3.68921 + p_cell_model->GetMDuration();
        double SG2MDuration3 = p_cell_model->GetSDuration() + 4.54725 + p_cell_model->GetMDuration();
        double g1_duration = 5.971;

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            if (time >= g1_duration+SG2MDuration1)
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

            bool parent_ready = p_wnt_cell->ReadyToDivide();
            bool daughter_ready = p_wnt_cell2->ReadyToDivide();

            if (time >= g1_duration+SG2MDuration2+time_of_birth)
            {
                TS_ASSERT_EQUALS(parent_ready, true);
            }
            else
            {
                TS_ASSERT_EQUALS(parent_ready, false);
            }
            if (time >= g1_duration+SG2MDuration3+time_of_birth2)
            {
                TS_ASSERT_EQUALS(daughter_ready, true);
            }
            else
            {
                TS_ASSERT_EQUALS(daughter_ready, false);
            }
        }

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the CellPtrs work with the Wnt
     * cell-cycle models here. This just tests the set-up and checks that
     * the functions can all be called (not what they return).
     *
     * For more in depth tests see TestNightlyCellPtr.hpp
     * (these test that the cell cycle times are correct for the
     * various mutant cells)
     */
    void TestWntMutantVariantsAndLabelling() throw(Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 10;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        boost::shared_ptr<AbstractCellProperty> p_wt_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc_one_hit_state(CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc_two_hit_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat_one_hit_state(CellPropertyRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        p_cell_cycle_model1->SetCellProliferativeType(TRANSIT);
        CellPtr p_wnt_cell(new Cell(p_apc_one_hit_state, p_cell_cycle_model1));
        p_wnt_cell->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model2 = new WntCellCycleModel();
        p_cell_cycle_model2->SetDimension(2);
        p_cell_cycle_model2->SetCellProliferativeType(TRANSIT);
        CellPtr p_wnt_cell2(new Cell(p_bcat_one_hit_state, p_cell_cycle_model2));
        p_wnt_cell2->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model3 = new WntCellCycleModel();
        p_cell_cycle_model3->SetDimension(2);
        p_cell_cycle_model3->SetCellProliferativeType(TRANSIT);
        CellPtr p_wnt_cell3(new Cell(p_apc_two_hit_state, p_cell_cycle_model3));
        p_wnt_cell3->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model4 = new WntCellCycleModel();
        p_cell_cycle_model4->SetDimension(2);
        p_cell_cycle_model4->SetCellProliferativeType(TRANSIT);
        CellPropertyCollection collection;
        collection.AddProperty(p_label);
        CellPtr p_wnt_cell4(new Cell(p_wt_state, p_cell_cycle_model4, false, collection));
        p_wnt_cell4->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell3->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell4->ReadyToDivide(), false);

        //Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTCELL_HPP_*/
