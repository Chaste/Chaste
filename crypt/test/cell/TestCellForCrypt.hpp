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

#ifndef TESTCELL_HPP_
#define TESTCELL_HPP_

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <iostream>

#include "Cell.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestCell: public AbstractCellBasedTestSuite
{
public:

    /*
     * ReadyToDivide() now calls UpdateCellProliferativeType() where appropriate.
     * (at the moment in Wnt-dependent cells).
     */
    void TestUpdateCellProliferativeTypes()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200, 20);
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();
        p_stem_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_stem_cell->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);

        p_stem_cell->SetCellProliferativeType(p_transit_type);

        p_stem_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_stem_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        // Test a Wnt dependent cell
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(0.0);

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);

        CellPtr p_wnt_cell(new Cell(p_healthy_state, p_cell_cycle_model1));
        p_wnt_cell->SetCellProliferativeType(p_transit_type);

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        p_wnt_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);

        p_wnt_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);

        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(1.0);

        // Go forward through time
        for (unsigned i=0; i<20; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }

        p_wnt_cell->ReadyToDivide();

        TS_ASSERT_EQUALS(p_wnt_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /*
     * We are checking that the CellPtrs work with the Wnt cell-cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithWntCellCycleModel()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        CellPtr p_wnt_cell(new Cell(p_healthy_state, p_cell_cycle_model1));
        p_wnt_cell->SetCellProliferativeType(p_transit_type);
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
    void TestWithStochasticWntCellCycleModel()
    {
        // If random number generation changes, then print these three lines to get the numbers to go in below
//        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(4, 0.9) << "\n";
//        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(4, 0.9) << "\n";
//        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(4, 0.9) << "\n";

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 101;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(50.0, num_steps+1);

        double wnt_stimulus = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_stimulus);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        StochasticWntCellCycleModel* p_cell_model = new StochasticWntCellCycleModel();
        p_cell_model->SetDimension(2);
        CellPtr p_wnt_cell(new Cell(p_healthy_state, p_cell_model));
        p_wnt_cell->SetCellProliferativeType(p_transit_type);
        p_wnt_cell->InitialiseCellCycleModel();

        // These are the first three normal random with mean of usual G2 Duration (4hrs), s.d. 0.9 and this seed (0)
        // Insert new random numbers in the middle of here!
        double SG2MDuration1 = p_cell_model->GetSDuration() + 3.17399 + p_cell_model->GetMDuration();
        double SG2MDuration2 = p_cell_model->GetSDuration() + 5.09655 + p_cell_model->GetMDuration();
        double SG2MDuration3 = p_cell_model->GetSDuration() + 5.55918 + p_cell_model->GetMDuration();
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
    void TestWntMutantVariantsAndLabelling()
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
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        WntCellCycleModel* p_cell_cycle_model1 = new WntCellCycleModel();
        p_cell_cycle_model1->SetDimension(2);
        CellPtr p_wnt_cell(new Cell(p_apc_one_hit_state, p_cell_cycle_model1));
        p_wnt_cell->SetCellProliferativeType(p_transit_type);
        p_wnt_cell->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model2 = new WntCellCycleModel();
        p_cell_cycle_model2->SetDimension(2);
        CellPtr p_wnt_cell2(new Cell(p_bcat_one_hit_state, p_cell_cycle_model2));
        p_wnt_cell2->SetCellProliferativeType(p_transit_type);
        p_wnt_cell2->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model3 = new WntCellCycleModel();
        p_cell_cycle_model3->SetDimension(2);
        CellPtr p_wnt_cell3(new Cell(p_apc_two_hit_state, p_cell_cycle_model3));
        p_wnt_cell3->SetCellProliferativeType(p_transit_type);
        p_wnt_cell3->InitialiseCellCycleModel();

        WntCellCycleModel* p_cell_cycle_model4 = new WntCellCycleModel();
        p_cell_cycle_model4->SetDimension(2);
        CellPropertyCollection collection;
        collection.AddProperty(p_label);
        CellPtr p_wnt_cell4(new Cell(p_wt_state, p_cell_cycle_model4, NULL, false, collection));
        p_wnt_cell4->SetCellProliferativeType(p_transit_type);
        p_wnt_cell4->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_wnt_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell3->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_wnt_cell4->ReadyToDivide(), false);

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTCELL_HPP_*/
