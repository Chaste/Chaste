/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTCELLCYCLEMODELSNOTFORRELEASE_HPP_
#define TESTCELLCYCLEMODELSNOTFORRELEASE_HPP_


#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "StochasticDivisionRuleCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"

/**
 * This class contains tests for methods on cell
 * cycle models that are not yet ready for release.
 */
class TestCellCycleModelsNotForRelease : public AbstractCellBasedTestSuite
{
public:

    void TestStochasticDivisionRuleCellCycleModel() throw(Exception)
    {
        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 72.0;
        unsigned num_timesteps = 1000*(unsigned)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Create a cell
        StochasticDivisionRuleCellCycleModel* p_cycle_model1 = new StochasticDivisionRuleCellCycleModel;
        p_cycle_model1->SetCellProliferativeType(STEM);
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TissueCellPtr p_cell1(new TissueCell(p_state, p_cycle_model1));
        p_cell1->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cycle_model1->GetGeneration(), 0u);

        /**
         * Test with asymmetric division
         */

        p_cycle_model1->SetSymmetricDivisionProbability(0.0);

        // Increment time
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model1, 13.0676);
        }

        // This cell must have divided asymmetrically
        TS_ASSERT_EQUALS(p_cycle_model1->DividedSymmetrically(), false);
        TS_ASSERT_EQUALS(p_cell1->GetCellCycleModel()->GetCellProliferativeType(), STEM);
        TS_ASSERT_EQUALS(p_cycle_model1->GetGeneration(), 0u);

        TS_ASSERT_EQUALS(p_cell1->ReadyToDivide(), true);
        TissueCellPtr p_cell2 = p_cell1->Divide();

        TS_ASSERT_EQUALS(p_cell2->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);

        StochasticDivisionRuleCellCycleModel* p_cycle_model2 = static_cast<StochasticDivisionRuleCellCycleModel*> (p_cell2->GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model2->GetGeneration(), 1u);

        /**
         * Test with symmetric division
         */

        p_cycle_model1->SetSymmetricDivisionProbability(1.0);
        TissueConfig::Instance()->SetMaxTransitGenerations(1);

        StochasticDivisionRuleCellCycleModel* p_cycle_model3 = new StochasticDivisionRuleCellCycleModel;
        p_cycle_model3->SetCellProliferativeType(STEM);        
        p_cycle_model3->SetSymmetricDivisionProbability(1.0);

        TissueCellPtr p_cell3(new TissueCell(p_state, p_cycle_model3));
        p_cell3->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cycle_model3->GetGeneration(), 0u);
        TS_ASSERT_DELTA(p_cycle_model3->GetSymmetricDivisionProbability(), 1.0, 1e-6);

        // Increment time
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model1, 13.2712);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, 1.22037);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model3, 12.747);
        }

        // The stem cell cell1 must have divided symmetrically, and it
        // happens to have divided into two stem cells
        TS_ASSERT_EQUALS(p_cell1->ReadyToDivide(), true);
        TissueCellPtr p_cell4 = p_cell1->Divide();

        TS_ASSERT_EQUALS(p_cycle_model1->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(p_cell1->GetCellCycleModel()->GetCellProliferativeType(), STEM);
        TS_ASSERT_EQUALS(p_cycle_model1->GetGeneration(), 0u);

        StochasticDivisionRuleCellCycleModel* p_cycle_model4 = static_cast<StochasticDivisionRuleCellCycleModel*> (p_cell4->GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model4->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(p_cell4->GetCellCycleModel()->GetCellProliferativeType(), STEM);
        TS_ASSERT_EQUALS(p_cycle_model4->GetGeneration(), 0u);
        TS_ASSERT_DELTA(p_cycle_model4->GetSymmetricDivisionProbability(), 1.0, 1e-6);

        // The stem cell cell3 must have divided symmetrically. For coverage,
        // we iterate the random number generator so that it divides into two
        // transit cells
        RandomNumberGenerator::Instance()->ranf();

        TS_ASSERT_EQUALS(p_cell3->ReadyToDivide(), true);
        TissueCellPtr p_cell5 = p_cell3->Divide();

        TS_ASSERT_EQUALS(p_cycle_model3->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(p_cell3->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_EQUALS(p_cycle_model3->GetGeneration(), 1u);

        StochasticDivisionRuleCellCycleModel* p_cycle_model5 = static_cast<StochasticDivisionRuleCellCycleModel*> (p_cell5->GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model5->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(p_cell5->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_EQUALS(p_cycle_model5->GetGeneration(), 1u);

        // The transit cell cell2 divides into two differentiated cells
        TS_ASSERT_EQUALS(p_cell2->ReadyToDivide(), true);
        TissueCellPtr p_cell6 = p_cell2->Divide();

        TS_ASSERT_EQUALS(p_cycle_model2->DividedSymmetrically(), false);
        TS_ASSERT_EQUALS(p_cell2->GetCellCycleModel()->GetCellProliferativeType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(p_cycle_model2->GetGeneration(), 2u);

        StochasticDivisionRuleCellCycleModel* p_cycle_model6 = static_cast<StochasticDivisionRuleCellCycleModel*> (p_cell6->GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model6->DividedSymmetrically(), false);
        TS_ASSERT_EQUALS(p_cell6->GetCellCycleModel()->GetCellProliferativeType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(p_cycle_model6->GetGeneration(), 2u);

        // For coverage
        StochasticDivisionRuleCellCycleModel* p_cycle_model7 = new StochasticDivisionRuleCellCycleModel;
        p_cycle_model7->SetCellProliferativeType(DIFFERENTIATED);
        TissueCellPtr p_cell7(new TissueCell(p_state, p_cycle_model7));
        p_cell7->InitialiseCellCycleModel();
        TS_ASSERT_EQUALS(p_cycle_model7->GetGeneration(), 0u);
    }


    void TestArchiveStochasticDivisionRuleCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "stoch_div_rule_cell_cycle.arch";

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            StochasticDivisionRuleCellCycleModel model(true);
            model.SetCellProliferativeType(STEM);
            model.SetSymmetricDivisionProbability(0.5);

            p_simulation_time->IncrementTimeOneStep();

            model.SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << static_cast<const StochasticDivisionRuleCellCycleModel&>(model);

            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            StochasticDivisionRuleCellCycleModel model;
            model.SetBirthTime(-2.0);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> model;

            // Check that archiving worked correctly
            TS_ASSERT_EQUALS(model.GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(model.GetCellProliferativeType(), STEM);
            TS_ASSERT_EQUALS(model.DividedSymmetrically(), true);

            TS_ASSERT_DELTA(model.GetSymmetricDivisionProbability(), 0.5, 1e-6);
            TS_ASSERT_DELTA(model.GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(model.GetAge(), 1.5, 1e-12);
        }
    }

};

#endif /*TESTCELLCYCLEMODELSNOTFORRELEASE_HPP_*/
