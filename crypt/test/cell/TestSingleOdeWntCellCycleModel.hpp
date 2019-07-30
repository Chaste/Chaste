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

#ifndef TESTSINGLEODEWNTCELLCYCLEMODEL_HPP_
#define TESTSINGLEODEWNTCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "SingleOdeWntCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSingleOdeWntCellCycleModel : public AbstractCellBasedTestSuite
{
private:
    // Just update these if random number generation changes (error messages will show new values)
    static const double mFirstRandomNumber;
    static const double mSecondRandomNumber;
    static const double mThirdRandomNumber;
    static const double mFourthRandomNumber;
    static const double mFifthRandomNumber;

public:

    void xTestSetTheRandomNumbersForAllTheseTests()
    {
        // If random numbers change then copy the output of these to the static definitions at the bottom.
        // (Some beta-catenin concentrations might need changing manually too)

        double mean = 2.0;
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
    }

    void TestCorrectBehaviour()
    {
        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 1200.0;
        unsigned num_timesteps = 10*(unsigned)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Set up Wnt concentration
        double wnt_level = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell mutation state
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

        // Create cell-cycle model
        SingleOdeWntCellCycleModel* p_cycle_model = new SingleOdeWntCellCycleModel();
        p_cycle_model->SetDimension(2);

        // Construct a cell with this cell-cycle model and cell mutation state
        CellPtr p_cell(new Cell(p_state, p_cycle_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

        // Test the cell-cycle model is behaving correctly
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Stem cell should have been changed into a transit cell by wnt cell-cycle model
            TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

            // The number for the G1 duration is taken from the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, mFirstRandomNumber);
        }

        double steady_beta_cat_at_wnt_equals_1 = p_cycle_model->GetBetaCateninConcentration();

#ifdef CHASTE_CVODE
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, 143.4973, 2e-3);
        // (Different CVODE versions give slightly different answers here...)
#else
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, 143.5119, 1e-4);
#endif

        // Divide the cell
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
        CellPtr p_cell2 = p_cell->Divide();
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);
        p_cell->AddCellProperty(p_label);

        SingleOdeWntCellCycleModel* p_cycle_model2 = static_cast<SingleOdeWntCellCycleModel*> (p_cell2->GetCellCycleModel());

        // Test that both cells have inherited the same Beta-catenin concentration.
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, p_cycle_model->GetBetaCateninConcentration(), 1e-12);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_1, p_cycle_model2->GetBetaCateninConcentration(), 1e-12);

        // Now reduce the Wnt concentration
        wnt_level = 0.2;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // The numbers for the G1 durations are taken from
        // the first two random numbers generated
        double new_g1_duration = mSecondRandomNumber;
        double new_g1_duration2 = mThirdRandomNumber;
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }

        TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
        TS_ASSERT_EQUALS(p_cell2->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        // Test that both cells have inherited the same beta-catenin concentration
        double steady_beta_cat_at_wnt_equals_0_2 = 113.4683;
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);

        p_cycle_model->ResetForDivision();
        p_cycle_model2->ResetForDivision();

        // Test that both cells have still inherited the same ceta-catenin concentration
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
        TS_ASSERT_DELTA(steady_beta_cat_at_wnt_equals_0_2, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);

        // Now reduce the Wnt concentration so only mutant cells divide
        wnt_level = 0.1;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Introduce a mutation (no immediate effect)
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        p_cell2->SetMutationState(p_apc1);
        TS_ASSERT_EQUALS(p_cycle_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_cycle_model2->ReadyToDivide(), false);

        // Coverage
        p_cell2->SetMutationState(p_apc2);
        TS_ASSERT_EQUALS(p_cycle_model->ReadyToDivide(), false);

        p_cell2->SetMutationState(p_bcat1);
        TS_ASSERT_EQUALS(p_cycle_model->ReadyToDivide(), false);

        // The numbers for the G1 durations are taken from the next two random numbers generated
        new_g1_duration = mFourthRandomNumber;
        new_g1_duration2 = mFifthRandomNumber;

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }

        TS_ASSERT_DELTA(91.6693, p_cycle_model->GetBetaCateninConcentration(), 1e-3);
#ifdef CHASTE_CVODE
        TS_ASSERT_DELTA(361.5455, p_cycle_model2->GetBetaCateninConcentration(), 5e-1);
        // Different versions of CVODE giving different answers here.
#else
        TS_ASSERT_DELTA(361.4677, p_cycle_model2->GetBetaCateninConcentration(), 1e-3);
#endif

        TS_ASSERT_DELTA(p_cycle_model->GetBetaCateninDivisionThreshold(), 100, 1e-9);
        TS_ASSERT_DELTA(p_cycle_model2->GetBetaCateninDivisionThreshold(), 100, 1e-9);

        TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);
        TS_ASSERT_EQUALS(p_cell2->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        // Coverage of 1D

        // Set up SimulationTime
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(30.0, 2);

        // Instantiate 1D Wnt concentration
        wnt_level = 1.0;
        WntConcentration<1>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell-cycle model and cell and test behaviour
        SingleOdeWntCellCycleModel* p_cell_model_1d = new SingleOdeWntCellCycleModel;
        p_cell_model_1d->SetDimension(1);
        p_cell_model_1d->SetUseCellProliferativeTypeDependentG1Duration();
        TS_ASSERT_EQUALS(p_cell_model_1d->GetDimension(), 1u);

        CellPtr p_stem_cell_1d(new Cell(p_state, p_cell_model_1d));
        p_stem_cell_1d->SetCellProliferativeType(p_stem_type);
        p_stem_cell_1d->InitialiseCellCycleModel();

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell_1d->ReadyToDivide(), false);

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell_1d->ReadyToDivide(), true);

        p_stem_cell_1d->Divide();

        // Coverage of 3D

        // Set up SimulationTime
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(20.0, 2);

        // Instantiate 3D Wnt concentration
        WntConcentration<3>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Create cell-cycle model and cell and test behaviour
        SingleOdeWntCellCycleModel* p_cell_model_3d = new SingleOdeWntCellCycleModel;
        p_cell_model_3d->SetDimension(3);

        TS_ASSERT_EQUALS(p_cell_model_3d->GetDimension(), 3u);

        CellPtr p_stem_cell_3d(new Cell(p_state, p_cell_model_3d));
        p_stem_cell_3d->SetCellProliferativeType(p_stem_type);
        p_stem_cell_3d->InitialiseCellCycleModel();

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell_3d->ReadyToDivide(), false);

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell_3d->ReadyToDivide(), true);

        p_stem_cell_3d->Divide();

        // Tidy up
        WntConcentration<1>::Destroy();
        WntConcentration<2>::Destroy();
        WntConcentration<3>::Destroy();
    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "single_ode_wnt.arch";

        double random_number_test = 0;

        double beta_catenin_threshold = 0.76;

        // Create an output archive
        {
            // The number for the G1 duration is taken from the first random number generated
            double g1_duration = mFirstRandomNumber;

            // Set up the Wnt concentration for testing
            WntConcentration<2>::Instance()->SetConstantWntValueForTesting(0.7);

            // Create cell-cycle model and associated cell
            SingleOdeWntCellCycleModel* p_cell_model = new SingleOdeWntCellCycleModel;
            p_cell_model->SetDimension(2);
            p_cell_model->SetBirthTime(-1.0);

            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

            CellPtr p_stem_cell(new Cell(p_healthy_state, p_cell_model));
            p_stem_cell->SetCellProliferativeType(p_stem_type);
            p_stem_cell->InitialiseCellCycleModel();

            // Set up the simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            double end_time = g1_duration + p_cell_model->GetSG2MDuration() + 5.0;
            unsigned num_timesteps = 50;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

            while (p_cell_model->GetAge() < g1_duration + p_cell_model->GetSG2MDuration()
                    - p_simulation_time->GetTimeStep()) // minus one to match birth time.
            {
                p_simulation_time->IncrementTimeOneStep();
                CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, g1_duration);
            }

            // Wnt should change this to a transit cell
            TS_ASSERT_EQUALS(p_stem_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
            TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(static_cast<SingleOdeWntCellCycleModel*>(p_stem_cell->GetCellCycleModel())->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            p_cell_model->SetBetaCateninDivisionThreshold(beta_catenin_threshold);
            TS_ASSERT_DELTA(p_cell_model->GetBetaCateninDivisionThreshold(), beta_catenin_threshold, 1e-12);

            // Archive the cell
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            CellPtr const p_const_cell = p_stem_cell;
            output_arch << p_const_cell;

            TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->ReadyToDivide(), true);

            // Tidy up
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }

        {
            // Set up SimulationTime
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            CellPtr p_cell;

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(36);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check
            SingleOdeWntCellCycleModel* p_cell_model = static_cast<SingleOdeWntCellCycleModel*>(p_cell->GetCellCycleModel());
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_DELTA(p_cell_model->GetBetaCateninDivisionThreshold(), beta_catenin_threshold, 1e-12);

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);

            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetSG2MDuration(), 10.0, 1e-12);

            TS_ASSERT_DELTA(p_gen->ranf(), random_number_test, 1e-7);
            TS_ASSERT_EQUALS(p_cell_model->GetDimension(), 2u);

            // Tidy up
            SimulationTime::Destroy();
        }
    }

    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with SingleOdeWntCellCycleModel
        SingleOdeWntCellCycleModel single_ode_wnt_cell_cycle_model;
        TS_ASSERT_EQUALS(single_ode_wnt_cell_cycle_model.GetIdentifier(), "SingleOdeWntCellCycleModel");

        out_stream single_ode_wnt_parameter_file = output_file_handler.OpenOutputFile("single_ode_wnt_results.parameters");
        single_ode_wnt_cell_cycle_model.OutputCellCycleModelParameters(single_ode_wnt_parameter_file);
        single_ode_wnt_parameter_file->close();

        std::string single_ode_wnt_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison( single_ode_wnt_results_dir + "single_ode_wnt_results.parameters", "crypt/test/data/TestCellCycleModels/single_ode_wnt_results.parameters").CompareFiles();
    }
};

// Member initialisation
const double TestSingleOdeWntCellCycleModel::mFirstRandomNumber = 1.08221;
const double TestSingleOdeWntCellCycleModel::mSecondRandomNumber = 3.21839;
const double TestSingleOdeWntCellCycleModel::mThirdRandomNumber = 3.73243;
const double TestSingleOdeWntCellCycleModel::mFourthRandomNumber = 2.83804;
const double TestSingleOdeWntCellCycleModel::mFifthRandomNumber = 1.7031;

#endif /* TESTSINGLEODEWNTCELLCYCLEMODEL_HPP_ */
