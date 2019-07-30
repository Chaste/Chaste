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

#ifndef TESTSIMPLECELLCYCLEMODELSFORCRYPT_HPP_
#define TESTSIMPLECELLCYCLEMODELSFORCRYPT_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "Cell.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSimpleCellCycleModelsForCrypt : public AbstractCellBasedTestSuite
{
private:
    static const double mFirstRandomNumber;
    static const double mSecondRandomNumber;
    static const double mThirdRandomNumber;
    static const double mFourthRandomNumber;
    static const double mFifthRandomNumber;

public:

    void xTestSetTheRandomNumbersForAllTheseTests()
    {
        // If random numbers change then copy the output of these to the static definitions at the bottom.

        double mean = 2.0;
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
        std::cout << RandomNumberGenerator::Instance()->NormalRandomDeviate(mean, 1.0) << "\n";
    }

    void TestSimpleWntCellCycleModel()
    {
        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        double end_time = 60.0;
        unsigned num_timesteps = 1000*(unsigned)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2*end_time, 2*num_timesteps);

        // Set up the Wnt concentration
        double wnt_level = 1.0;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        SimpleWntCellCycleModel* p_cycle_model = new SimpleWntCellCycleModel;
        TS_ASSERT_EQUALS(p_cycle_model->GetDimension(), UNSIGNED_UNSET);
        TS_ASSERT_EQUALS(p_cycle_model->CanCellTerminallyDifferentiate(), false);

        // Test the dimension must be 1, 2, 3 or UNSIGNED_UNSET
        TS_ASSERT_THROWS_THIS(p_cycle_model->SetDimension(4), "Dimension must be 1, 2, 3 or UNSIGNED_UNSET");

        // Test the set/get dimension methods
        p_cycle_model->SetDimension(2);
        TS_ASSERT_EQUALS(p_cycle_model->GetDimension(), 2u);

        TS_ASSERT_DELTA(p_cycle_model->GetWntStemThreshold(), 0.8, 1e-6);
        TS_ASSERT_DELTA(p_cycle_model->GetWntTransitThreshold(), 0.65, 1e-6);
        TS_ASSERT_DELTA(p_cycle_model->GetWntLabelledThreshold(), 0.65, 1e-6);

        p_cycle_model->SetWntStemThreshold(0.4);
        p_cycle_model->SetWntTransitThreshold(0.5);
        p_cycle_model->SetWntLabelledThreshold(0.3);

        TS_ASSERT_DELTA(p_cycle_model->GetWntStemThreshold(), 0.4, 1e-6);
        TS_ASSERT_DELTA(p_cycle_model->GetWntTransitThreshold(), 0.5, 1e-6);
        TS_ASSERT_DELTA(p_cycle_model->GetWntLabelledThreshold(), 0.3, 1e-6);

        p_cycle_model->SetWntStemThreshold(0.8);
        p_cycle_model->SetWntTransitThreshold(0.65);
        p_cycle_model->SetWntLabelledThreshold(0.65);

        boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

        CellPtr p_cell(new Cell(p_healthy_state, p_cycle_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The number for the G1 duration is taken from
            // the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, mFirstRandomNumber);
        }

        // Stem cell should have been changed into a transit cell by wnt cell-cycle model
        TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        // Divide the cell
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
        CellPtr p_cell2 = p_cell->Divide();
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);
        p_cell->AddCellProperty(p_label);

        SimpleWntCellCycleModel* p_cycle_model2 = static_cast<SimpleWntCellCycleModel*> (p_cell2->GetCellCycleModel());

        // Now reduce the Wnt concentration
        wnt_level = 0.7;
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

        p_cycle_model->ResetForDivision();
        p_cycle_model2->ResetForDivision();

        // Now reduce the Wnt concentration so only beta-cat or APC2 hit cells divide.
        wnt_level = 0.15;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        boost::shared_ptr<AbstractCellMutationState> p_apc1_mutation(new ApcOneHitCellMutationState);
        p_cell->SetMutationState(p_apc1_mutation);
        boost::shared_ptr<AbstractCellMutationState> p_bcat_mutation(new BetaCateninOneHitCellMutationState);
        p_cell2->SetMutationState(p_bcat_mutation);

        TS_ASSERT_EQUALS(p_cycle_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_cycle_model2->ReadyToDivide(), false);

        // Coverage...
        boost::shared_ptr<AbstractCellMutationState> p_apc2_mutation(new ApcTwoHitCellMutationState);
        p_cell->SetMutationState(p_apc2_mutation);
        TS_ASSERT_EQUALS(p_cycle_model->ReadyToDivide(), false);
        p_cell->SetMutationState(p_apc1_mutation);
        TS_ASSERT_EQUALS(p_cycle_model->ReadyToDivide(), false);

        // The numbers for the G1 durations are taken from
        // the next two random numbers generated
        new_g1_duration = mFourthRandomNumber;
        new_g1_duration2 = mFifthRandomNumber;

        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model, new_g1_duration);
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model2, new_g1_duration2);
        }

        TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);
        TS_ASSERT_EQUALS(p_cell2->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        // For coverage...
        SimpleWntCellCycleModel* p_cycle_model1 = new SimpleWntCellCycleModel;
        p_cycle_model1->SetDimension(2);

        CellPtr p_cell1(new Cell(p_healthy_state, p_cycle_model1));
        p_cell1->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        p_cell1->InitialiseCellCycleModel();

        SimpleWntCellCycleModel* p_another_cycle_model = new SimpleWntCellCycleModel;
        p_another_cycle_model->SetDimension(2);

        CellPtr p_another_cell(new Cell(p_healthy_state, p_another_cycle_model));
        p_another_cell->SetCellProliferativeType(p_stem_type);
        p_another_cell->InitialiseCellCycleModel();
        // ...end of coverage

        // Test the case of a radial Wnt concentration

        RandomNumberGenerator::Instance()->Reseed(0);

        // Set up the Wnt concentration
        wnt_level = 0.81;
        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // Set up a cell-cycle model and cell
        SimpleWntCellCycleModel* p_cycle_model4 = new SimpleWntCellCycleModel;
        p_cycle_model4->SetDimension(2);

        CellPtr p_cell4(new Cell(p_healthy_state,  p_cycle_model4));
        p_cell4->SetCellProliferativeType(p_stem_type);
        p_cell4->InitialiseCellCycleModel();

        // Test the GetCurrentCellCyclePhase() and ReadyToDivide() methods
        double first_g1_duration = mFirstRandomNumber;
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The number for the G1 duration is taken from
            // the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model4, first_g1_duration);
        }

        // We should still have a stem cell since the WntConcentration exceeds mRadialWntThreshold
        TS_ASSERT_EQUALS(p_cell4->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);

        // Divide the cell
        TS_ASSERT_EQUALS(p_cell4->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_cell4->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);
        CellPtr p_cell5 = p_cell4->Divide();
        TS_ASSERT_EQUALS(p_cell4->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);
        TS_ASSERT_EQUALS(p_cell5->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        p_cell2->AddCellProperty(p_label);

        // Now reduce the Wnt concentration
        wnt_level = 0.79;
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        // The numbers for the G1 durations are taken from
        // the first two random numbers generated
        new_g1_duration = mSecondRandomNumber;
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model4, new_g1_duration);
        }

        TS_ASSERT_DELTA(WntConcentration<2>::Instance()->GetWntLevel(p_cell4), wnt_level, 1e-12);
        TS_ASSERT_EQUALS(p_cell4->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
        TS_ASSERT_EQUALS(p_cell5->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

        // Coverage of 1D

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(30.0, 2);

        WntConcentration<1>::Instance()->SetConstantWntValueForTesting(wnt_level);
        SimpleWntCellCycleModel* p_cell_model_1d = new SimpleWntCellCycleModel;
        p_cell_model_1d->SetDimension(1);
        p_cell_model_1d->SetUseCellProliferativeTypeDependentG1Duration();

        TS_ASSERT_EQUALS(p_cell_model_1d->GetDimension(), 1u);
        TS_ASSERT_EQUALS(p_cell_model_1d->GetUseCellProliferativeTypeDependentG1Duration(), true);

        CellPtr p_stem_cell_1d(new Cell(p_healthy_state, p_cell_model_1d));
        p_stem_cell_1d->SetCellProliferativeType(p_stem_type);
        p_stem_cell_1d->InitialiseCellCycleModel();

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell_1d->ReadyToDivide(), false);

        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell_1d->ReadyToDivide(), true);

        p_stem_cell_1d->Divide();

        // Coverage of 3D

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(20.0, 2);

        WntConcentration<3>::Instance()->SetConstantWntValueForTesting(wnt_level);
        SimpleWntCellCycleModel* p_cell_model_3d = new SimpleWntCellCycleModel;
        p_cell_model_3d->SetDimension(3);
        TS_ASSERT_EQUALS(p_cell_model_3d->GetDimension(), 3u);

        CellPtr p_stem_cell_3d(new Cell(p_healthy_state, p_cell_model_3d));
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

    void noTestArchiveSimpleWntCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "simple_wnt_cell_cycle.arch";

        // Set up the Wnt concentration
        double wnt_level = 1.0;
        WntConcentration<1>::Instance()->SetConstantWntValueForTesting(wnt_level);

        double random_number_test = 0;

        // Create an output archive
        {
            // Set up the Wnt concentration for testing
            WntConcentration<1>::Instance()->SetConstantWntValueForTesting(0.7);

            // Create cell-cycle model
            SimpleWntCellCycleModel* p_cell_model = new SimpleWntCellCycleModel;
            p_cell_model->SetDimension(1);
            p_cell_model->SetBirthTime(-1.0);

            // Set up the simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();

            // The number for the G1 duration is taken from
            // the first random number generated
            double g1_duration = mFirstRandomNumber;
            double end_time = g1_duration + p_cell_model->GetSG2MDuration() + 5.0;
            unsigned num_timesteps = 50;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

            // Set up associated cell
            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

            CellPtr p_stem_cell(new Cell(p_healthy_state, p_cell_model));
            p_stem_cell->SetCellProliferativeType(p_stem_type);
            p_stem_cell->InitialiseCellCycleModel();

            while (p_cell_model->GetAge() < g1_duration + p_cell_model->GetSG2MDuration()
                    - p_simulation_time->GetTimeStep()) // minus one to match birth time.
            {
                p_simulation_time->IncrementTimeOneStep();
                CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, g1_duration);
            }

            // Wnt should change this to a transit cell
            TS_ASSERT_EQUALS(p_stem_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
            TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(static_cast<SimpleWntCellCycleModel*>(p_stem_cell->GetCellCycleModel())->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            std::ofstream ofs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_oarchive output_arch(ofs);

            CellPtr const p_const_cell = p_stem_cell;
            output_arch << p_const_cell;

            TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->ReadyToDivide(), true);

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }
        {
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
            SimpleWntCellCycleModel* p_cell_model = static_cast<SimpleWntCellCycleModel*>(p_cell->GetCellCycleModel());
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);

            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetSG2MDuration(), 10.0, 1e-12);

            TS_ASSERT_DELTA(p_gen->ranf(), random_number_test, 1e-7);
            TS_ASSERT_EQUALS((static_cast<SimpleWntCellCycleModel*>(p_cell_model))->GetDimension(), 1u);

            // Tidy up
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }

        /*
         * Test the case of a radial Wnt concentration
         */

        RandomNumberGenerator::Instance()->Reseed(0);

        OutputFileHandler handler2("archive", false);
        archive_filename = handler2.GetOutputDirectoryFullPath() + "crypt_projection_cell_cycle.arch";

        // Set up the Wnt concentration
        wnt_level = 0.79;
        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(wnt_level);

        random_number_test = 0;

        // Create an output archive
        {
            // Set up the simulation time note here it needs to be done before we create a cell-cycle model.
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);

            // Create cell-cycle model
            SimpleWntCellCycleModel* p_cell_model = new SimpleWntCellCycleModel;
            p_cell_model->SetDimension(2);
            p_cell_model->SetBirthTime(-1.0);

            // Set end time for simulation
            // The number for the G1 duration is taken from
            // the first random number generated
            double g1_duration = mFirstRandomNumber;

            double end_time = g1_duration + p_cell_model->GetSG2MDuration() + 5.0;
            unsigned num_timesteps = 50;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

            boost::shared_ptr<AbstractCellMutationState> p_healthy_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->InitialiseCellCycleModel();

            // Run to division age minus one time step to match birth time
            while (p_cell_model->GetAge() < g1_duration + p_cell_model->GetSG2MDuration()
                                            - p_simulation_time->GetTimeStep())
            {
                p_simulation_time->IncrementTimeOneStep();
                CheckReadyToDivideAndPhaseIsUpdated(p_cell_model, g1_duration);
            }

            // Wnt should change this to a transit cell
            TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
            TS_ASSERT_EQUALS(p_cell->GetCellCycleModel()->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(static_cast<SimpleWntCellCycleModel*>(p_cell->GetCellCycleModel())->GetCurrentCellCyclePhase(), G_TWO_PHASE);

            std::ofstream ofs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_oarchive output_arch(ofs);

            CellPtr const p_const_cell = p_cell;
            output_arch << p_const_cell;

            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->GetCellCycleModel()->ReadyToDivide(), true);

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();

            // Tidy Up
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }
        {
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
            SimpleWntCellCycleModel* p_cell_model = static_cast<SimpleWntCellCycleModel*>(p_cell->GetCellCycleModel());
            TS_ASSERT_EQUALS(p_cell, p_cell_model->GetCell());

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

            p_simulation_time->IncrementTimeOneStep();

            TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true);
            TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);

            TS_ASSERT_DELTA(p_cell_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_cell_model->GetSG2MDuration(), 10.0, 1e-12);

            TS_ASSERT_DELTA(p_gen->ranf(), random_number_test, 1e-7);

            // Tidy up
            RandomNumberGenerator::Destroy();
            SimulationTime::Destroy();
        }

        // Tidy up
        WntConcentration<1>::Destroy();
        WntConcentration<2>::Destroy();
    }

    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with SimpleWntCellCycleModel
        SimpleWntCellCycleModel simple_wnt_cell_cycle_model;
        TS_ASSERT_EQUALS(simple_wnt_cell_cycle_model.GetIdentifier(), "SimpleWntCellCycleModel");

        out_stream simple_wnt_parameter_file = output_file_handler.OpenOutputFile("simple_wnt_results.parameters");
        simple_wnt_cell_cycle_model.OutputCellCycleModelParameters(simple_wnt_parameter_file);
        simple_wnt_parameter_file->close();

        std::string simple_wnt_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison( simple_wnt_results_dir + "simple_wnt_results.parameters", "crypt/test/data/TestCellCycleModels/simple_wnt_results.parameters").CompareFiles();
    }
};

// Member initialisation
const double TestSimpleCellCycleModelsForCrypt::mFirstRandomNumber = 1.08221;
const double TestSimpleCellCycleModelsForCrypt::mSecondRandomNumber = 3.21839;
const double TestSimpleCellCycleModelsForCrypt::mThirdRandomNumber = 3.73243;
const double TestSimpleCellCycleModelsForCrypt::mFourthRandomNumber = 2.83804;
const double TestSimpleCellCycleModelsForCrypt::mFifthRandomNumber = 1.7031;

#endif /*TESTSIMPLECELLCYCLEMODELSFORCRYPT_HPP_*/
