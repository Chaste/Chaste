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

#ifndef TESTSIMPLECELLCYCLEMODELS_HPP_
#define TESTSIMPLECELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellCycleTimesGenerator.hpp"
#include "CellLabel.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "FileComparison.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "FixedSequenceCellCycleModel.hpp"
#include "GammaG1CellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSimpleCellCycleModels : public AbstractCellBasedTestSuite
{
public:
    void TestNoCellCycleModel()
    {
        // Test constructor
        TS_ASSERT_THROWS_NOTHING(NoCellCycleModel cell_cycle_model);

        // Test methods
        NoCellCycleModel* p_model = new NoCellCycleModel;
        TS_ASSERT_DELTA(p_model->GetAverageStemCellCycleTime(), DBL_MAX, 1e-6);
        TS_ASSERT_DELTA(p_model->GetAverageTransitCellCycleTime(), DBL_MAX, 1e-6);
        TS_ASSERT_EQUALS(p_model->ReadyToDivide(), false);

        // Test the cell-cycle model works correctly with a cell
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        for (unsigned i = 0; i < 10; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
        }
    }

    void TestFixedSequenceCellCycleModelMethods()
    {
        CellCycleTimesGenerator* p_cell_cycle_times_generator = CellCycleTimesGenerator::Instance();

        // Set up the required singleton
        // Make sure we can generate this model
        TS_ASSERT_THROWS_NOTHING(FixedSequenceCellCycleModel cell_model);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        // Get a pointer to a cell cycle model of this kind
        FixedSequenceCellCycleModel* p_stem_model = new FixedSequenceCellCycleModel;

        // Test set and get method for the rate parameter
        TS_ASSERT_DELTA(p_stem_model->GetRate(), 0.5, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 2.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 2.0, 1e-10);

        p_stem_model->SetRate(10.0);
        TS_ASSERT_DELTA(p_stem_model->GetRate(), 10.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 0.1, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 0.1, 1e-10);

        p_stem_model->SetRate(0.25);
        TS_ASSERT_DELTA(p_stem_model->GetRate(), 0.25, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 4.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 4.0, 1e-10);

        p_stem_model->SetRate(0.25);
        p_cell_cycle_times_generator->SetRandomSeed(0);

        TS_ASSERT_THROWS_THIS(p_cell_cycle_times_generator->GetNextCellCycleTime(),
                              "When using FixedSequenceCellCycleModel one must call CellCycleTimesGenerator::Instance()->GenerateCellCycleTimeSequence()"
                              " before the start of the simulation.");

        p_cell_cycle_times_generator->GenerateCellCycleTimeSequence();

        TS_ASSERT_THROWS_THIS(p_stem_model->SetTransitCellG1Duration(8.0),
                              "This cell cycle model does not differentiate stem cells and transit cells, please use SetRate() instead");

        TS_ASSERT_THROWS_THIS(p_stem_model->SetStemCellG1Duration(8.0),
                              "This cell cycle model does not differentiate stem cells and transit cells, please use SetRate() instead");

        TS_ASSERT_THROWS_THIS(p_stem_model->SetRate(8.0),
                              "You cannot reset the rate after cell cycle times are created.");

        TS_ASSERT_THROWS_THIS(p_cell_cycle_times_generator->GenerateCellCycleTimeSequence(),
                              "Trying to generate the cell cycle times twice. Need to call CellCycleTimesGenerator::Destroy() first.");
        // When we set the rate parameter we also reset the TransitCellG1Duration and StemCellG1Duration such that
        // average cell cycle times are calculated correctly
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 4.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 4.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetAverageTransitCellCycleTime(), 14.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetAverageStemCellCycleTime(), 14.0, 1e-10);

        // Make a stem cell with the model
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        // Make another cell cycle model of this kind and give it to a transit cell
        FixedSequenceCellCycleModel* p_transit_model = new FixedSequenceCellCycleModel;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);

        // And finally a cell cycle model for a differentiated cell
        FixedSequenceCellCycleModel* p_diff_model = new FixedSequenceCellCycleModel;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        // The random times should be the same across platforms (We won't bother testing the distributions)
        double exponentially_generated_g1 = 3.4590;
        TS_ASSERT_DELTA(p_stem_model->GetG1Duration(), exponentially_generated_g1, 1e-4);
        TS_ASSERT_EQUALS(p_diff_model->GetG1Duration(), DBL_MAX);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(14.0, 100);
        for (unsigned i = 0; i < 100; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The actual testing of the cell cycle model
            // The numbers for the G1 durations below are taken from the first random number generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, exponentially_generated_g1);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132); // any old number
        }

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

        p_stem_model->ResetForDivision();
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), M_PHASE);

        FixedSequenceCellCycleModel* p_stem_model2 = static_cast<FixedSequenceCellCycleModel*>(p_stem_model->CreateCellCycleModel());
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);

        CellPtr p_stem_cell2(new Cell(p_healthy_state, p_stem_model2));
        p_stem_cell2->SetCellProliferativeType(p_stem_type);
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);
        CellCycleTimesGenerator::Destroy();
    }

    void TestArchiveFixedSequenceCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "FixedSequenceCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;
        double fixed_sequence_test_number = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new FixedSequenceCellCycleModel;
            static_cast<FixedSequenceCellCycleModel*>(p_model)->SetRate(13.42);

            CellCycleTimesGenerator* p_cell_cycle_times_generator = CellCycleTimesGenerator::Instance();

            p_cell_cycle_times_generator->SetRandomSeed(12u);
            p_cell_cycle_times_generator->GenerateCellCycleTimeSequence();

            p_cell_cycle_times_generator->GetNextCellCycleTime();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            fixed_sequence_test_number = p_cell_cycle_times_generator->GetNextCellCycleTime();

            delete p_model;
            SimulationTime::Destroy();

            random_number_test = RandomNumberGenerator::Instance()->ranf();
            RandomNumberGenerator::Destroy();
            CellCycleTimesGenerator::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);
            TS_ASSERT_DELTA(CellCycleTimesGenerator::Instance()->GetNextCellCycleTime(), fixed_sequence_test_number, 1e-6);

            // Check private data has been restored correctly
            TS_ASSERT_DELTA(static_cast<FixedSequenceCellCycleModel*>(p_model2)->GetRate(), 13.42, 1e-12);

            // Avoid memory leaks
            delete p_model2;
            CellCycleTimesGenerator::Destroy();
        }
    }

    void TestCellCycleTimesGeneratorSingleton()
    {
        // how to test this: first, test that it is a singleton, then a test for all methods, then a test for fixed sequence
        {
            CellCycleTimesGenerator::Instance();
            CellCycleTimesGenerator::Instance()->SetRate(298.0);
        }

        {
            CellCycleTimesGenerator::Instance();
            TS_ASSERT_DELTA(CellCycleTimesGenerator::Instance()->GetRate(), 298.0, 1e-7)
        }

        // well, looks like it's a singleton, why not test the methods next

        CellCycleTimesGenerator::Instance()->SetRate(87.0);
        TS_ASSERT_DELTA(CellCycleTimesGenerator::Instance()->GetRate(), 87.0, 1e-7);

        CellCycleTimesGenerator::Destroy();
        TS_ASSERT_DELTA(CellCycleTimesGenerator::Instance()->GetRate(), 1.0 / 2.0, 1e-7);

        CellCycleTimesGenerator::Instance()->SetRate(45.0);
        TS_ASSERT_EQUALS(CellCycleTimesGenerator::Instance()->GetRandomSeed(), 0u);
        CellCycleTimesGenerator::Instance()->SetRandomSeed(78u);
        TS_ASSERT_EQUALS(CellCycleTimesGenerator::Instance()->GetRandomSeed(), 78u);

        CellCycleTimesGenerator::Instance()->GenerateCellCycleTimeSequence();

        std::vector<double> random_sequence;

        for (unsigned index = 0u; index < 10; index++)
        {
            double next_random_number = CellCycleTimesGenerator::Instance()->GetNextCellCycleTime();
            random_sequence.push_back(next_random_number);
            if (index > 0)
            {
                assert(next_random_number != random_sequence[index - 1]);
            }
        }
        CellCycleTimesGenerator::Destroy();

        CellCycleTimesGenerator::Instance()->SetRate(45.0);
        CellCycleTimesGenerator::Instance()->SetRandomSeed(78u);

        RandomNumberGenerator::Instance()->Reseed(12u);

        CellCycleTimesGenerator::Instance()->GenerateCellCycleTimeSequence();

        RandomNumberGenerator::Instance()->Reseed(13u);

        for (unsigned index = 0u; index < 10; index++)
        {
            unsigned my_seed = 13u * index;
            RandomNumberGenerator::Instance()->Reseed(my_seed);
            double my_rate = double(index) * 0.5 + 1.0; //make sure we don't accidentally set a rate of 0, gives nasty boost errors.
            RandomNumberGenerator::Instance()->ExponentialRandomDeviate(my_rate);
            double next_cell_cycle_time = CellCycleTimesGenerator::Instance()->GetNextCellCycleTime();
            TS_ASSERT_DELTA(next_cell_cycle_time, random_sequence[index], 1e-7);
        }
        CellCycleTimesGenerator::Destroy();
    }

    void TestBernoulliTrialCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(BernoulliTrialCellCycleModel cell_model3);

        BernoulliTrialCellCycleModel* p_diff_model = new BernoulliTrialCellCycleModel;
        BernoulliTrialCellCycleModel* p_transit_model = new BernoulliTrialCellCycleModel;

        TS_ASSERT_DELTA(p_transit_model->GetDivisionProbability(), 0.1, 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetMinimumDivisionAge(), 1.0, 1e-9);

        // Change parameters for this model
        p_transit_model->SetDivisionProbability(0.5);
        p_transit_model->SetMinimumDivisionAge(0.1);
        TS_ASSERT_DELTA(p_transit_model->GetDivisionProbability(), 0.5, 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetMinimumDivisionAge(), 0.1, 1e-9);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, num_steps);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The division time below is taken from the first random number generated
            if (i < 33)
            {
                TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), false);
            }
            else
            {
                TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), true);
            }
            TS_ASSERT_EQUALS(p_diff_cell->ReadyToDivide(), false);
        }
        TS_ASSERT_DELTA(p_transit_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_DELTA(p_diff_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

        // Check that cell division correctly resets the cell cycle phase
        CellPtr p_transit_cell2 = p_transit_cell->Divide();
        BernoulliTrialCellCycleModel* p_transit_model2 = static_cast<BernoulliTrialCellCycleModel*>(p_transit_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_transit_model2->ReadyToDivide(), false);
        TS_ASSERT_DELTA(p_transit_model2->GetDivisionProbability(), 0.5, 1e-9);
        TS_ASSERT_DELTA(p_transit_model2->GetMinimumDivisionAge(), 0.1, 1e-9);

        TS_ASSERT_DELTA(p_transit_model2->GetAverageTransitCellCycleTime(), 2.0, 1e-9);
        TS_ASSERT_DELTA(p_transit_model2->GetAverageStemCellCycleTime(), 2.0, 1e-9);
    }

    void TestFixedG1GenerationalCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(FixedG1GenerationalCellCycleModel model3);

        FixedG1GenerationalCellCycleModel* p_stem_model = new FixedG1GenerationalCellCycleModel;

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100; //Per cycle
        unsigned num_cycles = 2;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            num_cycles * (p_stem_model->GetStemCellG1Duration() + p_stem_model->GetSG2MDuration()),
            num_cycles * num_steps);

        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), M_PHASE);
        TS_ASSERT_EQUALS(p_stem_model->GetGeneration(), 0u);
        TS_ASSERT_EQUALS(p_stem_model->GetMaxTransitGenerations(), 3u);
        TS_ASSERT_EQUALS(p_stem_model->CanCellTerminallyDifferentiate(), true);

        p_stem_model->SetMaxTransitGenerations(6);
        TS_ASSERT_EQUALS(p_stem_model->GetMaxTransitGenerations(), 6u);
        p_stem_model->SetMaxTransitGenerations(3);

        TS_ASSERT_EQUALS(p_stem_cell->GetCellProliferativeType(), p_stem_type);

        FixedG1GenerationalCellCycleModel* p_transit_model = new FixedG1GenerationalCellCycleModel;

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_transit_cell->GetCellProliferativeType(), p_transit_type);
        TS_ASSERT_EQUALS(p_transit_model->GetGeneration(), 0u);

        FixedG1GenerationalCellCycleModel* p_diff_model = new FixedG1GenerationalCellCycleModel;

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_diff_cell->GetCellProliferativeType(), p_diff_type);
        TS_ASSERT_EQUALS(p_diff_model->GetGeneration(), 0u);

        // First cycle
        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, p_stem_model->GetStemCellG1Duration());
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, p_transit_model->GetTransitCellG1Duration());
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 100);
        }

        TS_ASSERT_DELTA(p_stem_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_DELTA(p_diff_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

        double hepa_one_cell_birth_time = p_simulation_time->GetTime();

        FixedG1GenerationalCellCycleModel* p_hepa_one_model = new FixedG1GenerationalCellCycleModel;

        // Change G1 duration for this model
        p_hepa_one_model->SetStemCellG1Duration(8.0);
        p_hepa_one_model->SetTransitCellG1Duration(8.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->InitialiseCellCycleModel();

        // Second cycle
        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 8.0);
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge() + hepa_one_cell_birth_time, p_simulation_time->GetTime(), 1e-9);

        // Test get and set methods
        TS_ASSERT_DELTA(p_stem_model->GetSDuration(), 5.0, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetG2Duration(), 4.0, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetMDuration(), 1.0, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 2.0, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 14.0, 1e-9);

        p_stem_model->SetSDuration(7.4);
        p_stem_model->SetG2Duration(1.4);
        p_stem_model->SetMDuration(0.72);
        p_stem_model->SetTransitCellG1Duration(9.4);
        p_stem_model->SetStemCellG1Duration(9.4);

        TS_ASSERT_DELTA(p_stem_model->GetSDuration(), 7.4, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetG2Duration(), 1.4, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetMDuration(), 0.72, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 9.4, 1e-9);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 9.4, 1e-9);
    }

    void TestUniformG1GenerationalCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(UniformG1GenerationalCellCycleModel cell_model3);

        UniformG1GenerationalCellCycleModel* p_stem_model = new UniformG1GenerationalCellCycleModel;

        // Change G1 duration for this model
        p_stem_model->SetStemCellG1Duration(1.0);

        UniformG1GenerationalCellCycleModel* p_transit_model = new UniformG1GenerationalCellCycleModel;

        // Change G1 duration for this model
        p_transit_model->SetTransitCellG1Duration(1.0);

        UniformG1GenerationalCellCycleModel* p_diff_model = new UniformG1GenerationalCellCycleModel;

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0 * (p_stem_model->GetStemCellG1Duration() + p_stem_model->GetSG2MDuration()), 2 * num_steps);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three
            // random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 3.19525);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 2.18569);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132); // any old number
        }

        UniformG1GenerationalCellCycleModel* p_hepa_one_model = new UniformG1GenerationalCellCycleModel;

        // Change G1 duration for this model
        p_hepa_one_model->SetStemCellG1Duration(1.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->InitialiseCellCycleModel();

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 3.86076);
        }
    }

    void TestUniformCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(UniformCellCycleModel cell_model3);

        UniformCellCycleModel* p_stem_model = new UniformCellCycleModel;

        // Change min and max cell cycle duration for this model
        p_stem_model->SetMinCellCycleDuration(18.0);
        p_stem_model->SetMaxCellCycleDuration(20.0);

        UniformCellCycleModel* p_transit_model = new UniformCellCycleModel;

        // Use default min and max cell cycle duration for this model
        p_transit_model->SetMinCellCycleDuration(12.0);
        p_transit_model->SetMaxCellCycleDuration(14.0);

        UniformCellCycleModel* p_diff_model = new UniformCellCycleModel;

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0 * (p_stem_model->GetAverageStemCellCycleTime()), 2 * num_steps);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the cell cycle durations below are taken from the first three
            // random numbers generated
            CheckReadyToDivideIsUpdated(p_stem_model, 19.09763);
            CheckReadyToDivideIsUpdated(p_transit_model, 13.18569);
            CheckReadyToDivideIsUpdated(p_diff_model, 132); // any old number
        }

        // Check with a mutation.
        UniformCellCycleModel* p_hepa_one_model = new UniformCellCycleModel;

        // Use default Min and Max cell cycle duration for this model
        p_hepa_one_model->SetMinCellCycleDuration(14.0);
        p_hepa_one_model->SetMaxCellCycleDuration(16.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->InitialiseCellCycleModel();

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideIsUpdated(p_hepa_one_model, 15.43038);
        }
    }

    void TestGammaG1CellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(GammaG1CellCycleModel cell_model);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        GammaG1CellCycleModel* p_stem_model = new GammaG1CellCycleModel;
        p_stem_model->SetShape(3.517);
        p_stem_model->SetScale(2.986);

        TS_ASSERT_DELTA(p_stem_model->GetShape(), 3.517, 1e-4);
        TS_ASSERT_DELTA(p_stem_model->GetScale(), 2.986, 1e-4);

        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        GammaG1CellCycleModel* p_transit_model = new GammaG1CellCycleModel;
        p_transit_model->SetShape(3.5);
        p_transit_model->SetScale(2.9);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        GammaG1CellCycleModel* p_diff_model = new GammaG1CellCycleModel;
        p_diff_model->SetShape(3.5);
        p_diff_model->SetScale(2.9);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_stem_model->GetG1Duration(), 3.6104, 1e-4);
        TS_ASSERT_DELTA(p_transit_model->GetG1Duration(), 3.8511, 1e-4);
        TS_ASSERT_EQUALS(p_diff_model->GetG1Duration(), DBL_MAX);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(14.0, 100);
        for (unsigned i = 0; i < 100; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 3.61046);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 3.8511);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132); // any old number
        }

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

        p_stem_model->ResetForDivision();
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), M_PHASE);

        GammaG1CellCycleModel* p_stem_model2 = static_cast<GammaG1CellCycleModel*>(p_stem_model->CreateCellCycleModel());
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);

        CellPtr p_stem_cell2(new Cell(p_healthy_state, p_stem_model2));
        p_stem_cell2->SetCellProliferativeType(p_stem_type);
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);
    }

    void TestExponentialG1GenerationalCellCycleModel()
    {
        // Make sure we can generate this model
        TS_ASSERT_THROWS_NOTHING(ExponentialG1GenerationalCellCycleModel cell_model);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        // Get a pointer to a cell cycle model of this kind
        ExponentialG1GenerationalCellCycleModel* p_stem_model = new ExponentialG1GenerationalCellCycleModel;

        // Test set and get method for the rate parameter
        TS_ASSERT_DELTA(p_stem_model->GetRate(), 0.5, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 14.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 2.0, 1e-10);

        p_stem_model->SetStemCellG1Duration(8.0);
        TS_ASSERT_DELTA(p_stem_model->GetRate(), 0.125, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 8.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 2.0, 1e-10);

        p_stem_model->SetTransitCellG1Duration(0.1);
        TS_ASSERT_DELTA(p_stem_model->GetRate(), 10.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 8.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 0.1, 1e-10);

        p_stem_model->SetRate(0.25);
        TS_ASSERT_DELTA(p_stem_model->GetRate(), 0.25, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 4.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 4.0, 1e-10);

        // When we set the rate parameter we also reset the TransitCellG1Duration and StemCellG1Duration such that
        // average cell cycle times are calculated correctly
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 4.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 4.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetAverageTransitCellCycleTime(), 14.0, 1e-10);
        TS_ASSERT_DELTA(p_stem_model->GetAverageStemCellCycleTime(), 14.0, 1e-10);

        // Make a stem cell with the model
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        // Make another cell cycle model of this kind and give it to a transit cell
        ExponentialG1GenerationalCellCycleModel* p_transit_model = new ExponentialG1GenerationalCellCycleModel;
        p_transit_model->SetRate(1.0);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        // And finally a cell cycle model for a differentiated cell
        ExponentialG1GenerationalCellCycleModel* p_diff_model = new ExponentialG1GenerationalCellCycleModel;
        p_diff_model->SetRate(200.0);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        // The random times should be the same across platforms (We won't bother testing the distributions)
        double exponentially_generated_g1 = 3.4590;
        TS_ASSERT_DELTA(p_stem_model->GetG1Duration(), exponentially_generated_g1, 1e-4);
        double second_exponentially_generated_g1 = 1.3666;
        TS_ASSERT_DELTA(p_transit_model->GetG1Duration(), second_exponentially_generated_g1, 1e-4);
        TS_ASSERT_EQUALS(p_diff_model->GetG1Duration(), DBL_MAX);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(14.0, 100);
        for (unsigned i = 0; i < 100; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The actual testing of the cell cycle model
            // The numbers for the G1 durations below are taken from the first two random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, exponentially_generated_g1);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, second_exponentially_generated_g1);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132); // any old number
        }

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

        p_stem_model->ResetForDivision();
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), M_PHASE);

        ExponentialG1GenerationalCellCycleModel* p_stem_model2 = static_cast<ExponentialG1GenerationalCellCycleModel*>(p_stem_model->CreateCellCycleModel());
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);

        CellPtr p_stem_cell2(new Cell(p_healthy_state, p_stem_model2));
        p_stem_cell2->SetCellProliferativeType(p_stem_type);
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);
    }

    void TestSimpleOxygenBasedCellCycleModel()
    {
        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
        p_model->SetDimension(2);

        p_model->SetStemCellG1Duration(8.0);
        p_model->SetTransitCellG1Duration(8.0);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

        // Set up oxygen_concentration
        double lo_oxygen_concentration = 0.0;
        double hi_oxygen_concentration = 1.0;
        p_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_cell->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);

        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        p_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);
        p_simulation_time->IncrementTimeOneStep(); // t=3.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(4.0 * 18.0, num_steps);

        TS_ASSERT_THROWS_NOTHING(SimpleOxygenBasedCellCycleModel model);

        // Create cell-cycle models and cells
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model = new SimpleOxygenBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);

        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);
        p_hepa_one_cell->InitialiseCellCycleModel();

        SimpleOxygenBasedCellCycleModel* p_diff_model = new SimpleOxygenBasedCellCycleModel;
        p_diff_model->SetDimension(2);

        // Coverage
        TS_ASSERT_DELTA(p_diff_model->GetCriticalHypoxicDuration(), 2.0, 1e-6);
        p_diff_model->SetCriticalHypoxicDuration(0.5);
        TS_ASSERT_DELTA(p_diff_model->GetCriticalHypoxicDuration(), 0.5, 1e-6);

        CellPtr p_diff_cell(new Cell(p_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();
        p_diff_cell->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);

        // Check that the cell cycle phase and ready to divide
        // are updated correctly
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(), M_PHASE);

        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(), G_ZERO_PHASE);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_hepa_one_model->GetG1Duration());
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), true);

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_hepa_one_cell->ReadyToDivide(), true);
        CellPtr p_hepa_one_cell2 = p_hepa_one_cell->Divide();
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast<SimpleOxygenBasedCellCycleModel*>(p_hepa_one_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0 * p_hepa_one_model2->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell-cycle model
        SimpleOxygenBasedCellCycleModel* p_cell_model = new SimpleOxygenBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        CellPtr p_apoptotic_cell(new Cell(p_state, p_cell_model));
        p_apoptotic_cell->SetCellProliferativeType(p_stem_type);

        // Set up oxygen_concentration
        p_apoptotic_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);

        // Force the cell to be apoptotic
        for (unsigned i = 0; i < num_steps; i++)
        {
            TS_ASSERT(!(p_apoptotic_cell->HasCellProperty<ApoptoticCellProperty>())
                      || p_simulation_time->GetTime() >= p_cell_model->GetCriticalHypoxicDuration());
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            p_apoptotic_cell->ReadyToDivide();
        }

        // Test that the cell is updated to be apoptotic
        TS_ASSERT_EQUALS(p_apoptotic_cell->HasCellProperty<ApoptoticCellProperty>(), true);
        TS_ASSERT_EQUALS(p_cell_model->GetCurrentHypoxicDuration(), 2.04);
    }

    void TestContactInhibitionCellCycleModel()
    {
        // Check that mQuiescentVolumeFraction and mEquilibriumVolume are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
        p_model->SetDimension(2);

        p_model->SetStemCellG1Duration(8.0);
        p_model->SetTransitCellG1Duration(8.0);

        TS_ASSERT_THROWS_THIS(p_model->UpdateCellCyclePhase(),
                              "The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");

        p_model->SetQuiescentVolumeFraction(0.5);
        p_model->SetEquilibriumVolume(1.0);

        // Set the birth time such that at t=0, the cell has just entered G1 phase
        p_model->SetBirthTime(-1.0);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

        double lo_volume = 0.0;
        double hi_volume = 1.0;
        p_cell->GetCellData()->SetItem("volume", lo_volume);

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentOnsetTime(), 0.0, 1e-12);

        p_cell->GetCellData()->SetItem("volume", hi_volume);
        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentOnsetTime(), 2.0, 1e-12);

        p_cell->GetCellData()->SetItem("volume", lo_volume);
        p_simulation_time->IncrementTimeOneStep(); // t=3.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentOnsetTime(), 2.0, 1e-12);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0 * 24.0, num_steps);

        // Create cell-cycle models and cells
        ContactInhibitionCellCycleModel* p_hepa_one_model = new ContactInhibitionCellCycleModel;
        p_hepa_one_model->SetDimension(2);
        p_hepa_one_model->SetBirthTime(0.0);
        p_hepa_one_model->SetQuiescentVolumeFraction(0.5);
        p_hepa_one_model->SetEquilibriumVolume(1.0);

        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);

        // Set up cell volume
        p_hepa_one_cell->GetCellData()->SetItem("volume", hi_volume);
        p_hepa_one_cell->InitialiseCellCycleModel();

        ContactInhibitionCellCycleModel* p_diff_model = new ContactInhibitionCellCycleModel;
        p_diff_model->SetDimension(2);
        p_diff_model->SetQuiescentVolumeFraction(0.5);
        p_diff_model->SetEquilibriumVolume(1.0);

        CellPtr p_diff_cell(new Cell(p_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->GetCellData()->SetItem("volume", hi_volume);
        p_diff_cell->InitialiseCellCycleModel();

        // Check that the cell cycle phase and ready to divide are updated correctly
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(), M_PHASE);

        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(), G_ZERO_PHASE);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_hepa_one_model->GetG1Duration());
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_hepa_one_cell->ReadyToDivide(), true);

        CellPtr p_hepa_one_cell2 = p_hepa_one_cell->Divide();
        ContactInhibitionCellCycleModel* p_hepa_one_model2 = static_cast<ContactInhibitionCellCycleModel*>(p_hepa_one_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(), M_PHASE);
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);
    }

    void TestStochasticOxygenBasedCellCycleModel()
    {
        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel;
        p_model->SetDimension(2);

        p_model->SetStemCellG1Duration(8.0);
        p_model->SetTransitCellG1Duration(8.0);

        // Coverage
        TS_ASSERT_DELTA(p_model->GetHypoxicConcentration(), 0.4, 1e-6);
        TS_ASSERT_DELTA(p_model->GetQuiescentConcentration(), 1.0, 1e-6);
        TS_ASSERT_DELTA(p_model->GetCriticalHypoxicDuration(), 2.0, 1e-6);

        p_model->SetHypoxicConcentration(0.5);
        p_model->SetQuiescentConcentration(0.5);
        p_model->SetCriticalHypoxicDuration(3.0);

        TS_ASSERT_DELTA(p_model->GetHypoxicConcentration(), 0.5, 1e-6);
        TS_ASSERT_DELTA(p_model->GetQuiescentConcentration(), 0.5, 1e-6);
        TS_ASSERT_DELTA(p_model->GetCriticalHypoxicDuration(), 3.0, 1e-6);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        // Set up oxygen_concentration
        double lo_oxygen_concentration = 0.0;
        double hi_oxygen_concentration = 1.0;

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);
        p_cell->InitialiseCellCycleModel();

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_cell->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);
        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        p_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);
        p_simulation_time->IncrementTimeOneStep(); // t=3.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        // Set up simulation time
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();

        unsigned num_steps = 100;
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(4.0 * 18.0, num_steps);

        TS_ASSERT_THROWS_NOTHING(StochasticOxygenBasedCellCycleModel model);

        // Create cell-cycle model
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model = new StochasticOxygenBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);

        StochasticOxygenBasedCellCycleModel* p_diff_model = new StochasticOxygenBasedCellCycleModel;
        p_diff_model->SetDimension(2);

        // Create cell
        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);
        p_hepa_one_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);
        p_diff_cell->InitialiseCellCycleModel();

        // Check that the cell cycle phase and ready to divide
        // are updated correctly
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(), M_PHASE);

        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(), G_ZERO_PHASE);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_hepa_one_model->GetG1Duration(), p_hepa_one_model->GetG2Duration());
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), true);

        // Coverage
        TS_ASSERT_EQUALS(p_hepa_one_cell->ReadyToDivide(), true);
        p_hepa_one_cell->Divide();

        // Check that cell division correctly resets the cell cycle phase
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast<StochasticOxygenBasedCellCycleModel*>(p_hepa_one_model->CreateCellCycleModel());
        CellPtr p_hepa_one_cell2(new Cell(p_state, p_hepa_one_model2));
        p_hepa_one_cell2->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell2->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);
        p_hepa_one_cell2->InitialiseCellCycleModel();
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0 * p_hepa_one_model2->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell-cycle model
        StochasticOxygenBasedCellCycleModel* p_cell_model = new StochasticOxygenBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        CellPtr p_apoptotic_cell(new Cell(p_state, p_cell_model));
        p_apoptotic_cell->SetCellProliferativeType(p_stem_type);
        p_apoptotic_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);
        p_apoptotic_cell->InitialiseCellCycleModel();

        // Force the cell to be apoptotic
        for (unsigned i = 0; i < num_steps; i++)
        {
            TS_ASSERT(!(p_apoptotic_cell->HasCellProperty<ApoptoticCellProperty>())
                      || p_simulation_time->GetTime() >= p_cell_model->GetCriticalHypoxicDuration());
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            p_apoptotic_cell->ReadyToDivide();
        }

        // Test that the cell is updated to be apoptotic
        TS_ASSERT_EQUALS(p_apoptotic_cell->HasCellProperty<ApoptoticCellProperty>(), true);
        TS_ASSERT_EQUALS(p_cell_model->GetCurrentHypoxicDuration(), 2.04);

        StochasticOxygenBasedCellCycleModel* p_cell_model2 = new StochasticOxygenBasedCellCycleModel;
        p_cell_model2->SetDimension(2);

        // Coverage
        p_cell_model2->SetMinimumGapDuration(1e20);
        TS_ASSERT_DELTA(p_cell_model2->GetMinimumGapDuration(), 1e20, 1e-4);

        CellPtr p_cell2(new Cell(p_state, p_cell_model2));
        p_cell2->SetCellProliferativeType(p_stem_type);
        p_cell2->InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model2->GetG2Duration(), 1e20, 1e-4);
    }

    void TestArchiveNoCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "NoCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new NoCellCycleModel;

            p_model->SetDimension(2);
            p_model->SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            // Check private data has been restored correctly
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.0, 1e-12);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);
            TS_ASSERT_EQUALS(p_model2->ReadyToDivide(), false);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveBernoulliTrialCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BernoulliTrialCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new BernoulliTrialCellCycleModel;

            p_model->SetDimension(2);
            p_model->SetBirthTime(-1.0);
            static_cast<BernoulliTrialCellCycleModel*>(p_model)->SetDivisionProbability(0.5);
            static_cast<BernoulliTrialCellCycleModel*>(p_model)->SetMinimumDivisionAge(0.1);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            // Check private data has been restored correctly
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.0, 1e-12);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);
            TS_ASSERT_DELTA(static_cast<BernoulliTrialCellCycleModel*>(p_model2)->GetDivisionProbability(), 0.5, 1e-9);
            TS_ASSERT_DELTA(static_cast<BernoulliTrialCellCycleModel*>(p_model2)->GetMinimumDivisionAge(), 0.1, 1e-9);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveFixedG1GenerationalCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "FixedG1GenerationalCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new FixedG1GenerationalCellCycleModel;

            p_model->SetDimension(2);
            p_model->SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            // Check private data has been restored correctly
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.0, 1e-12);
            TS_ASSERT_EQUALS(static_cast<AbstractPhaseBasedCellCycleModel*>(p_model2)->GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveUniformG1GenerationalCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "UniformG1GenerationalCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new UniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);
            static_cast<UniformG1GenerationalCellCycleModel*>(p_model)->SetTransitCellG1Duration(1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();

            random_number_test = RandomNumberGenerator::Instance()->ranf();
            RandomNumberGenerator::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveUniformCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "UniformCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new UniformCellCycleModel;
            p_model->SetDimension(2);
            dynamic_cast<UniformCellCycleModel*>(p_model)->SetMinCellCycleDuration(1.0);
            dynamic_cast<UniformCellCycleModel*>(p_model)->SetMaxCellCycleDuration(2.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();

            random_number_test = RandomNumberGenerator::Instance()->ranf();
            RandomNumberGenerator::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);
            TS_ASSERT_DELTA(dynamic_cast<UniformCellCycleModel*>(p_model2)->GetMinCellCycleDuration(), 1.0, 1e-5);
            TS_ASSERT_DELTA(dynamic_cast<UniformCellCycleModel*>(p_model2)->GetMaxCellCycleDuration(), 2.0, 1e-5);

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveGammaG1CellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "GammaG1CellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new GammaG1CellCycleModel;
            static_cast<GammaG1CellCycleModel*>(p_model)->SetShape(2.45);
            static_cast<GammaG1CellCycleModel*>(p_model)->SetScale(13.42);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();

            random_number_test = RandomNumberGenerator::Instance()->ranf();
            RandomNumberGenerator::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Check private data has been restored correctly
            TS_ASSERT_DELTA(static_cast<GammaG1CellCycleModel*>(p_model2)->GetShape(), 2.45, 1e-12);
            TS_ASSERT_DELTA(static_cast<GammaG1CellCycleModel*>(p_model2)->GetScale(), 13.42, 1e-12);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveExponentialG1GenerationalCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ExponentialG1GenerationalCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new ExponentialG1GenerationalCellCycleModel;
            static_cast<ExponentialG1GenerationalCellCycleModel*>(p_model)->SetRate(13.42);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();

            random_number_test = RandomNumberGenerator::Instance()->ranf();
            RandomNumberGenerator::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Check private data has been restored correctly
            TS_ASSERT_DELTA(static_cast<ExponentialG1GenerationalCellCycleModel*>(p_model2)->GetRate(), 13.42, 1e-12);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveSimpleOxygenBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SimpleOxygenBasedCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(3);

            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetHypoxicConcentration(0.8);
            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetQuiescentConcentration(0.7);
            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetCriticalHypoxicDuration(2.5);
            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetCurrentHypoxiaOnsetTime(3.1);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetHypoxicConcentration(), 0.8, 1e-6);
            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetQuiescentConcentration(), 0.7, 1e-6);
            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetCriticalHypoxicDuration(), 2.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetCurrentHypoxiaOnsetTime(), 3.1, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveStochasticOxygenBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StochasticOxygenBasedCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new StochasticOxygenBasedCellCycleModel;
            p_model->SetDimension(3);

            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetHypoxicConcentration(0.8);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetQuiescentConcentration(0.7);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetCriticalHypoxicDuration(2.5);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetCurrentHypoxiaOnsetTime(3.1);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->GenerateStochasticG2Duration();

            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->GetG2Duration(), 3.0822, 1e-4);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();

            random_number_test = RandomNumberGenerator::Instance()->ranf();
            RandomNumberGenerator::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetG2Duration(), 3.0822, 1e-4);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetHypoxicConcentration(), 0.8, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetQuiescentConcentration(), 0.7, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetCriticalHypoxicDuration(), 2.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetCurrentHypoxiaOnsetTime(), 3.1, 1e-6);

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveContactInhibitionCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ContactInhibitionCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new ContactInhibitionCellCycleModel();

            p_model->SetDimension(1);
            p_model->SetBirthTime(-1.5);

            static_cast<ContactInhibitionCellCycleModel*>(p_model)->SetQuiescentVolumeFraction(0.5);
            static_cast<ContactInhibitionCellCycleModel*>(p_model)->SetEquilibriumVolume(1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_EQUALS(p_model2->GetDimension(), 1u);
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.5, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.5, 1e-12);

            TS_ASSERT_EQUALS(static_cast<AbstractPhaseBasedCellCycleModel*>(p_model2)->GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_DELTA(static_cast<ContactInhibitionCellCycleModel*>(p_model2)->GetQuiescentVolumeFraction(), 0.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<ContactInhibitionCellCycleModel*>(p_model2)->GetEquilibriumVolume(), 1.0, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with NoCellCycleModel
        NoCellCycleModel no_model;
        TS_ASSERT_EQUALS(no_model.GetIdentifier(), "NoCellCycleModel");

        out_stream no_model_parameter_file = output_file_handler.OpenOutputFile("no_model_results.parameters");
        no_model.OutputCellCycleModelParameters(no_model_parameter_file);
        no_model_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("no_model_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/no_model_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with BernoulliTrialCellCycleModel
        BernoulliTrialCellCycleModel random_division_cell_cycle_model;
        TS_ASSERT_EQUALS(random_division_cell_cycle_model.GetIdentifier(), "BernoulliTrialCellCycleModel");

        out_stream random_division_parameter_file = output_file_handler.OpenOutputFile("random_division_results.parameters");
        random_division_cell_cycle_model.OutputCellCycleModelParameters(random_division_parameter_file);
        random_division_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("random_division_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/random_division_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with FixedG1GenerationalCellCycleModel
        FixedG1GenerationalCellCycleModel fixed_duration_generation_based_cell_cycle_model;
        TS_ASSERT_EQUALS(fixed_duration_generation_based_cell_cycle_model.GetIdentifier(), "FixedG1GenerationalCellCycleModel");

        out_stream fixed_duration_generation_based_parameter_file = output_file_handler.OpenOutputFile("fixed_duration_generation_based_results.parameters");
        fixed_duration_generation_based_cell_cycle_model.OutputCellCycleModelParameters(fixed_duration_generation_based_parameter_file);
        fixed_duration_generation_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("fixed_duration_generation_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/fixed_duration_generation_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with UniformG1GenerationalCellCycleModel
        UniformG1GenerationalCellCycleModel uniform_distributed_generation_based_cell_cycle_model;
        TS_ASSERT_EQUALS(uniform_distributed_generation_based_cell_cycle_model.GetIdentifier(), "UniformG1GenerationalCellCycleModel");

        out_stream uniform_distributed_generation_based_parameter_file = output_file_handler.OpenOutputFile("uniform_distributed_generation_based_results.parameters");
        uniform_distributed_generation_based_cell_cycle_model.OutputCellCycleModelParameters(uniform_distributed_generation_based_parameter_file);
        uniform_distributed_generation_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("uniform_distributed_generation_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/uniform_distributed_generation_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with UniformCellCycleModel
        UniformCellCycleModel uniform_distributed_cell_cycle_model;
        TS_ASSERT_EQUALS(uniform_distributed_cell_cycle_model.GetIdentifier(), "UniformCellCycleModel");

        out_stream uniform_distributed_parameter_file = output_file_handler.OpenOutputFile("uniform_distributed_results.parameters");
        uniform_distributed_cell_cycle_model.OutputCellCycleModelParameters(uniform_distributed_parameter_file);
        uniform_distributed_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("uniform_distributed_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/uniform_distributed_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with SimpleOxygenBasedCellCycleModel
        SimpleOxygenBasedCellCycleModel simple_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(simple_oxygen_based_cell_cycle_model.GetIdentifier(), "SimpleOxygenBasedCellCycleModel");

        out_stream simple_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("simple_oxygen_based_results.parameters");
        simple_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(simple_oxygen_based_parameter_file);
        simple_oxygen_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("simple_oxygen_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/simple_oxygen_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with StochasticOxygenBasedCellCycleModel
        StochasticOxygenBasedCellCycleModel stochastic_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(stochastic_oxygen_based_cell_cycle_model.GetIdentifier(), "StochasticOxygenBasedCellCycleModel");

        out_stream stochastic_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("stochastic_oxygen_based_results.parameters");
        stochastic_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(stochastic_oxygen_based_parameter_file);
        stochastic_oxygen_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("stochastic_oxygen_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/stochastic_oxygen_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with ContactInhibitionCellCycleModel
        ContactInhibitionCellCycleModel contact_inhibition_cell_cycle_model;
        contact_inhibition_cell_cycle_model.SetQuiescentVolumeFraction(0.5);
        contact_inhibition_cell_cycle_model.SetEquilibriumVolume(1.0);
        TS_ASSERT_EQUALS(contact_inhibition_cell_cycle_model.GetIdentifier(), "ContactInhibitionCellCycleModel");

        out_stream contact_inhibition_parameter_file = output_file_handler.OpenOutputFile("contact_inhibition_results.parameters");
        contact_inhibition_cell_cycle_model.OutputCellCycleModelParameters(contact_inhibition_parameter_file);
        contact_inhibition_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("contact_inhibition_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/contact_inhibition_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with GammaG1CellCycleModel
        GammaG1CellCycleModel gamma_distributed_cell_cycle_model;
        gamma_distributed_cell_cycle_model.SetShape(0.85);
        gamma_distributed_cell_cycle_model.SetScale(1.23);
        TS_ASSERT_EQUALS(gamma_distributed_cell_cycle_model.GetIdentifier(), "GammaG1CellCycleModel");

        out_stream gamma_distributed_parameter_file = output_file_handler.OpenOutputFile("gamma_distributed_results.parameters");
        gamma_distributed_cell_cycle_model.OutputCellCycleModelParameters(gamma_distributed_parameter_file);
        gamma_distributed_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("gamma_distributed_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/gamma_distributed_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with ExponentialG1GenerationalCellCycleModel
        ExponentialG1GenerationalCellCycleModel exponential_cell_cycle_model;
        exponential_cell_cycle_model.SetRate(1.23);
        TS_ASSERT_EQUALS(exponential_cell_cycle_model.GetIdentifier(), "ExponentialG1GenerationalCellCycleModel");

        out_stream exponential_parameter_file = output_file_handler.OpenOutputFile("exponential_results.parameters");
        exponential_cell_cycle_model.OutputCellCycleModelParameters(exponential_parameter_file);
        exponential_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("exponential_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/exponential_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with FixedSequenceCellCycleModel
        FixedSequenceCellCycleModel fixed_sequence_cell_cycle_model;
        fixed_sequence_cell_cycle_model.SetRate(1.23);
        TS_ASSERT_EQUALS(fixed_sequence_cell_cycle_model.GetIdentifier(), "FixedSequenceCellCycleModel");

        out_stream fixed_sequence_parameter_file = output_file_handler.OpenOutputFile("fixed_sequence_results.parameters");
        fixed_sequence_cell_cycle_model.OutputCellCycleModelParameters(fixed_sequence_parameter_file);
        fixed_sequence_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("fixed_sequence_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/fixed_sequence_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        CellCycleTimesGenerator::Destroy();
    }
};

#endif /*TESTSIMPLECELLCYCLEMODELS_HPP_*/
