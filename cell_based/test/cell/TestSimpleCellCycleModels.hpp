/*

Copyright (c) 2005-2016, University of Oxford.
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

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "FixedDurationGenerationBasedPhaseBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedPhaseBasedCellCycleModel.hpp"
#include "GammaDistributedStochasticDurationPhaseBasedCellCycleModel.hpp"
#include "ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel.hpp"
#include "StochasticDurationPhaseBasedCellCycleModel.hpp"
#include "SimpleOxygenBasedPhaseBasedCellCycleModel.hpp"
#include "StochasticOxygenBasedPhaseBasedCellCycleModel.hpp"
#include "ContactInhibitionPhaseBasedCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSimplePhaseBasedCellCycleModels : public AbstractCellBasedTestSuite
{
public:

    void TestFixedDurationGenerationBasedPhaseBasedCellCycleModel() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(FixedDurationGenerationBasedPhaseBasedCellCycleModel model3);

        FixedDurationGenerationBasedPhaseBasedCellCycleModel* p_stem_model = new FixedDurationGenerationBasedPhaseBasedCellCycleModel;

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
            num_cycles*(p_stem_model->GetStemCellG1Duration() + p_stem_model->GetSG2MDuration()),
            num_cycles*num_steps);

        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(),M_PHASE);
        TS_ASSERT_EQUALS(p_stem_model->GetGeneration(), 0u);
        TS_ASSERT_EQUALS(p_stem_model->GetMaxTransitGenerations(), 3u);
        TS_ASSERT_EQUALS(p_stem_model->CanCellTerminallyDifferentiate(), true);

        p_stem_model->SetMaxTransitGenerations(6);
        TS_ASSERT_EQUALS(p_stem_model->GetMaxTransitGenerations(), 6u);
        p_stem_model->SetMaxTransitGenerations(3);

        TS_ASSERT_EQUALS(p_stem_cell->GetCellProliferativeType(), p_stem_type);

        FixedDurationGenerationBasedPhaseBasedCellCycleModel* p_transit_model = new FixedDurationGenerationBasedPhaseBasedCellCycleModel;

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_transit_cell->GetCellProliferativeType(), p_transit_type);
        TS_ASSERT_EQUALS(p_transit_model->GetGeneration(), 0u);

        FixedDurationGenerationBasedPhaseBasedCellCycleModel* p_diff_model = new FixedDurationGenerationBasedPhaseBasedCellCycleModel;

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_diff_cell->GetCellProliferativeType(), p_diff_type);
        TS_ASSERT_EQUALS(p_diff_model->GetGeneration(), 0u);

        //First cycle
        for (unsigned i=0; i<num_steps; i++)
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

        FixedDurationGenerationBasedPhaseBasedCellCycleModel* p_hepa_one_model = new FixedDurationGenerationBasedPhaseBasedCellCycleModel;

        // Change G1 Duration for this model
        p_hepa_one_model->SetStemCellG1Duration(8.0);
        p_hepa_one_model->SetTransitCellG1Duration(8.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->InitialiseCellCycleModel();

        //Second cycle
        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 8.0);
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge() + hepa_one_cell_birth_time, p_simulation_time->GetTime(), 1e-9);

        // Test Get and Set Methods
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

    void TestStochasticDurationGenerationBasedPhaseBasedCellCycleModel() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(StochasticDurationGenerationBasedPhaseBasedCellCycleModel cell_model3);

        StochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_stem_model = new StochasticDurationGenerationBasedPhaseBasedCellCycleModel;

        // Change G1 Duration for this model
        p_stem_model->SetStemCellG1Duration(1.0);

        StochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_transit_model = new StochasticDurationGenerationBasedPhaseBasedCellCycleModel;

        // Change G1 Duration for this model
        p_transit_model->SetTransitCellG1Duration(1.0);

        StochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_diff_model = new StochasticDurationGenerationBasedPhaseBasedCellCycleModel;

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
                        4.0*(p_stem_model->GetStemCellG1Duration() + p_stem_model->GetSG2MDuration()), 2*num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three
            // random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 3.19525);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 2.18569);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        StochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_hepa_one_model = new StochasticDurationGenerationBasedPhaseBasedCellCycleModel;

        // Change G1 Duration for this model
        p_hepa_one_model->SetStemCellG1Duration(1.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i< num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 3.86076);
        }
    }

    void TestStochasticDurationPhaseBasedCellCycleModel() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(StochasticDurationPhaseBasedCellCycleModel cell_model3);

        StochasticDurationPhaseBasedCellCycleModel* p_stem_model = new StochasticDurationPhaseBasedCellCycleModel;

        // Change G1 duration for this model
        p_stem_model->SetStemCellG1Duration(8.0);

        StochasticDurationPhaseBasedCellCycleModel* p_transit_model = new StochasticDurationPhaseBasedCellCycleModel;

        // Change G1 duration for this model
        p_stem_model->SetTransitCellG1Duration(8.0);

        StochasticDurationPhaseBasedCellCycleModel* p_diff_model = new StochasticDurationPhaseBasedCellCycleModel;

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
                        4.0*(p_stem_model->GetStemCellG1Duration() + p_stem_model->GetSG2MDuration()), 2*num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three
            // random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 9.09763);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 3.18569);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        StochasticDurationPhaseBasedCellCycleModel* p_hepa_one_model = new StochasticDurationPhaseBasedCellCycleModel;
        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i< num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 15.4304);
        }
    }

    void TestGammaDistributedStochasticDurationPhaseBasedCellCycleModel() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(GammaDistributedStochasticDurationPhaseBasedCellCycleModel cell_model);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        GammaDistributedStochasticDurationPhaseBasedCellCycleModel* p_stem_model = new GammaDistributedStochasticDurationPhaseBasedCellCycleModel;
        p_stem_model->SetShape(3.517);
        p_stem_model->SetScale(2.986);

        TS_ASSERT_DELTA(p_stem_model->GetShape(), 3.517, 1e-4);
        TS_ASSERT_DELTA(p_stem_model->GetScale(), 2.986, 1e-4);

        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        GammaDistributedStochasticDurationPhaseBasedCellCycleModel* p_transit_model = new GammaDistributedStochasticDurationPhaseBasedCellCycleModel;
        p_transit_model->SetShape(3.5);
        p_transit_model->SetScale(2.9);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        GammaDistributedStochasticDurationPhaseBasedCellCycleModel* p_diff_model = new GammaDistributedStochasticDurationPhaseBasedCellCycleModel;
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
        for (unsigned i=0; i<100; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 3.61046);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 3.8511);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

        p_stem_model->ResetForDivision();
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), M_PHASE);

        GammaDistributedStochasticDurationPhaseBasedCellCycleModel* p_stem_model2 = static_cast <GammaDistributedStochasticDurationPhaseBasedCellCycleModel*> (p_stem_model->CreateCellCycleModel());
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);

        CellPtr p_stem_cell2(new Cell(p_healthy_state, p_stem_model2));
        p_stem_cell2->SetCellProliferativeType(p_stem_type);
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);
    }

    void TestExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel() throw(Exception)
    {
        // Make sure we can generate this model
        TS_ASSERT_THROWS_NOTHING(ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel cell_model);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        // Get a pointer to a cell cycle model of this kind
        ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_stem_model =
                new ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel;

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
        ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_transit_model = new ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel;
        p_transit_model->SetRate(1.0);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        // And finally a cell cycle model for a differentiated cell
        ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_diff_model = new ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel;
        p_diff_model->SetRate(200.0);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        // The random times should be the same across platforms (We won't bother testing the distributions)
        TS_ASSERT_DELTA(p_stem_model->GetG1Duration(), 3.1834, 1e-4);
        TS_ASSERT_DELTA(p_transit_model->GetG1Duration(),0.8985, 1e-4);
        TS_ASSERT_EQUALS(p_diff_model->GetG1Duration(), DBL_MAX);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(14.0, 100);
        for (unsigned i=0; i<100; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The actual testing of the cell cycle model
            // The numbers for the G1 durations below are taken from the first three random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model,3.1834);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 0.8985);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), G_TWO_PHASE);

        p_stem_model->ResetForDivision();
        TS_ASSERT_EQUALS(p_stem_model->GetCurrentCellCyclePhase(), M_PHASE);

        ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel* p_stem_model2 =
                        static_cast <ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel*> (p_stem_model->CreateCellCycleModel());
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);

        CellPtr p_stem_cell2(new Cell(p_healthy_state, p_stem_model2));
        p_stem_cell2->SetCellProliferativeType(p_stem_type);
        TS_ASSERT_EQUALS(p_stem_model2->GetCurrentCellCyclePhase(), M_PHASE);
    }

    void TestSimpleOxygenBasedPhaseBasedCellCycleModel() throw(Exception)
    {
        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        SimpleOxygenBasedPhaseBasedCellCycleModel* p_model = new SimpleOxygenBasedPhaseBasedCellCycleModel;
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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(4.0*18.0, num_steps);

        TS_ASSERT_THROWS_NOTHING(SimpleOxygenBasedPhaseBasedCellCycleModel model);

        // Create cell-cycle models and cells
        SimpleOxygenBasedPhaseBasedCellCycleModel* p_hepa_one_model = new SimpleOxygenBasedPhaseBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);

        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->GetCellData()->SetItem("oxygen", hi_oxygen_concentration);
        p_hepa_one_cell->InitialiseCellCycleModel();

        SimpleOxygenBasedPhaseBasedCellCycleModel* p_diff_model = new SimpleOxygenBasedPhaseBasedCellCycleModel;
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
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(),M_PHASE);

        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(),G_ZERO_PHASE);

        for (unsigned i=0; i<num_steps; i++)
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
        SimpleOxygenBasedPhaseBasedCellCycleModel* p_hepa_one_model2 = static_cast <SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_hepa_one_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_hepa_one_model2->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell-cycle model
        SimpleOxygenBasedPhaseBasedCellCycleModel* p_cell_model = new SimpleOxygenBasedPhaseBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        CellPtr p_apoptotic_cell(new Cell(p_state, p_cell_model));
        p_apoptotic_cell->SetCellProliferativeType(p_stem_type);

        // Set up oxygen_concentration
        p_apoptotic_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);

        // Force the cell to be apoptotic
        for (unsigned i=0; i<num_steps; i++)
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

    void TestContactInhibitionPhaseBasedCellCycleModel() throw(Exception)
    {
        // Check that mQuiescentVolumeFraction and mEquilibriumVolume are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        ContactInhibitionPhaseBasedCellCycleModel* p_model = new ContactInhibitionPhaseBasedCellCycleModel;
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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0*24.0, num_steps);

        // Create cell-cycle models and cells
        ContactInhibitionPhaseBasedCellCycleModel* p_hepa_one_model = new ContactInhibitionPhaseBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);
        p_hepa_one_model->SetBirthTime(0.0);
        p_hepa_one_model->SetQuiescentVolumeFraction(0.5);
        p_hepa_one_model->SetEquilibriumVolume(1.0);

        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);

        // Set up cell volume
        p_hepa_one_cell->GetCellData()->SetItem("volume", hi_volume);
        p_hepa_one_cell->InitialiseCellCycleModel();

        ContactInhibitionPhaseBasedCellCycleModel* p_diff_model = new ContactInhibitionPhaseBasedCellCycleModel;
        p_diff_model->SetDimension(2);
        p_diff_model->SetQuiescentVolumeFraction(0.5);
        p_diff_model->SetEquilibriumVolume(1.0);

        CellPtr p_diff_cell(new Cell(p_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->GetCellData()->SetItem("volume", hi_volume);
        p_diff_cell->InitialiseCellCycleModel();

        // Check that the cell cycle phase and ready to divide are updated correctly
        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(),M_PHASE);

        TS_ASSERT_EQUALS(p_diff_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_diff_model->GetCurrentCellCyclePhase(),G_ZERO_PHASE);

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, p_hepa_one_model->GetG1Duration());
        }

        TS_ASSERT_DELTA(p_hepa_one_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

        // Check that cell division correctly resets the cell cycle phase
        TS_ASSERT_EQUALS(p_hepa_one_cell->ReadyToDivide(), true);

        CellPtr p_hepa_one_cell2 = p_hepa_one_cell->Divide();
        ContactInhibitionPhaseBasedCellCycleModel* p_hepa_one_model2 = static_cast <ContactInhibitionPhaseBasedCellCycleModel*>(p_hepa_one_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_hepa_one_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model->GetCurrentCellCyclePhase(), M_PHASE);
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);
    }

    void TestStochasticOxygenBasedPhaseBasedCellCycleModel() throw(Exception)
    {
        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        StochasticOxygenBasedPhaseBasedCellCycleModel* p_model = new StochasticOxygenBasedPhaseBasedCellCycleModel;
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
        double lo_oxygen_concentration=0.0;
        double hi_oxygen_concentration=1.0;

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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(4.0*18.0, num_steps);

        TS_ASSERT_THROWS_NOTHING(StochasticOxygenBasedPhaseBasedCellCycleModel model);

        // Create cell-cycle model
        StochasticOxygenBasedPhaseBasedCellCycleModel* p_hepa_one_model = new StochasticOxygenBasedPhaseBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);

        StochasticOxygenBasedPhaseBasedCellCycleModel* p_diff_model = new StochasticOxygenBasedPhaseBasedCellCycleModel;
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

        for (unsigned i=0; i<num_steps; i++)
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
        StochasticOxygenBasedPhaseBasedCellCycleModel* p_hepa_one_model2 = static_cast <StochasticOxygenBasedPhaseBasedCellCycleModel*> (p_hepa_one_model->CreateCellCycleModel());
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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_hepa_one_model2->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell-cycle model
        StochasticOxygenBasedPhaseBasedCellCycleModel* p_cell_model = new StochasticOxygenBasedPhaseBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        CellPtr p_apoptotic_cell(new Cell(p_state, p_cell_model));
        p_apoptotic_cell->SetCellProliferativeType(p_stem_type);
        p_apoptotic_cell->GetCellData()->SetItem("oxygen", lo_oxygen_concentration);
        p_apoptotic_cell->InitialiseCellCycleModel();

        // Force the cell to be apoptotic
        for (unsigned i=0; i<num_steps; i++)
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

        StochasticOxygenBasedPhaseBasedCellCycleModel* p_cell_model2 = new StochasticOxygenBasedPhaseBasedCellCycleModel;
        p_cell_model2->SetDimension(2);

        // Coverage
        p_cell_model2->SetMinimumGapDuration(1e20);

        CellPtr p_cell2(new Cell(p_state, p_cell_model2));
        p_cell2->SetCellProliferativeType(p_stem_type);
        p_cell2->InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model2->GetG2Duration(), 1e20, 1e-4);
    }

    void TestArchiveFixedDurationGenerationBasedPhaseBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "FixedDurationGenerationBasedPhaseBasedCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new FixedDurationGenerationBasedPhaseBasedCellCycleModel;

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
            TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveStochasticDurationGenerationBasedPhaseBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StochasticDurationGenerationBasedPhaseBasedCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new StochasticDurationGenerationBasedPhaseBasedCellCycleModel;

            p_model->SetDimension(2);
            static_cast<AbstractPhaseBasedCellCycleModel*>(p_model)->SetTransitCellG1Duration(1.0);

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

    void TestArchiveStochasticDurationPhaseBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StochasticDurationPhaseBasedCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new StochasticDurationPhaseBasedCellCycleModel;
            p_model->SetDimension(2);

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

    void TestArchiveGammaDistributedStochasticDurationPhaseBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "GammaDistributedStochasticDurationPhaseBasedCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new GammaDistributedStochasticDurationPhaseBasedCellCycleModel;
            static_cast<GammaDistributedStochasticDurationPhaseBasedCellCycleModel*>(p_model)->SetShape(2.45);
            static_cast<GammaDistributedStochasticDurationPhaseBasedCellCycleModel*>(p_model)->SetScale(13.42);

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
            TS_ASSERT_DELTA(static_cast<GammaDistributedStochasticDurationPhaseBasedCellCycleModel*>(p_model2)->GetShape(), 2.45, 1e-12);
            TS_ASSERT_DELTA(static_cast<GammaDistributedStochasticDurationPhaseBasedCellCycleModel*>(p_model2)->GetScale(), 13.42, 1e-12);

            // Avoid memory leaks
            delete p_model2;
       }
    }

    void TestArchiveExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() +
                "ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel;
            static_cast<ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel*>(p_model)->SetRate(13.42);

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
            TS_ASSERT_DELTA(static_cast<ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel*>(p_model2)->GetRate(), 13.42, 1e-12);

            // Avoid memory leaks
            delete p_model2;
       }
    }

    void TestArchiveSimpleOxygenBasedPhaseBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SimpleOxygenBasedPhaseBasedCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new SimpleOxygenBasedPhaseBasedCellCycleModel;
            p_model->SetDimension(3);

            static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetHypoxicConcentration(0.8);
            static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetQuiescentConcentration(0.7);
            static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetCriticalHypoxicDuration(2.5);
            static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetCurrentHypoxiaOnsetTime(3.1);

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

            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetHypoxicConcentration(), 0.8, 1e-6);
            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetQuiescentConcentration(), 0.7, 1e-6);
            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetCriticalHypoxicDuration(), 2.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetCurrentHypoxiaOnsetTime(), 3.1, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveStochasticOxygenBasedPhaseBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StochasticOxygenBasedPhaseBasedCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new StochasticOxygenBasedPhaseBasedCellCycleModel;
            p_model->SetDimension(3);

            static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetHypoxicConcentration(0.8);
            static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetQuiescentConcentration(0.7);
            static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetCriticalHypoxicDuration(2.5);
            static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model)->SetCurrentHypoxiaOnsetTime(3.1);
            static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model)->GenerateStochasticG2Duration();

            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model)->GetG2Duration(), 3.0822, 1e-4);

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

            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetG2Duration(), 3.0822, 1e-4);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetHypoxicConcentration(), 0.8, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetQuiescentConcentration(), 0.7, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetCriticalHypoxicDuration(), 2.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedPhaseBasedCellCycleModel*>(p_model2)->GetCurrentHypoxiaOnsetTime(), 3.1, 1e-6);

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveContactInhibitionPhaseBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ContactInhibitionPhaseBasedCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new ContactInhibitionPhaseBasedCellCycleModel();

            p_model->SetDimension(1);
            p_model->SetBirthTime(-1.5);

            static_cast<ContactInhibitionPhaseBasedCellCycleModel*>(p_model)->SetQuiescentVolumeFraction(0.5);
            static_cast<ContactInhibitionPhaseBasedCellCycleModel*>(p_model)->SetEquilibriumVolume(1.0);

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

            TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 1u);
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.5, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.5, 1e-12);

            TS_ASSERT_DELTA(static_cast<ContactInhibitionPhaseBasedCellCycleModel*>(p_model2)->GetQuiescentVolumeFraction(), 0.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<ContactInhibitionPhaseBasedCellCycleModel*>(p_model2)->GetEquilibriumVolume(), 1.0, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with FixedDurationGenerationBasedPhaseBasedCellCycleModel
        FixedDurationGenerationBasedPhaseBasedCellCycleModel fixed_duration_generation_based_cell_cycle_model;
        TS_ASSERT_EQUALS(fixed_duration_generation_based_cell_cycle_model.GetIdentifier(), "FixedDurationGenerationBasedPhaseBasedCellCycleModel");

        out_stream fixed_duration_generation_based_parameter_file = output_file_handler.OpenOutputFile("fixed_duration_generation_based_results.parameters");
        fixed_duration_generation_based_cell_cycle_model.OutputCellCycleModelParameters(fixed_duration_generation_based_parameter_file);
        fixed_duration_generation_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("fixed_duration_generation_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/fixed_duration_generation_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with StochasticDurationGenerationBasedPhaseBasedCellCycleModel
        StochasticDurationGenerationBasedPhaseBasedCellCycleModel stochastic_duration_generation_based_cell_cycle_model;
        TS_ASSERT_EQUALS(stochastic_duration_generation_based_cell_cycle_model.GetIdentifier(), "StochasticDurationGenerationBasedPhaseBasedCellCycleModel");

        out_stream stochastic_duration_generation_based_parameter_file = output_file_handler.OpenOutputFile("stochastic_duration_generation_based_results.parameters");
        stochastic_duration_generation_based_cell_cycle_model.OutputCellCycleModelParameters(stochastic_duration_generation_based_parameter_file);
        stochastic_duration_generation_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("stochastic_duration_generation_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/stochastic_duration_generation_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with StochasticDurationPhaseBasedCellCycleModel
        StochasticDurationPhaseBasedCellCycleModel stochastic_duration_cell_cycle_model;
        TS_ASSERT_EQUALS(stochastic_duration_cell_cycle_model.GetIdentifier(), "StochasticDurationPhaseBasedCellCycleModel");

        out_stream stochastic_duration_parameter_file = output_file_handler.OpenOutputFile("stochastic_duration_results.parameters");
        stochastic_duration_cell_cycle_model.OutputCellCycleModelParameters(stochastic_duration_parameter_file);
        stochastic_duration_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("stochastic_duration_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/stochastic_duration_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with SimpleOxygenBasedPhaseBasedCellCycleModel
        SimpleOxygenBasedPhaseBasedCellCycleModel simple_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(simple_oxygen_based_cell_cycle_model.GetIdentifier(), "SimpleOxygenBasedPhaseBasedCellCycleModel");

        out_stream simple_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("simple_oxygen_based_results.parameters");
        simple_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(simple_oxygen_based_parameter_file);
        simple_oxygen_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("simple_oxygen_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/simple_oxygen_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with StochasticOxygenBasedPhaseBasedCellCycleModel
        StochasticOxygenBasedPhaseBasedCellCycleModel stochastic_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(stochastic_oxygen_based_cell_cycle_model.GetIdentifier(), "StochasticOxygenBasedPhaseBasedCellCycleModel");

        out_stream stochastic_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("stochastic_oxygen_based_results.parameters");
        stochastic_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(stochastic_oxygen_based_parameter_file);
        stochastic_oxygen_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("stochastic_oxygen_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/stochastic_oxygen_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with ContactInhibitionPhaseBasedCellCycleModel
        ContactInhibitionPhaseBasedCellCycleModel contact_inhibition_cell_cycle_model;
        contact_inhibition_cell_cycle_model.SetQuiescentVolumeFraction(0.5);
        contact_inhibition_cell_cycle_model.SetEquilibriumVolume(1.0);
        TS_ASSERT_EQUALS(contact_inhibition_cell_cycle_model.GetIdentifier(), "ContactInhibitionPhaseBasedCellCycleModel");

        out_stream contact_inhibition_parameter_file = output_file_handler.OpenOutputFile("contact_inhibition_results.parameters");
        contact_inhibition_cell_cycle_model.OutputCellCycleModelParameters(contact_inhibition_parameter_file);
        contact_inhibition_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("contact_inhibition_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/contact_inhibition_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with GammaDistributedStochasticDurationPhaseBasedCellCycleModel
        GammaDistributedStochasticDurationPhaseBasedCellCycleModel gamma_cell_cycle_model;
        gamma_cell_cycle_model.SetShape(0.85);
        gamma_cell_cycle_model.SetScale(1.23);
        TS_ASSERT_EQUALS(gamma_cell_cycle_model.GetIdentifier(), "GammaDistributedStochasticDurationPhaseBasedCellCycleModel");

        out_stream gamma_parameter_file = output_file_handler.OpenOutputFile("gamma_results.parameters");
        gamma_cell_cycle_model.OutputCellCycleModelParameters(gamma_parameter_file);
        gamma_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("gamma_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/gamma_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel
        ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel exponential_cell_cycle_model;
        exponential_cell_cycle_model.SetRate(1.23);
        TS_ASSERT_EQUALS(exponential_cell_cycle_model.GetIdentifier(), "ExponentiallyDistributedStochasticDurationGenerationBasedPhaseBasedCellCycleModel");

        out_stream exponential_parameter_file = output_file_handler.OpenOutputFile("exponential_results.parameters");
        exponential_cell_cycle_model.OutputCellCycleModelParameters(exponential_parameter_file);
        exponential_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("exponential_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/exponential_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTSIMPLECELLCYCLEMODELS_HPP_*/
