/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTSIMPLECELLCYCLEMODELS_HPP_
#define TESTSIMPLECELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"

class TestSimpleCellCycleModels : public AbstractCellBasedTestSuite
{
public:

    void TestFixedDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(FixedDurationGenerationBasedCellCycleModel model3);

        FixedDurationGenerationBasedCellCycleModel* p_stem_model = new FixedDurationGenerationBasedCellCycleModel;
        p_stem_model->SetCellProliferativeType(STEM);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
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

        TS_ASSERT_EQUALS(p_stem_cell->GetCellCycleModel()->GetCellProliferativeType(),STEM);

        FixedDurationGenerationBasedCellCycleModel* p_transit_model = new FixedDurationGenerationBasedCellCycleModel;
        p_transit_model->SetCellProliferativeType(TRANSIT);

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_transit_cell->GetCellCycleModel()->GetCellProliferativeType(),TRANSIT);
        TS_ASSERT_EQUALS(p_transit_model->GetGeneration(), 0u);

        FixedDurationGenerationBasedCellCycleModel* p_diff_model = new FixedDurationGenerationBasedCellCycleModel;
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_diff_cell->GetCellCycleModel()->GetCellProliferativeType(),DIFFERENTIATED);
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

        FixedDurationGenerationBasedCellCycleModel* p_hepa_one_model = new FixedDurationGenerationBasedCellCycleModel;
        p_hepa_one_model->SetCellProliferativeType(STEM);

        // Change G1 Duration for this model
        p_hepa_one_model->SetStemCellG1Duration(8.0);
        p_hepa_one_model->SetTransitCellG1Duration(8.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
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

    void TestStochasticDurationGenerationBasedCellCycleModel() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(StochasticDurationGenerationBasedCellCycleModel cell_model3);

        StochasticDurationGenerationBasedCellCycleModel* p_stem_model = new StochasticDurationGenerationBasedCellCycleModel;
        p_stem_model->SetCellProliferativeType(STEM);

        // Change G1 Duration for this model
        p_stem_model->SetStemCellG1Duration(1.0);

        StochasticDurationGenerationBasedCellCycleModel* p_transit_model = new StochasticDurationGenerationBasedCellCycleModel;
        p_transit_model->SetCellProliferativeType(TRANSIT);

        // Change G1 Duration for this model
        p_transit_model->SetTransitCellG1Duration(1.0);

        StochasticDurationGenerationBasedCellCycleModel* p_diff_model = new StochasticDurationGenerationBasedCellCycleModel;
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->InitialiseCellCycleModel();

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
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
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 4.36075);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 1.78877);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        StochasticDurationGenerationBasedCellCycleModel* p_hepa_one_model = new StochasticDurationGenerationBasedCellCycleModel;
        p_hepa_one_model->SetCellProliferativeType(STEM);
        // Change G1 Duration for this model
        p_hepa_one_model->SetStemCellG1Duration(1.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i< num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 4.1324);
        }
    }


    void TestStochasticDurationCellCycleModel() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(StochasticDurationCellCycleModel cell_model3);

        StochasticDurationCellCycleModel* p_stem_model = new StochasticDurationCellCycleModel;
        p_stem_model->SetCellProliferativeType(STEM);

        // Change G1 Duration for this model
        p_stem_model->SetStemCellG1Duration(8.0);

        StochasticDurationCellCycleModel* p_transit_model = new StochasticDurationCellCycleModel;
        p_transit_model->SetCellProliferativeType(TRANSIT);

        // Change G1 Duration for this model
        p_stem_model->SetTransitCellG1Duration(8.0);

        StochasticDurationCellCycleModel* p_diff_model = new StochasticDurationCellCycleModel;
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->InitialiseCellCycleModel();

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
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
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 9.68038);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 2.78877);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        StochasticDurationCellCycleModel* p_hepa_one_model = new StochasticDurationCellCycleModel;
        p_hepa_one_model->SetCellProliferativeType(STEM);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->InitialiseCellCycleModel();

        for (unsigned i=0; i< num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 15.5662);
        }
    }


    void TestSimpleOxygenBasedCellCycleModel() throw(Exception)
    {
        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
        p_model->SetDimension(2);
        p_model->SetCellProliferativeType(STEM);

        p_model->SetStemCellG1Duration(8.0);
        p_model->SetTransitCellG1Duration(8.0);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        CellPtr p_cell(new Cell(p_state, p_model));

        p_cell->InitialiseCellCycleModel();

        // Set up constant oxygen_concentration
        std::vector<double> low_oxygen_concentration;
        std::vector<double> high_oxygen_concentration;
        low_oxygen_concentration.push_back(0.0);
        high_oxygen_concentration.push_back(1.0);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(high_oxygen_concentration);

        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);
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

        // Set up constant oxygen concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        TS_ASSERT_THROWS_NOTHING(SimpleOxygenBasedCellCycleModel model);

        // Create cell-cycle models and cells
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model = new SimpleOxygenBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);
        p_hepa_one_model->SetCellProliferativeType(STEM);

        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->InitialiseCellCycleModel();

        SimpleOxygenBasedCellCycleModel* p_diff_model = new SimpleOxygenBasedCellCycleModel;
        p_diff_model->SetDimension(2);
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);

        // Coverage
        TS_ASSERT_DELTA(p_diff_model->GetCriticalHypoxicDuration(), 2.0, 1e-6);
        p_diff_model->SetCriticalHypoxicDuration(0.5);
        TS_ASSERT_DELTA(p_diff_model->GetCriticalHypoxicDuration(), 0.5, 1e-6);

        CellPtr p_diff_cell(new Cell(p_state, p_diff_model));
        p_diff_cell->InitialiseCellCycleModel();

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
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast <SimpleOxygenBasedCellCycleModel*>(p_hepa_one_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_hepa_one_model2->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell-cycle model
        SimpleOxygenBasedCellCycleModel* p_cell_model = new SimpleOxygenBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        p_cell_model->SetCellProliferativeType(STEM);
        CellPtr p_apoptotic_cell(new Cell(p_state, p_cell_model));

        // Set up constant oxygen_concentration
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

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

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        SimpleOxygenBasedCellCycleModel* p_cell_model1d = new SimpleOxygenBasedCellCycleModel;
        p_cell_model1d->SetDimension(1);
        p_cell_model1d->SetCellProliferativeType(STEM);
        CellPtr p_cell1d(new Cell(p_state, p_cell_model1d));

        p_cell1d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();

        // For coverage, create a 3D model
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        SimpleOxygenBasedCellCycleModel* p_cell_model3d = new SimpleOxygenBasedCellCycleModel;
        p_cell_model3d->SetDimension(3);
        p_cell_model3d->SetCellProliferativeType(STEM);
        CellPtr p_cell3d(new Cell(p_state, p_cell_model3d));

        p_cell3d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<3>::Destroy();
    }

    void TestContactInhibitionCellCycleModel() throw(Exception)
    {
        // Check that mQuiescentVolumeFraction and mEquilibriumVolume are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
        p_model->SetDimension(2);
        p_model->SetCellProliferativeType(STEM);

        p_model->SetStemCellG1Duration(8.0);
        p_model->SetTransitCellG1Duration(8.0);

        TS_ASSERT_THROWS_THIS(p_model->UpdateCellCyclePhase(),
            "The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");

        p_model->SetQuiescentVolumeFraction(0.5);
        p_model->SetEquilibriumVolume(1.0);

        // Set the birth time such that at t=0, the cell has just entered G1 phase
        p_model->SetBirthTime(-1.0);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        CellPtr p_cell(new Cell(p_state, p_model));

        p_cell->InitialiseCellCycleModel();

        // Set up constant volume for inhibition test
        std::vector<double> low_volume;
        std::vector<double> high_volume;
        low_volume.push_back(0.0);
        high_volume.push_back(1.0);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_volume);

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentOnsetTime(), 0.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(high_volume);

        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentQuiescentOnsetTime(), 2.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_volume);
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

        // Set up constant cell volume
        std::vector<double> cell_volume;
        cell_volume.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(cell_volume);

        // Create cell-cycle models and cells
        ContactInhibitionCellCycleModel* p_hepa_one_model = new ContactInhibitionCellCycleModel;
        p_hepa_one_model->SetDimension(2);
        p_hepa_one_model->SetCellProliferativeType(STEM);
        p_hepa_one_model->SetBirthTime(0.0);
        p_hepa_one_model->SetQuiescentVolumeFraction(0.5);
        p_hepa_one_model->SetEquilibriumVolume(1.0);

        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->InitialiseCellCycleModel();

        ContactInhibitionCellCycleModel* p_diff_model = new ContactInhibitionCellCycleModel;
        p_diff_model->SetDimension(2);
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);
        p_diff_model->SetQuiescentVolumeFraction(0.5);
        p_diff_model->SetEquilibriumVolume(1.0);

        CellPtr p_diff_cell(new Cell(p_state, p_diff_model));
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
        ContactInhibitionCellCycleModel* p_hepa_one_model2 = static_cast <ContactInhibitionCellCycleModel*>(p_hepa_one_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(cell_volume);

        ContactInhibitionCellCycleModel* p_cell_model1d = new ContactInhibitionCellCycleModel;
        p_cell_model1d->SetDimension(1);
        p_cell_model1d->SetCellProliferativeType(STEM);
        p_cell_model1d->SetQuiescentVolumeFraction(0.5);
        p_cell_model1d->SetEquilibriumVolume(1.0);

        CellPtr p_cell1d(new Cell(p_state, p_cell_model1d));
        p_cell1d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();

        // For coverage, create a 3D model
        CellwiseData<3>::Instance()->SetConstantDataForTesting(cell_volume);

        ContactInhibitionCellCycleModel* p_cell_model3d = new ContactInhibitionCellCycleModel;
        p_cell_model3d->SetDimension(3);
        p_cell_model3d->SetCellProliferativeType(STEM);
        p_cell_model3d->SetQuiescentVolumeFraction(0.5);
        p_cell_model3d->SetEquilibriumVolume(1.0);

        CellPtr p_cell3d(new Cell(p_state, p_cell_model3d));
        p_cell3d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<3>::Destroy();
    }

    void TestStochasticOxygenBasedCellCycleModel() throw(Exception)
    {
        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel;
        p_model->SetDimension(2);
        p_model->SetCellProliferativeType(STEM);

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
        CellPtr p_cell(new Cell(p_state, p_model));

        p_cell->InitialiseCellCycleModel();

        // Set up constant oxygen_concentration
        std::vector<double> low_oxygen_concentration;
        std::vector<double> high_oxygen_concentration;
        low_oxygen_concentration.push_back(0.0);
        high_oxygen_concentration.push_back(1.0);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=1.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 0.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(high_oxygen_concentration);

        p_simulation_time->IncrementTimeOneStep(); // t=2.0
        p_model->ReadyToDivide();
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxicDuration(), 0.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetCurrentHypoxiaOnsetTime(), 2.0, 1e-12);

        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);
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

        // Set up constant oxygen_concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        TS_ASSERT_THROWS_NOTHING(StochasticOxygenBasedCellCycleModel model);

        // Create cell-cycle model
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model = new StochasticOxygenBasedCellCycleModel;
        p_hepa_one_model->SetDimension(2);
        p_hepa_one_model->SetCellProliferativeType(STEM);

        StochasticOxygenBasedCellCycleModel* p_diff_model = new StochasticOxygenBasedCellCycleModel;
        p_diff_model->SetDimension(2);
        p_diff_model->SetCellProliferativeType(DIFFERENTIATED);

        // Create cell
        CellPtr p_hepa_one_cell(new Cell(p_state, p_hepa_one_model));
        p_hepa_one_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_state, p_diff_model));
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
        CellPtr p_hepa_one_cell_divide = p_hepa_one_cell->Divide();

        // Check that cell division correctly resets the cell cycle phase
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast <StochasticOxygenBasedCellCycleModel*> (p_hepa_one_model->CreateCellCycleModel());
        p_hepa_one_model2->SetCellProliferativeType(STEM);
        CellPtr p_hepa_one_cell2(new Cell(p_state, p_hepa_one_model2));
        p_hepa_one_cell2->InitialiseCellCycleModel();
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*p_hepa_one_model2->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell-cycle model
        StochasticOxygenBasedCellCycleModel* p_cell_model = new StochasticOxygenBasedCellCycleModel;
        p_cell_model->SetDimension(2);
        p_cell_model->SetCellProliferativeType(STEM);
        CellPtr p_apoptotic_cell(new Cell(p_state, p_cell_model));

        p_apoptotic_cell->InitialiseCellCycleModel();

        // Set up constant oxygen_concentration
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

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

        StochasticOxygenBasedCellCycleModel* p_cell_model2 = new StochasticOxygenBasedCellCycleModel;
        p_cell_model2->SetDimension(2);
        p_cell_model2->SetCellProliferativeType(STEM);

        // Coverage
        p_cell_model2->SetMinimumGapDuration(1e20);

        CellPtr p_cell2(new Cell(p_state, p_cell_model2));
        p_cell2->InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model2->GetG2Duration(), 1e20, 1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        StochasticOxygenBasedCellCycleModel* p_cell_model1d = new StochasticOxygenBasedCellCycleModel;
        p_cell_model1d->SetDimension(1);
        p_cell_model1d->SetCellProliferativeType(STEM);
        CellPtr p_cell1d(new Cell(p_state, p_cell_model1d));

        p_cell1d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();

        // For coverage, create a 3D model
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        StochasticOxygenBasedCellCycleModel* p_cell_model3d = new StochasticOxygenBasedCellCycleModel;
        p_cell_model3d->SetDimension(3);
        p_cell_model3d->SetCellProliferativeType(STEM);
        CellPtr p_cell3d(new Cell(p_state, p_cell_model3d));

        p_cell3d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<3>::Destroy();
    }

    void TestArchiveFixedDurationGenerationBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "FixedDurationGenerationBasedCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);


            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new FixedDurationGenerationBasedCellCycleModel;

            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(TRANSIT);
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

            // Check private data has been restored correctly.
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.0, 1e-12);
            TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(p_model2->GetCellProliferativeType(), TRANSIT);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveStochasticDurationGenerationBasedCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StochasticDurationGenerationBasedCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new StochasticDurationGenerationBasedCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetTransitCellG1Duration(1.0);

            static_cast<StochasticDurationGenerationBasedCellCycleModel*>(p_model)->SetG1Duration();
            TS_ASSERT_DELTA(p_model->GetG1Duration(), 2.6803, 1e-4);

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

            TS_ASSERT_DELTA(p_model2->GetG1Duration(), 2.6803, 1e-4);
            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveStochasticDurationCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StochasticDurationCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new StochasticDurationCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(TRANSIT);

            static_cast<StochasticDurationCellCycleModel*>(p_model)->SetG1Duration();
            TS_ASSERT_DELTA(p_model->GetG1Duration(), 3.6803, 1e-4);

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

            TS_ASSERT_DELTA(p_model2->GetG1Duration(), 3.6803, 1e-4);
            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
       }
    }

    void TestArchiveSimpleOxygenBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SimpleOxygenBasedCellCycleModel.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(3);
            p_model->SetCellProliferativeType(STEM);

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
            CellwiseData<3>::Destroy();
        }
    }

    void TestArchiveStochasticOxygenBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StochasticOxygenBasedCellCycleModel.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new StochasticOxygenBasedCellCycleModel;
            p_model->SetDimension(3);
            p_model->SetCellProliferativeType(STEM);

            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetHypoxicConcentration(0.8);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetQuiescentConcentration(0.7);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetCriticalHypoxicDuration(2.5);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->SetCurrentHypoxiaOnsetTime(3.1);
            static_cast<StochasticOxygenBasedCellCycleModel*>(p_model)->GenerateStochasticG2Duration();

            TS_ASSERT_DELTA(p_model->GetG2Duration(), 5.1122, 1e-4);

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

            TS_ASSERT_DELTA(p_model2->GetG2Duration(), 5.1122, 1e-4);

            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetHypoxicConcentration(), 0.8, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetQuiescentConcentration(), 0.7, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetCriticalHypoxicDuration(), 2.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<StochasticOxygenBasedCellCycleModel*>(p_model2)->GetCurrentHypoxiaOnsetTime(), 3.1, 1e-6);

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
            CellwiseData<3>::Destroy();
        }
    }

    void TestArchiveContactInhibitionCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ContactInhibitionCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new ContactInhibitionCellCycleModel();

            p_model->SetDimension(1);
            p_model->SetCellProliferativeType(STEM);
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

            TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 1u);
            TS_ASSERT_EQUALS(p_model2->GetCellProliferativeType(), STEM);
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.5, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.5, 1e-12);

            TS_ASSERT_DELTA(static_cast<ContactInhibitionCellCycleModel*>(p_model2)->GetQuiescentVolumeFraction(), 0.5, 1e-6);
            TS_ASSERT_DELTA(static_cast<ContactInhibitionCellCycleModel*>(p_model2)->GetEquilibriumVolume(), 1.0, 1e-6);

            // Avoid memory leaks
            delete p_model2;
            CellwiseData<1>::Destroy();
        }
    }

    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with FixedDurationGenerationBasedCellCycleModel
        FixedDurationGenerationBasedCellCycleModel fixed_duration_generation_based_cell_cycle_model;
        TS_ASSERT_EQUALS(fixed_duration_generation_based_cell_cycle_model.GetIdentifier(), "FixedDurationGenerationBasedCellCycleModel");

        out_stream fixed_duration_generation_based_parameter_file = output_file_handler.OpenOutputFile("fixed_duration_generation_based_results.parameters");
        fixed_duration_generation_based_cell_cycle_model.OutputCellCycleModelParameters(fixed_duration_generation_based_parameter_file);
        fixed_duration_generation_based_parameter_file->close();

        std::string fixed_duration_generation_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + fixed_duration_generation_based_results_dir + "fixed_duration_generation_based_results.parameters cell_based/test/data/TestCellCycleModels/fixed_duration_generation_based_results.parameters").c_str()), 0);

        // Test with StochasticDurationGenerationBasedCellCycleModel
        StochasticDurationGenerationBasedCellCycleModel stochastic_duration_generation_based_cell_cycle_model;
        TS_ASSERT_EQUALS(stochastic_duration_generation_based_cell_cycle_model.GetIdentifier(), "StochasticDurationGenerationBasedCellCycleModel");

        out_stream stochastic_duration_generation_based_parameter_file = output_file_handler.OpenOutputFile("stochastic_duration_generation_based_results.parameters");
        stochastic_duration_generation_based_cell_cycle_model.OutputCellCycleModelParameters(stochastic_duration_generation_based_parameter_file);
        stochastic_duration_generation_based_parameter_file->close();

        std::string stochastic_duration_generation_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + stochastic_duration_generation_based_results_dir + "stochastic_duration_generation_based_results.parameters cell_based/test/data/TestCellCycleModels/stochastic_duration_generation_based_results.parameters").c_str()), 0);

        // Test with StochasticDurationCellCycleModel
        StochasticDurationCellCycleModel stochastic_duration_cell_cycle_model;
        TS_ASSERT_EQUALS(stochastic_duration_cell_cycle_model.GetIdentifier(), "StochasticDurationCellCycleModel");

        out_stream stochastic_duration_parameter_file = output_file_handler.OpenOutputFile("stochastic_duration_results.parameters");
        stochastic_duration_cell_cycle_model.OutputCellCycleModelParameters(stochastic_duration_parameter_file);
        stochastic_duration_parameter_file->close();

        std::string stochastic_duration_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + stochastic_duration_results_dir + "stochastic_duration_results.parameters cell_based/test/data/TestCellCycleModels/stochastic_duration_results.parameters").c_str()), 0);

        // Test with SimpleOxygenBasedCellCycleModel
        SimpleOxygenBasedCellCycleModel simple_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(simple_oxygen_based_cell_cycle_model.GetIdentifier(), "SimpleOxygenBasedCellCycleModel");

        out_stream simple_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("simple_oxygen_based_results.parameters");
        simple_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(simple_oxygen_based_parameter_file);
        simple_oxygen_based_parameter_file->close();

        std::string simple_oxygen_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + simple_oxygen_based_results_dir + "simple_oxygen_based_results.parameters cell_based/test/data/TestCellCycleModels/simple_oxygen_based_results.parameters").c_str()), 0);

        // Test with StochasticOxygenBasedCellCycleModel
        StochasticOxygenBasedCellCycleModel stochastic_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(stochastic_oxygen_based_cell_cycle_model.GetIdentifier(), "StochasticOxygenBasedCellCycleModel");

        out_stream stochastic_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("stochastic_oxygen_based_results.parameters");
        stochastic_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(stochastic_oxygen_based_parameter_file);
        stochastic_oxygen_based_parameter_file->close();

        std::string stochastic_oxygen_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + stochastic_oxygen_based_results_dir + "stochastic_oxygen_based_results.parameters cell_based/test/data/TestCellCycleModels/stochastic_oxygen_based_results.parameters").c_str()), 0);

        // Test with ContactInhibitionCellCycleModel
        ContactInhibitionCellCycleModel contact_inhibition_cell_cycle_model;
        contact_inhibition_cell_cycle_model.SetQuiescentVolumeFraction(0.5);
        contact_inhibition_cell_cycle_model.SetEquilibriumVolume(1.0);
        TS_ASSERT_EQUALS(contact_inhibition_cell_cycle_model.GetIdentifier(), "ContactInhibitionCellCycleModel");

        out_stream contact_inhibition_parameter_file = output_file_handler.OpenOutputFile("contact_inhibition_results.parameters");
        contact_inhibition_cell_cycle_model.OutputCellCycleModelParameters(contact_inhibition_parameter_file);
        contact_inhibition_parameter_file->close();

        std::string contact_inhibition_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + contact_inhibition_results_dir + "contact_inhibition_results.parameters cell_based/test/data/TestCellCycleModels/contact_inhibition_results.parameters").c_str()), 0);
    }
};

#endif /*TESTSIMPLECELLCYCLEMODELS_HPP_*/
