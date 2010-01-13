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

#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StochasticDivisionRuleCellCycleModel.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"


/**
 * This class contains tests for methods on cell
 * cycle models that are not yet ready for release.
 */
class TestCellCycleModelsNotForRelease : public AbstractCellBasedTestSuite
{
public:

    void TestAlarcon2004OxygenBasedCellCycleModel() throw(Exception)
    {
        TissueConfig::Instance()->SetHepaOneParameters();

        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20.0, 2);

        // Set up oxygen_concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create cell cycle model and associated cell
        Alarcon2004OxygenBasedCellCycleModel* p_cell_model = new Alarcon2004OxygenBasedCellCycleModel(2);
        TissueCell cell(STEM, HEALTHY, p_cell_model);

        // Coverage of cell cycle model copying without an ODE system set up
        TissueCell stem_cell2 = cell;
        TS_ASSERT_EQUALS(stem_cell2.GetMutationState(), HEALTHY);

        cell.InitialiseCellCycleModel();

        // Check oxygen concentration is correct in cell cycle model
        TS_ASSERT_DELTA(p_cell_model->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false);

        // Divide a cell
        Alarcon2004OxygenBasedCellCycleModel* p_cell_model2 = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_cell_model->CreateCellCycleModel());

        TissueCell cell2(STEM, HEALTHY, p_cell_model2);

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), false)
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(), true)
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(), true);

        TS_ASSERT_THROWS_NOTHING(p_cell_model->ResetForDivision());

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        Alarcon2004OxygenBasedCellCycleModel* p_cell_model3 = new Alarcon2004OxygenBasedCellCycleModel(1);
        TissueCell cell3(STEM, HEALTHY, p_cell_model3);
        cell3.InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model3->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();
    }


    void TestSimpleOxygenBasedCellCycleModel() throw(Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();
        p_params->SetHepaOneParameters();

        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration are
        // updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel(2);
        TissueCell cell(STEM, HEALTHY, p_model);
        cell.InitialiseCellCycleModel();

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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0*(p_params->GetHepaOneCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);

        // Set up constant oxygen_concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        TS_ASSERT_THROWS_NOTHING(SimpleOxygenBasedCellCycleModel model(2));

        // Create cell cycle model
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model = new SimpleOxygenBasedCellCycleModel(2);
        SimpleOxygenBasedCellCycleModel* p_diff_model = new SimpleOxygenBasedCellCycleModel(2);

        // Create cell
        TissueCell hepa_one_cell(STEM, HEALTHY, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();

        TissueCell diff_cell(DIFFERENTIATED, HEALTHY, p_diff_model);
        diff_cell.InitialiseCellCycleModel();

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
        TS_ASSERT_EQUALS(hepa_one_cell.ReadyToDivide(), true);
        TissueCell hepa_one_cell2 = hepa_one_cell.Divide();
        SimpleOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast <SimpleOxygenBasedCellCycleModel*>(hepa_one_cell2.GetCellCycleModel());
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*TissueConfig::Instance()->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell cycle model
        SimpleOxygenBasedCellCycleModel* p_cell_model = new SimpleOxygenBasedCellCycleModel(2);
        TissueCell apoptotic_cell(STEM, HEALTHY, p_cell_model);

        // Set up constant oxygen_concentration
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        // Force the cell to be apoptotic
        for (unsigned i=0; i<num_steps; i++)
        {
            TS_ASSERT(apoptotic_cell.GetCellProliferativeType()!=APOPTOTIC ||
                      p_simulation_time->GetTime() >= TissueConfig::Instance()->GetCriticalHypoxicDuration());
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            apoptotic_cell.ReadyToDivide();
        }

        // Test that the cell type is updated to be APOPTOTIC
        TS_ASSERT(apoptotic_cell.GetCellProliferativeType()==APOPTOTIC);
        TS_ASSERT_EQUALS(p_cell_model->GetCurrentHypoxicDuration(), 2.04);

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        SimpleOxygenBasedCellCycleModel* p_cell_model1d = new SimpleOxygenBasedCellCycleModel(1);
        TissueCell cell1d(STEM, HEALTHY, p_cell_model1d);
        cell1d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();

        // For coverage, create a 3D model
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        SimpleOxygenBasedCellCycleModel* p_cell_model3d = new SimpleOxygenBasedCellCycleModel(3);
        TissueCell cell3d(STEM, HEALTHY, p_cell_model3d);
        cell3d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<3>::Destroy();
    }


    void TestStochasticDivisionRuleCellCycleModel() throw(Exception)
    {
        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double end_time = 72.0;
        unsigned num_timesteps = 1000*(unsigned)end_time;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_timesteps);

        // Create a cell
        StochasticDivisionRuleCellCycleModel* p_cycle_model1 = new StochasticDivisionRuleCellCycleModel;
        TissueCell cell1(STEM, HEALTHY, p_cycle_model1);
        cell1.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cycle_model1->GetGeneration(), 0u);

        /**
         * Test with asymmetric division
         */

        TissueConfig::Instance()->SetSymmetricDivisionProbability(0.0);

        // Increment time
        for (unsigned i=0; i<num_timesteps/3; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_cycle_model1, 13.0676);
        }

        // This cell must have divided asymmetrically
        TS_ASSERT_EQUALS(p_cycle_model1->DividedSymmetrically(), false);
        TS_ASSERT_EQUALS(cell1.GetCellProliferativeType(), STEM);
        TS_ASSERT_EQUALS(p_cycle_model1->GetGeneration(), 0u);

        TS_ASSERT_EQUALS(cell1.ReadyToDivide(), true);
        TissueCell cell2 = cell1.Divide();

        TS_ASSERT_EQUALS(cell2.GetCellProliferativeType(), TRANSIT);

        StochasticDivisionRuleCellCycleModel* p_cycle_model2 = static_cast <StochasticDivisionRuleCellCycleModel*> (cell2.GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model2->GetGeneration(), 1u);

        /**
         * Test with symmetric division
         */

        TissueConfig::Instance()->SetSymmetricDivisionProbability(1.0);
        TissueConfig::Instance()->SetMaxTransitGenerations(1);

        StochasticDivisionRuleCellCycleModel* p_cycle_model3 = new StochasticDivisionRuleCellCycleModel;
        TissueCell cell3(STEM, HEALTHY, p_cycle_model3);
        cell3.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cycle_model3->GetGeneration(), 0u);

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
        TS_ASSERT_EQUALS(cell1.ReadyToDivide(), true);
        TissueCell cell4 = cell1.Divide();

        TS_ASSERT_EQUALS(p_cycle_model1->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(cell1.GetCellProliferativeType(), STEM);
        TS_ASSERT_EQUALS(p_cycle_model1->GetGeneration(), 0u);

        StochasticDivisionRuleCellCycleModel* p_cycle_model4 = static_cast <StochasticDivisionRuleCellCycleModel*> (cell4.GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model4->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(cell4.GetCellProliferativeType(), STEM);
        TS_ASSERT_EQUALS(p_cycle_model4->GetGeneration(), 0u);

        // The stem cell cell3 must have divided symmetrically. For coverage,
        // we iterate the random number generator so that it divides into two
        // transit cells
        RandomNumberGenerator::Instance()->ranf();

        TS_ASSERT_EQUALS(cell3.ReadyToDivide(), true);
        TissueCell cell5 = cell3.Divide();

        TS_ASSERT_EQUALS(p_cycle_model3->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(cell3.GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_EQUALS(p_cycle_model3->GetGeneration(), 1u);

        StochasticDivisionRuleCellCycleModel* p_cycle_model5 = static_cast <StochasticDivisionRuleCellCycleModel*> (cell5.GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model5->DividedSymmetrically(), true);
        TS_ASSERT_EQUALS(cell5.GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_EQUALS(p_cycle_model5->GetGeneration(), 1u);

        // The transit cell cell2 divides into two differentiated cells
        TS_ASSERT_EQUALS(cell2.ReadyToDivide(), true);
        TissueCell cell6 = cell2.Divide();

        TS_ASSERT_EQUALS(p_cycle_model2->DividedSymmetrically(), false);
        TS_ASSERT_EQUALS(cell2.GetCellProliferativeType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(p_cycle_model2->GetGeneration(), 2u);

        StochasticDivisionRuleCellCycleModel* p_cycle_model6 = static_cast <StochasticDivisionRuleCellCycleModel*> (cell6.GetCellCycleModel());
        TS_ASSERT_EQUALS(p_cycle_model6->DividedSymmetrically(), false);
        TS_ASSERT_EQUALS(cell6.GetCellProliferativeType(), DIFFERENTIATED);
        TS_ASSERT_EQUALS(p_cycle_model6->GetGeneration(), 2u);

        // For coverage
        StochasticDivisionRuleCellCycleModel* p_cycle_model7 = new StochasticDivisionRuleCellCycleModel;
        TissueCell cell7(DIFFERENTIATED, HEALTHY, p_cycle_model7);
        cell7.InitialiseCellCycleModel();
        TS_ASSERT_EQUALS(p_cycle_model7->GetGeneration(), 0u);
    }


    void TestStochasticOxygenBasedCellCycleModel() throw(Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();
        p_params->SetHepaOneParameters();

        // Check that mCurrentHypoxiaOnsetTime and mCurrentHypoxicDuration
        // are updated correctly
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

        StochasticOxygenBasedCellCycleModel* p_model = new StochasticOxygenBasedCellCycleModel(2);
        TissueCell cell(STEM, HEALTHY, p_model);
        cell.InitialiseCellCycleModel();

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
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0*(p_params->GetHepaOneCellG1Duration()
                  +p_params->GetSG2MDuration()     ), num_steps);

        // Set up constant oxygen_concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        TS_ASSERT_THROWS_NOTHING(StochasticOxygenBasedCellCycleModel model(2));

        // Create cell cycle model
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model = new StochasticOxygenBasedCellCycleModel(2);
        StochasticOxygenBasedCellCycleModel* p_diff_model = new StochasticOxygenBasedCellCycleModel(2);

        // Create cell
        TissueCell hepa_one_cell(STEM, HEALTHY, p_hepa_one_model);
        hepa_one_cell.InitialiseCellCycleModel();

        TissueCell diff_cell(DIFFERENTIATED, HEALTHY, p_diff_model);
        diff_cell.InitialiseCellCycleModel();

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
        TS_ASSERT_EQUALS(hepa_one_cell.ReadyToDivide(), true);
        TissueCell hepa_one_cell_divide = hepa_one_cell.Divide();

        // Check that cell division correctly resets the cell cycle phase
        StochasticOxygenBasedCellCycleModel* p_hepa_one_model2 = static_cast <StochasticOxygenBasedCellCycleModel*> (p_hepa_one_model->CreateCellCycleModel());

        TissueCell hepa_one_cell2(STEM, HEALTHY, p_hepa_one_model2);
        TS_ASSERT_EQUALS(p_hepa_one_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_hepa_one_model2->GetCurrentCellCyclePhase(), M_PHASE);

        // Set up SimulationTime
        SimulationTime::Destroy();
        p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0*TissueConfig::Instance()->GetCriticalHypoxicDuration(), num_steps);

        // Create a cell with a simple oxygen-based cell cycle model
        StochasticOxygenBasedCellCycleModel* p_cell_model = new StochasticOxygenBasedCellCycleModel(2);
        TissueCell apoptotic_cell(STEM, HEALTHY, p_cell_model);
        apoptotic_cell.InitialiseCellCycleModel();

        // Set up constant oxygen_concentration
        CellwiseData<2>::Instance()->SetConstantDataForTesting(low_oxygen_concentration);

        // Force the cell to be apoptotic
        for (unsigned i=0; i<num_steps; i++)
        {
            TS_ASSERT(apoptotic_cell.GetCellProliferativeType()!=APOPTOTIC ||
                      p_simulation_time->GetTime() >= TissueConfig::Instance()->GetCriticalHypoxicDuration());
            p_simulation_time->IncrementTimeOneStep();

            // Note that we need to pass in the updated G1 duration
            apoptotic_cell.ReadyToDivide();
        }

        // Test that the cell type is updated to be APOPTOTIC
        TS_ASSERT_EQUALS(apoptotic_cell.GetCellProliferativeType(), APOPTOTIC);
        TS_ASSERT_EQUALS(p_cell_model->GetCurrentHypoxicDuration(), 2.04);

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        StochasticOxygenBasedCellCycleModel* p_cell_model1d = new StochasticOxygenBasedCellCycleModel(1);
        TissueCell cell1d(STEM, HEALTHY, p_cell_model1d);
        cell1d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();

        // For coverage, create a 3D model
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        StochasticOxygenBasedCellCycleModel* p_cell_model3d = new StochasticOxygenBasedCellCycleModel(3);
        TissueCell cell3d(STEM, HEALTHY, p_cell_model3d);
        cell3d.InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<3>::Destroy();
    }


    void TestArchiveAlarcon2004OxygenBasedCellCycleModels()
    {
        // Set up
        TissueConfig::Instance()->SetHepaOneParameters();

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "alarcon_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 2);

            // Create cell cycle model and associated cell
            Alarcon2004OxygenBasedCellCycleModel* p_cell_model = new Alarcon2004OxygenBasedCellCycleModel(3);
            TissueCell cell(STEM, HEALTHY, p_cell_model);

            cell.InitialiseCellCycleModel();
            cell.GetCellCycleModel()->SetBirthTime(-10.0);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(), false);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(cell.GetCellCycleModel()->ReadyToDivide(), true);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive cell
            TissueCell* const p_cell = &cell;
            output_arch << p_cell;

            // Tidy up
            SimulationTime::Destroy();
        }

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            TissueConfig* inst1 = TissueConfig::Instance();

            inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check that archiving worked correctly
            Alarcon2004OxygenBasedCellCycleModel* p_model = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_cell->GetCellCycleModel());

            TS_ASSERT_EQUALS(p_cell, p_model->GetCell());
            TS_ASSERT_EQUALS(p_model->GetDimension(), 3u);
            TS_ASSERT_EQUALS(p_model->ReadyToDivide(), true);

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -10.0, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 20.0, 1e-12);
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(), 10.0, 1e-12);

            // Tidy up
            delete p_cell;
            CellwiseData<3>::Destroy();
        }
    }


    void TestArchiveSimpleOxygenBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "oxygen_based_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            SimpleOxygenBasedCellCycleModel model(1);

            p_simulation_time->IncrementTimeOneStep();

            model.SetBirthTime(-1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << static_cast<const SimpleOxygenBasedCellCycleModel&>(model);

            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            SimpleOxygenBasedCellCycleModel model(2);
            model.SetBirthTime(-2.0);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> model;

            // Check that archiving worked correctly
            TS_ASSERT_EQUALS(model.GetCurrentCellCyclePhase(), M_PHASE);
            TS_ASSERT_EQUALS(model.GetDimension(), 1u);

            TS_ASSERT_DELTA(model.GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(model.GetAge(), 1.5, 1e-12);

            // Tidy up
            CellwiseData<1>::Destroy();
        }
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
            TS_ASSERT_EQUALS(model.DividedSymmetrically(), true);

            TS_ASSERT_DELTA(model.GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(model.GetAge(), 1.5, 1e-12);
        }
    }


    void TestArchiveStochasticOxygenBasedCellCycleModel() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "stochastic_oxygen_based_cell_cycle.arch";

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // Create cell cycle model and associated cell
            StochasticOxygenBasedCellCycleModel* p_cell_model = new StochasticOxygenBasedCellCycleModel(3);
            TissueCell cell(STEM, HEALTHY, p_cell_model);

            cell.InitialiseCellCycleModel();
            cell.GetCellCycleModel()->SetBirthTime(-1.0);

            p_simulation_time->IncrementTimeOneStep();

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive cell
            TissueCell* const p_cell = &cell;
            output_arch << p_cell;

            // Tidy up
            SimulationTime::Destroy();
        }

        {
            // Set up simulation time
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            TissueConfig* inst1 = TissueConfig::Instance();

            inst1->SetSDuration(101.0);

            TissueCell* p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            // Check that archiving worked correctly
            StochasticOxygenBasedCellCycleModel* p_model = static_cast<StochasticOxygenBasedCellCycleModel*> (p_cell->GetCellCycleModel());
;
            TS_ASSERT_EQUALS(p_cell, p_model->GetCell());
            TS_ASSERT_EQUALS(p_model->GetDimension(), 3u);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), M_PHASE);

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-4);
            TS_ASSERT_DELTA(p_model->GetAge(), 1.5, 1e-4);
            TS_ASSERT_DELTA(p_model->GetG2Duration(), 3.0676, 1e-4); // first random number generated

            // Tidy up
            delete p_cell;
        }
    }

};

#endif /*TESTCELLCYCLEMODELSNOTFORRELEASE_HPP_*/
