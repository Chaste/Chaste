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
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "Goldbeter1991SrnModel.hpp"
#include "NullSrnModel.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellAncestor.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCell: public AbstractCellBasedTestSuite
{
public:

    void TestCellConstructor()
    {
        // Coverage
        TS_ASSERT_THROWS_NOTHING(null_deleter());

        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200, 20);

        // Create a cell mutation state
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);

        // Create a cell-cycle model
        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

        // Create a cell
        CellPtr p_cell(new Cell(p_state, p_model));

        // Test coverage
        TS_ASSERT_THROWS_THIS(p_cell->SetCellProliferativeType(p_state), "Attempting to give cell a cell proliferative type that is not a subtype of AbstractCellProliferativeType");

        p_cell->SetCellProliferativeType(p_type);
        p_cell->SetBirthTime(-0.5);

        // Test members of cell directly
        TS_ASSERT_EQUALS(p_cell->GetCellId(), 0u);
        TS_ASSERT_DELTA(p_cell->GetBirthTime(), -0.5, 1e-6);

        // Test members of cell through cell-cycle model
        TS_ASSERT_EQUALS(p_cell->GetCellCycleModel()->GetCell()->GetCellId(), 0u);
        TS_ASSERT_DELTA(p_cell->GetCellCycleModel()->GetCell()->GetBirthTime(), -0.5, 1e-6);

        // Check using default SRN
        TS_ASSERT(dynamic_cast<NullSrnModel*>(p_cell->GetSrnModel()));
    }

    void TestCellConstructorWithSrn()
    {
        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(100, 1000);

        // Create a cell mutation state
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);

        // Create a cell-cycle model
        FixedG1GenerationalCellCycleModel* p_cc_model = new FixedG1GenerationalCellCycleModel();

        // Create a SRN model
        Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();

        // Create a cell
        CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->SetBirthTime(-0.5);

        // Test members of cell directly
        TS_ASSERT_EQUALS(p_cell->GetCellId(), 0u);
        TS_ASSERT_DELTA(p_cell->GetBirthTime(), -0.5, 1e-6);

        // Test members of cell through SRN model
        TS_ASSERT_EQUALS(p_cell->GetCellId(), p_cell->GetCellCycleModel()->GetCell()->GetCellId());
        TS_ASSERT_DELTA(p_cell->GetBirthTime(), p_cell->GetCellCycleModel()->GetCell()->GetBirthTime(), 1e-6);

        //For coverage change the SRN
        NullSrnModel* p_null_srn_model = new NullSrnModel();
        p_cell->SetSrnModel(p_null_srn_model);

        // Check SRN is updated
        TS_ASSERT(dynamic_cast<NullSrnModel*>(p_cell->GetSrnModel()));
    }

    void TestWithCellPropertyCollection()
    {
        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(24, 3);

        // Create a cell mutation state
        boost::shared_ptr<AbstractCellProperty> p_wild_type(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        // Create a cell-cycle model
        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellLabel, p_label_2);
        MAKE_PTR(ApcTwoHitCellMutationState, p_apc2_mutation);
        MAKE_PTR(StemCellProliferativeType, p_type);

        // Create a cell property collection (for the time being, populate with mutation states)
        CellPropertyCollection collection;
        collection.AddProperty(p_label);

        // Create a cell
        CellPtr p_cell(new Cell(p_wild_type, p_model, NULL, false, collection)); // NUll -> NullSrnModel
        p_cell->SetCellProliferativeType(p_type);
        p_cell->SetBirthTime(-1.0);
        p_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(static_cast<FixedG1GenerationalCellCycleModel*>(p_cell->GetCellCycleModel())->GetGeneration(), 0u);

        // Test cell property collection
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().GetSize(), 5u);
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasProperty(p_wild_type), true);
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasProperty(p_label), true);
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasProperty(p_apc2_mutation), false);
        TS_ASSERT_EQUALS(p_cell->HasCellProperty<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_cell->HasCellProperty<CellLabel>(), true);
        TS_ASSERT_EQUALS(p_cell->HasCellProperty<ApcTwoHitCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasPropertyType<AbstractCellProperty>(), true);
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasPropertyType<AbstractCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasPropertyType<CellId>(), true);

        // Test creating another cell with the same cell property collection
        TS_ASSERT_EQUALS(p_wild_type->GetCellCount(), 1u);

        FixedG1GenerationalCellCycleModel* p_model2 = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell2(new Cell(p_wild_type, p_model2, NULL, false, collection));
        p_cell2->SetCellProliferativeType(p_type);
        p_cell2->SetBirthTime(-1.0);
        p_cell2->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_wild_type->GetCellCount(), 2u);

        // Test adding a cell property
        boost::shared_ptr<AbstractCellProperty> p_apoptotic(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());

        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasProperty(p_apoptotic), false);
        TS_ASSERT_EQUALS(p_apoptotic->GetCellCount(), 0u);

        p_cell->AddCellProperty(p_apoptotic);

        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasProperty(p_apoptotic), true);
        TS_ASSERT_EQUALS(p_apoptotic->GetCellCount(), 1u);

        CellPropertyCollection apoptotic_collection = p_cell->rGetCellPropertyCollection().GetProperties<ApoptoticCellProperty>();
        TS_ASSERT_EQUALS(apoptotic_collection.GetProperty()->GetCellCount(), 1u);

        // Test removing a cell property
        p_cell->RemoveCellProperty<ApoptoticCellProperty>();

        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasProperty(p_apoptotic), false);
        TS_ASSERT_EQUALS(p_apoptotic->GetCellCount(), 0u);

        // Now age cell
        p_simulation_time->IncrementTimeOneStep(); // t=8
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep(); // t=16
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep(); // t=24
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        // Divide cell
        CellPtr p_daughter_cell = p_cell->Divide();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(static_cast<FixedG1GenerationalCellCycleModel*>(p_daughter_cell->GetCellCycleModel())->GetGeneration(), 1u);
        TS_ASSERT_EQUALS(p_daughter_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
        TS_ASSERT_DELTA(p_daughter_cell->GetAge(), 0.0, 1e-9);

        // Test cell property collection has been inherited correctly during division
        TS_ASSERT_EQUALS(p_daughter_cell->rGetCellPropertyCollection().GetSize(), 5u);

        TS_ASSERT_EQUALS(p_daughter_cell->rGetCellPropertyCollection().HasProperty(p_wild_type), true);
        TS_ASSERT_EQUALS(p_daughter_cell->rGetCellPropertyCollection().HasProperty(p_apc2_mutation), false);
        TS_ASSERT_EQUALS(p_daughter_cell->HasCellProperty<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_daughter_cell->HasCellProperty<ApcOneHitCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_daughter_cell->HasCellProperty<ApcTwoHitCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_daughter_cell->rGetCellPropertyCollection().HasPropertyType<AbstractCellProperty>(), true);
        TS_ASSERT_EQUALS(p_daughter_cell->rGetCellPropertyCollection().HasPropertyType<AbstractCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasPropertyType<CellId>(), true);

        TS_ASSERT(&(p_daughter_cell->rGetCellPropertyCollection()) != &(p_cell->rGetCellPropertyCollection()));

        TS_ASSERT_EQUALS(p_cell->HasCellProperty<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_daughter_cell->HasCellProperty<WildTypeCellMutationState>(), true);

        // Check that the cell property registry gets updated correctly when cells are killed
        unsigned num_wild_type = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>()->GetCellCount();
        TS_ASSERT_EQUALS(num_wild_type, 3u);

        p_daughter_cell->Kill();
        num_wild_type = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>()->GetCellCount();
        TS_ASSERT_EQUALS(num_wild_type, 2u);

        p_cell->Kill();
        num_wild_type = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>()->GetCellCount();
        TS_ASSERT_EQUALS(num_wild_type, 1u);

        p_cell2->Kill();
        num_wild_type = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>()->GetCellCount();
        TS_ASSERT_EQUALS(num_wild_type, 0u);
    }

    void TestCellsAgeingCorrectly()
    {
        // These lines are added to cover the exception case that a cell is
        // created without simulation time being set up...
        FixedG1GenerationalCellCycleModel fixed_model;
        SimulationTime::Destroy();

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        TS_ASSERT_THROWS_THIS(CellPtr p_bad_cell(new Cell(p_healthy_state, &fixed_model)),
                              "Cell is setting up a cell-cycle model but SimulationTime has not been set up");

        // Cell wasn't created - count should be zero
        TS_ASSERT_EQUALS(p_healthy_state->GetCellCount(), 0u);

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

        TS_ASSERT_THROWS_THIS(CellPtr p_stem_cell(new Cell(p_healthy_state, NULL)),
                              "Cell-cycle model is null");

        // Cell wasn't created - count should be zero
        TS_ASSERT_EQUALS(p_healthy_state->GetCellCount(), 0u);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

        // Create a SRN model
        NullSrnModel* p_srn_model = new NullSrnModel();

        TS_ASSERT_THROWS_THIS(CellPtr p_another_bad_cell(new Cell(p_label, &fixed_model,p_srn_model)),
                              "Attempting to create cell with a cell mutation state that is not a subtype of AbstractCellMutationState");

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model,p_srn_model));
        p_stem_cell->SetCellProliferativeType(p_type);
        p_stem_cell->InitialiseCellCycleModel();
        p_stem_cell->InitialiseSrnModel();

        TS_ASSERT_THROWS_THIS(p_stem_cell->SetMutationState(p_label),
                              "Attempting to give cell a cell mutation state that is not a subtype of AbstractCellMutationState");

        // Cell was created - count should be one
        TS_ASSERT_EQUALS(p_healthy_state->GetCellCount(), 1u);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_stem_cell->GetAge(), 0.5);

        p_simulation_time->IncrementTimeOneStep();
        p_stem_cell->SetBirthTime(p_simulation_time->GetTime());
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_stem_cell->GetAge(), 1.0);

        TS_ASSERT_EQUALS(p_stem_cell->IsDead(), false);
        p_stem_cell->Kill();
        TS_ASSERT_EQUALS(p_stem_cell->IsDead(), true);

        // Coverage of operator equals
        FixedG1GenerationalCellCycleModel* p_model2 = new FixedG1GenerationalCellCycleModel();
        NullSrnModel* p_srn_model2 = new NullSrnModel();
        CellPtr p_live_cell(new Cell(p_healthy_state, p_model2,p_srn_model2));
        p_live_cell->SetCellProliferativeType(p_type);

        TS_ASSERT_EQUALS(p_live_cell->IsDead(), false);
        p_live_cell = p_stem_cell;
        TS_ASSERT_EQUALS(p_live_cell->IsDead(), true);
    }

    void TestCellDivision()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);

        // We are going to start at t=0 and jump up in steps of 6.0
        p_simulation_time->IncrementTimeOneStep(); //t=6

        // Cover bad cell-cycle model
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_type_transit(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        TS_ASSERT_THROWS_THIS(CellPtr p_bad_cell2(new Cell(p_healthy_state, NULL)), "Cell-cycle model is null");

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model,p_srn_model));
        p_stem_cell->SetCellProliferativeType(p_type);
        p_stem_cell->InitialiseCellCycleModel();
        p_stem_cell->InitialiseSrnModel();

        // Set the time over which the cell would undergo apoptosis, if it were told to
        p_stem_cell->SetApoptosisTime(15.78);
        TS_ASSERT_DELTA(p_stem_cell->GetApoptosisTime(), 15.78, 1e-6);

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_model->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_model->GetSG2MDuration(), 10.0, 1e-12);

        // Test coverage of operator=
        FixedG1GenerationalCellCycleModel* p_model2 = new FixedG1GenerationalCellCycleModel();
        Goldbeter1991SrnModel* p_srn_model2 = new Goldbeter1991SrnModel();
        CellPtr p_other_cell(new Cell(p_healthy_state, p_model2, p_srn_model2));
        p_other_cell->SetCellProliferativeType(p_type_transit);
        p_other_cell->InitialiseCellCycleModel();
        p_other_cell->InitialiseSrnModel();

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_model2->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_model2->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_model2->GetSG2MDuration(), 10.0, 1e-12);

        TS_ASSERT_EQUALS(p_other_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
        p_other_cell = p_stem_cell;
        TS_ASSERT_EQUALS(p_other_cell->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);

        // Back to the test
        p_simulation_time->IncrementTimeOneStep(); //t=12
        p_simulation_time->IncrementTimeOneStep(); //t=18
        p_simulation_time->IncrementTimeOneStep(); //t=24

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep(); //t=30

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);

        // Create transit progeny of stem
        CellPtr p_daughter_cell = p_stem_cell->Divide();

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(static_cast<FixedG1GenerationalCellCycleModel*>(p_daughter_cell->GetCellCycleModel())->GetGeneration(), 1u);
        TS_ASSERT(dynamic_cast<Goldbeter1991SrnModel*>(p_daughter_cell->GetSrnModel()));
        TS_ASSERT_EQUALS(p_daughter_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(), true);
        TS_ASSERT_DELTA(p_daughter_cell->GetAge(), 0, 1e-9);

        // Test that the stem cell's progeny has correctly inherited the apoptosis time member variable
        TS_ASSERT_DELTA(p_stem_cell->GetApoptosisTime(), 15.78, 1e-6);
        TS_ASSERT_DELTA(p_daughter_cell->GetApoptosisTime(), 15.78, 1e-6);

        p_simulation_time->IncrementTimeOneStep(); //t=36

        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep(); //t=42

        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), true);

        // Create transit progeny of transit
        CellPtr p_grandaughter_cell = p_daughter_cell->Divide();

        p_simulation_time->IncrementTimeOneStep(); //t=48

        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_grandaughter_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), false);

        // Stem cell ready to divide again
        p_simulation_time->IncrementTimeOneStep(); //t=54
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);

        // Both grandaughter and daughter cells should be ready to divide
        TS_ASSERT_EQUALS(p_grandaughter_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_daughter_cell->ReadyToDivide(), true);
    }

    void Test0DBucket()
    {
        double end_time = 61.0;
        int time_steps = 61;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, time_steps);
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
        p_stem_cell->SetCellProliferativeType(p_type);
        p_stem_cell->InitialiseCellCycleModel();

        std::vector<CellPtr> cells;
        std::vector<CellPtr> newly_born;
        std::vector<unsigned> stem_cells(time_steps);
        std::vector<unsigned> transit_cells(time_steps);
        std::vector<unsigned> differentiated_cells(time_steps);
        std::vector<double> times(time_steps);

        cells.push_back(p_stem_cell);
        std::vector<CellPtr>::iterator cell_iterator;

        unsigned i=0;
        while (p_simulation_time->GetTime()< end_time)
        {
            // Produce the offspring of the cells
            p_simulation_time->IncrementTimeOneStep();
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->ReadyToDivide())
                {
                    newly_born.push_back((*cell_iterator)->Divide());
                }
                ++cell_iterator;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Update cell counts
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                {
                    stem_cells[i]++;
                }
                else if ((*cell_iterator)->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                {
                    transit_cells[i]++;
                }
                else if ((*cell_iterator)->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
                {
                    differentiated_cells[i]++;
                }

                ++cell_iterator;
            }
            times[i] = p_simulation_time->GetTime();
            i++;
        }
        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 8u);
    }

    void TestWithFixedG1GenerationalCellCycleModel()
    {
        // Simulation time is 6000 because we want to test that differentiated cells never divide.

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(6000.0, 1000);

        p_simulation_time->IncrementTimeOneStep();
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        //  Creating different types of cells with different cell-cycle models at SimulationTime = 6 hours
        FixedG1GenerationalCellCycleModel* p_stem_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_type);
        p_stem_cell->InitialiseCellCycleModel();

        // This test needs particular cell cycle times could also test for other cell-cycle models
        TS_ASSERT_DELTA(p_stem_model->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_stem_model->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_stem_model->GetSG2MDuration(), 10.0, 1e-12);

        UniformG1GenerationalCellCycleModel* p_stoch_model = new UniformG1GenerationalCellCycleModel();
        CellPtr p_stochastic_stem_cell(new Cell(p_healthy_state, p_stoch_model));
        p_stochastic_stem_cell->SetCellProliferativeType(p_type);
        p_stochastic_stem_cell->InitialiseCellCycleModel();

        FixedG1GenerationalCellCycleModel* p_diff_model = new FixedG1GenerationalCellCycleModel();
        p_diff_model->SetGeneration(6);
        CellPtr p_differentiated_cell(new Cell(p_healthy_state, p_diff_model));
        p_differentiated_cell->SetCellProliferativeType(p_diff_type);
        p_differentiated_cell->InitialiseCellCycleModel();

        UniformG1GenerationalCellCycleModel* p_stoch_diff_model = new UniformG1GenerationalCellCycleModel();
        p_stoch_diff_model->SetGeneration(6);
        CellPtr p_stochastic_differentiated_cell(new Cell(p_healthy_state, p_stoch_diff_model));
        p_stochastic_differentiated_cell->SetCellProliferativeType(p_diff_type);
        p_stochastic_differentiated_cell->InitialiseCellCycleModel();

        FixedG1GenerationalCellCycleModel* p_transit_model = new FixedG1GenerationalCellCycleModel();
        p_transit_model->SetGeneration(2);
        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        // SimulationTime = 6 hours
        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime = 18 hours
        TS_ASSERT_EQUALS(p_stochastic_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), true);

        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime = 24 hours
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_stochastic_stem_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime = 30 hours
        TS_ASSERT_EQUALS(p_stem_cell->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_stochastic_stem_cell->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep();

        // SimulationTime = 36 hours
        TS_ASSERT_EQUALS(p_stochastic_stem_cell->ReadyToDivide(), true);

        CellPtr p_daughter_cell1 = p_stem_cell->Divide();
        TS_ASSERT(typeid(p_daughter_cell1->GetCellCycleModel()) == typeid(p_stem_cell->GetCellCycleModel()));

        // Go to large time to ensure that differentiated cells can not divide
        for (unsigned i=0; i<990; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        TS_ASSERT_EQUALS(p_differentiated_cell->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_stochastic_differentiated_cell->ReadyToDivide(), false);
    }

    void TestStochasticCycleModel()
    {
        // Go up in steps of 0.01 to test stochasticity in cell-cycle models
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 5400);

        for (unsigned i=0; i<600; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }

        // Now at t=6.00
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
        CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_transit_model = new FixedG1GenerationalCellCycleModel();
        p_transit_model->SetGeneration(2);
        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        // This test needs particular cell cycle times
        TS_ASSERT_DELTA(p_transit_model->GetStemCellG1Duration(), 14.0, 1e-12);
        TS_ASSERT_DELTA(p_transit_model->GetTransitCellG1Duration(), 2.0, 1e-12);
        TS_ASSERT_DELTA(p_transit_model->GetSG2MDuration(), 10.0, 1e-12);

        for (unsigned i=0; i<1199; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }

        // Now at t = 17.99, cell is 11.99 old
        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), false);

        UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel;
        p_transit_cell->SetCellProliferativeType(p_transit_type);

        // This now resets the age of the cell to 0.0 so more time added in underneath
        p_transit_cell->SetCellCycleModel(p_cell_cycle_model);
        p_transit_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_transit_cell->GetCellCycleModel(), p_cell_cycle_model);
        for (unsigned i=0; i<1399; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), true);

        // Ensure transit cell divides
        while (!p_transit_cell->ReadyToDivide())
        {
            p_simulation_time->IncrementTimeOneStep();
        }
        CellPtr p_daughter_cell2 = p_transit_cell->Divide();
        TS_ASSERT(typeid(p_daughter_cell2->GetCellCycleModel()) == typeid(p_transit_cell->GetCellCycleModel()));
    }

    void Test0DBucketStochastic()
    {
        const double end_time = 70.0;
        const unsigned number_of_simulations = 1000;

        std::vector<CellPtr> cells;
        std::vector<CellPtr> newly_born;

        std::vector<unsigned> stem_cells(number_of_simulations);
        std::vector<unsigned> transit_cells(number_of_simulations);
        std::vector<unsigned> differentiated_cells(number_of_simulations);
        double stem_cell_mean = 0.0;
        double transit_cell_mean = 0.0;
        double differentiated_cell_mean = 0.0;

        for (unsigned simulation_number=0; simulation_number<number_of_simulations; simulation_number++)
        {
            SimulationTime::Destroy();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(70.0, 70);

            boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
            CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
            CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();

            UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();
            CellPtr p_stem_cell(new Cell(p_healthy_state, p_model));
            p_stem_cell->SetCellProliferativeType(p_type);
            p_stem_cell->InitialiseCellCycleModel();
            cells.push_back(p_stem_cell);

            // Produce the offspring of the cells
            std::vector<CellPtr>::iterator cell_iterator = cells.begin();

            while (p_simulation_time->GetTime()< end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                cell_iterator = cells.begin();
                while (cell_iterator < cells.end())
                {
                    if ((*cell_iterator)->ReadyToDivide())
                    {
                        newly_born.push_back((*cell_iterator)->Divide());
                    }
                    ++cell_iterator;
                }

                // Copy offspring in newly_born vector to cells vector
                cell_iterator = newly_born.begin();
                while (cell_iterator < newly_born.end())
                {
                    cells.push_back(*cell_iterator);
                    ++cell_iterator;
                }

                newly_born.clear();
            }
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                {
                    stem_cells[simulation_number]++;
                    stem_cell_mean++;
                }
                else if ((*cell_iterator)->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                {
                    transit_cells[simulation_number]++;
                    transit_cell_mean++;
                }
                else if ((*cell_iterator)->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
                {
                    differentiated_cells[simulation_number]++;
                    differentiated_cell_mean++;
                }
                ++cell_iterator;
            }
            cells.clear();
        }
        stem_cell_mean=stem_cell_mean/(double) number_of_simulations;
        transit_cell_mean=transit_cell_mean/(double) number_of_simulations;
        differentiated_cell_mean=differentiated_cell_mean/(double) number_of_simulations;

        TS_ASSERT_DELTA(stem_cell_mean, 1.0, 1e-12);
        TS_ASSERT_DELTA(transit_cell_mean, 2.0, 1.0);
        TS_ASSERT_DELTA(differentiated_cell_mean, 8.0, 1.0);
    }

    void TestInitialise0DBucket()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 60;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(60.0, num_steps);

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_stem_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_type);
        p_stem_cell->InitialiseCellCycleModel();
        cells.push_back(p_stem_cell);

        FixedG1GenerationalCellCycleModel* p_transit_model = new FixedG1GenerationalCellCycleModel();
        p_transit_model->SetGeneration(1);
        CellPtr p_transit_cell_1(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell_1->SetCellProliferativeType(p_transit_type);
        p_transit_cell_1->InitialiseCellCycleModel();
        cells.push_back(p_transit_cell_1);

        FixedG1GenerationalCellCycleModel* p_transit_model2 = new FixedG1GenerationalCellCycleModel();
        p_transit_model2->SetGeneration(2);
        CellPtr p_transit_cell_2(new Cell(p_healthy_state, p_transit_model2));
        p_transit_cell_2->SetCellProliferativeType(p_transit_type);
        p_transit_cell_2->InitialiseCellCycleModel();
        cells.push_back(p_transit_cell_2);

        FixedG1GenerationalCellCycleModel* p_transit_model3 = new FixedG1GenerationalCellCycleModel();
        p_transit_model3->SetGeneration(3);
        CellPtr p_transit_cell_3(new Cell(p_healthy_state, p_transit_model3));
        p_transit_cell_3->SetCellProliferativeType(p_transit_type);
        p_transit_cell_3->InitialiseCellCycleModel();
        cells.push_back(p_transit_cell_3);

        FixedG1GenerationalCellCycleModel* p_diff_model = new FixedG1GenerationalCellCycleModel();
        p_diff_model->SetGeneration(4);
        CellPtr p_differentiated_cell(new Cell(p_healthy_state, p_diff_model));
        p_differentiated_cell->SetCellProliferativeType(p_diff_type);
        p_differentiated_cell->InitialiseCellCycleModel();
        cells.push_back(p_differentiated_cell);

        std::vector<CellPtr> newly_born;
        std::vector<unsigned> stem_cells(num_steps);
        std::vector<unsigned> transit_cells(num_steps);
        std::vector<unsigned> differentiated_cells(num_steps);
        std::vector<double> times(num_steps);

        std::vector<CellPtr>::iterator cell_iterator;

        unsigned i = 0;
        while (!p_simulation_time->IsFinished())
        {
            p_simulation_time->IncrementTimeOneStep();

            // Produce the offspring of the cells
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                TS_ASSERT_EQUALS((*cell_iterator)->GetMutationState()->IsType<WildTypeCellMutationState>(), true);
                if ((*cell_iterator)->ReadyToDivide())
                {
                    newly_born.push_back((*cell_iterator)->Divide());
                }
                ++cell_iterator;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Count number of cells of each type
            cell_iterator = cells.begin();
            stem_cells[i] = 0;
            transit_cells[i] = 0;
            differentiated_cells[i] = 0;
            while (cell_iterator < cells.end())
            {
                if ((*cell_iterator)->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                {
                    stem_cells[i]++;
                }
                else if ((*cell_iterator)->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                {
                    transit_cells[i]++;
                }
                else if ((*cell_iterator)->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
                {
                    differentiated_cells[i]++;
                }
                ++cell_iterator;
            }

            times[i] = p_simulation_time->GetTime();
            i++;
        }

        TS_ASSERT_EQUALS(stem_cells[59], 1u);
        TS_ASSERT_EQUALS(transit_cells[59], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[59], 23u);
    }

    /*
     * We are checking that the CellPtrs work with the T&N cell-cycle models here
     * That division of wnt cells and stuff works OK.
     *
     * It checks that the cell division thing works nicely too.
     */
    void TestWithTysonNovakCellCycleModel()
    {
        double standard_tyson_duration = 1.242;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(200.0/60.0, num_steps+1);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel();
        CellPtr p_tn_cell(new Cell(p_healthy_state, p_cell_model));
        p_tn_cell->SetCellProliferativeType(p_transit_type);
        p_tn_cell->InitialiseCellCycleModel();
        p_cell_model->SetDt(0.1/60.0);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();
            if (time>standard_tyson_duration)
            {
                TS_ASSERT_EQUALS(p_tn_cell->ReadyToDivide(), true);
            }
            else
            {
                TS_ASSERT_EQUALS(p_tn_cell->ReadyToDivide(), false);
            }
        }

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_tn_cell->ReadyToDivide(), true);

        CellPtr p_tn_cell2 = p_tn_cell->Divide();

        double time_of_birth = p_tn_cell->GetBirthTime();
        double time_of_birth2 = p_tn_cell2->GetBirthTime();

        TS_ASSERT_DELTA(time_of_birth, time_of_birth2, 1e-9);

        for (unsigned i=0; i<num_steps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();
            bool result1 = p_tn_cell->ReadyToDivide();
            bool result2 = p_tn_cell2->ReadyToDivide();

            if (time >= standard_tyson_duration + time_of_birth)
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
    }

    void TestTysonNovakSteadyState()
    {
        // Keep dividing until we reach steady-state
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps=100000;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(20000.0/60.0, num_steps+1);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel();
        CellPtr p_tn_cell(new Cell(p_healthy_state, p_cell_model));
        p_tn_cell->SetCellProliferativeType(p_transit_type);
        p_tn_cell->InitialiseCellCycleModel();

        unsigned num_divisions = 0;

        while (!p_simulation_time->IsFinished())
        {
            while (!p_simulation_time->IsFinished() && !p_tn_cell->ReadyToDivide())
            {
                p_simulation_time->IncrementTimeOneStep();
            }
            if (p_tn_cell->ReadyToDivide())
            {
                p_tn_cell->Divide();
                ++num_divisions;
            }
        }
        TS_ASSERT_EQUALS(num_divisions, 268u);
    }

    void TestApoptosisAndDeath()
    {
        // We are going to start at t=0 and jump up in steps of 0.2
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.6, 3);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model,p_srn_model));
        p_cell->SetCellProliferativeType(p_transit_type);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), false);
        TS_ASSERT_EQUALS(p_cell->IsDead(), false);
        TS_ASSERT_THROWS_THIS(p_cell->GetTimeUntilDeath(),"Shouldn\'t be checking time until apoptosis as it isn\'t set");
        TS_ASSERT_DELTA(p_cell->GetApoptosisTime(), 0.25, 1e-6);

        p_simulation_time->IncrementTimeOneStep(); // t=0.2

        p_cell->StartApoptosis();
        TS_ASSERT_THROWS_THIS(p_cell->StartApoptosis(),"StartApoptosis() called when already undergoing apoptosis");

        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell->IsDead(), false);
        TS_ASSERT_DELTA(p_cell->GetTimeUntilDeath(),0.25,1e-12);

        // Check that we can copy a cell that has started apoptosis
        CellPtr p_cell2 = p_cell;

        p_simulation_time->IncrementTimeOneStep(); // t=0.4
        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell->IsDead(), false);
        TS_ASSERT_DELTA(p_cell->GetTimeUntilDeath(),0.05,1e-12);

        TS_ASSERT_EQUALS(p_cell2->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell2->IsDead(), false);
        TS_ASSERT_DELTA(p_cell2->GetTimeUntilDeath(),0.05,1e-12);

        p_simulation_time->IncrementTimeOneStep(); // t=0.6
        TS_ASSERT_EQUALS(p_cell->HasApoptosisBegun(), true);
        TS_ASSERT_EQUALS(p_cell->IsDead(), true);
    }

    void TestCantDivideIfUndergoingApoptosis()
    {
        // We are going to start at t=0 and jump up to t=25
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 1);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model,p_srn_model));
        p_cell->SetCellProliferativeType(p_transit_type);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();
        p_simulation_time->IncrementTimeOneStep(); // t=25

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
        p_cell->StartApoptosis();
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), false);
    }

    void Test0DBucketWithDeath()
    {
        double end_time = 92.0;
        unsigned num_time_steps = 92;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_time_steps);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_stem_cell(new Cell(p_healthy_state, p_cell_model));
        p_stem_cell->SetCellProliferativeType(p_type);
        p_stem_cell->InitialiseCellCycleModel();

        std::vector<CellPtr> cells;
        std::vector<CellPtr> newly_born;
        std::vector<unsigned> stem_cells(num_time_steps);
        std::vector<unsigned> transit_cells(num_time_steps);
        std::vector<unsigned> differentiated_cells(num_time_steps);
        std::vector<unsigned> dead_cells(num_time_steps);
        std::vector<double> times(num_time_steps);

        cells.push_back(p_stem_cell);
        std::vector<CellPtr>::iterator cell_iterator;

        unsigned i = 0;
        while (p_simulation_time->GetTime() < end_time)
        {
            // Produce the offspring of the cells

            p_simulation_time->IncrementTimeOneStep();
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (!(*cell_iterator)->IsDead())
                {
                    if ((*cell_iterator)->ReadyToDivide())
                    {
                        newly_born.push_back((*cell_iterator)->Divide());
                    }

                    if (((*cell_iterator)->GetAge() > 30))
                    {
                        (*cell_iterator)->StartApoptosis();
                    }
                }
                ++cell_iterator;
            }

            // Copy offspring in newly_born vector to cells vector
            cell_iterator = newly_born.begin();
            while (cell_iterator < newly_born.end())
            {
                cells.push_back(*cell_iterator);
                ++cell_iterator;
            }
            newly_born.clear();

            // Update cell counts
            cell_iterator = cells.begin();
            while (cell_iterator < cells.end())
            {
                if (!(*cell_iterator)->IsDead())
                {
                    if ((*cell_iterator)->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                    {
                        stem_cells[i]++;
                    }
                    else if ((*cell_iterator)->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                    {
                        transit_cells[i]++;
                    }
                    else if ((*cell_iterator)->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
                    {
                        differentiated_cells[i]++;
                    }
                }
                else
                {
                    dead_cells[i]++;
                }

                ++cell_iterator;
            }
            times[i] = p_simulation_time->GetTime();
            i++;
        }

        TS_ASSERT_EQUALS(stem_cells[num_time_steps-1], 1u);
        TS_ASSERT_EQUALS(transit_cells[num_time_steps-1], 2u);
        TS_ASSERT_EQUALS(differentiated_cells[num_time_steps-1], 8u);
        TS_ASSERT_EQUALS(dead_cells[num_time_steps-1], 8u);
    }

    void TestIsLogged()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 1);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell->IsLogged(), false);

        p_cell->SetLogged();

        TS_ASSERT_EQUALS(p_cell->IsLogged(), true);

        CellPtr p_copied_cell = p_cell;

        TS_ASSERT_EQUALS(p_copied_cell->IsLogged(), true);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_daughter_cell = p_cell->Divide();

        TS_ASSERT_EQUALS(p_cell->IsLogged(), true);
        TS_ASSERT_EQUALS(p_daughter_cell->IsLogged(), false);
    }

    void TestAncestors()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model, p_srn_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (2u));
        p_cell->SetAncestor(p_cell_ancestor);

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_cell2 = p_cell->Divide();

        TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);
        TS_ASSERT_EQUALS(p_cell2->GetAncestor(), 2u);

        //Coverage of SetAncesctor method
        TS_ASSERT_THROWS_THIS(p_cell->SetAncestor(p_healthy_state), "Attempting to give cell a cell ancestor which is not a CellAncestor");
    }

    void TestCellId()
    {
        // Resetting the Maximum cell Id to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        TS_ASSERT_EQUALS(p_cell->GetCellId(), 0u);

        CellPtr p_cell2 = p_cell->Divide();

        TS_ASSERT_EQUALS(p_cell->GetCellId(), 0u);
        TS_ASSERT_EQUALS(p_cell2->GetCellId(), 1u);
    }
};

#endif /*TESTCELL_HPP_*/
