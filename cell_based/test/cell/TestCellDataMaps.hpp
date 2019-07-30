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

#ifndef TESTCELLDATAMAPS_HPP_
#define TESTCELLDATAMAPS_HPP_

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <iostream>

#include "PetscVecTools.hpp"
#include "ReplicatableVector.hpp"

#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellVecData.hpp"
#include "CellData.hpp"
#include "SmartPointers.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCellDataMaps : public AbstractCellBasedTestSuite
{
public:

    void TestCellData()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();

        //Before adding the CellData to the cell
        TS_ASSERT_THROWS_NOTHING(p_cell->GetCellData());

        p_cell->GetCellData()->SetItem("something", 1.0);
        p_cell->GetCellData()->SetItem("some other thing", 2.0);

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_cell2 = p_cell->Divide();

        CellPropertyCollection cell_data_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellData>();
        boost::shared_ptr<CellData> p_parentcell_data = boost::static_pointer_cast<CellData>(cell_data_collection.GetProperty());
        p_parentcell_data->SetItem("something", 3.0);

        TS_ASSERT_EQUALS(p_cell->GetCellData()->GetItem("something"), 3.0);
        TS_ASSERT_EQUALS(p_cell->GetCellData()->GetItem("some other thing"), 2.0);

        CellPropertyCollection cell2_data_collection = p_cell2->rGetCellPropertyCollection().GetPropertiesType<CellData>();
        boost::shared_ptr<CellData> p_daughtercell_data = boost::static_pointer_cast<CellData>(cell2_data_collection.GetProperty());

        TS_ASSERT_EQUALS(p_daughtercell_data->GetItem("something"), 1.0);
        TS_ASSERT_EQUALS(p_daughtercell_data->GetItem("some other thing"), 2.0);
    }

    void TestCellVecData()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();

        // No CellVecData is available by default
        TS_ASSERT(!p_cell->HasCellVecData());

        // Add CellVecData
        boost::shared_ptr<AbstractCellProperty> p_vec_data(CellPropertyRegistry::Instance()->Get<CellVecData>());
        p_cell->AddCellProperty(p_vec_data);
        TS_ASSERT(p_cell->HasCellVecData());

        CellPropertyCollection parent_cell_property_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellVecData>();
        boost::shared_ptr<CellVecData> p_parent_cell_vec_data = boost::static_pointer_cast<CellVecData>(parent_cell_property_collection.GetProperty());

        Vec vec_item_1 = PetscTools::CreateAndSetVec(1, 17.3);
        std::string i1 = "item 1";
        p_parent_cell_vec_data->SetItem(i1, vec_item_1);

        Vec vec_item_2 = PetscTools::CreateAndSetVec(1, 3.58);
        p_parent_cell_vec_data->SetItem("item 2", vec_item_2);

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_daughter_cell = p_cell->Divide();

        Vec vec_replace_item_1 = PetscTools::CreateAndSetVec(1, -8.54);

        CellPropertyCollection parent_cell_property_collection_now = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellVecData>();
        boost::shared_ptr<CellVecData> p_parent_cell_vec_data_now = boost::static_pointer_cast<CellVecData>(parent_cell_property_collection_now.GetProperty());

        p_parent_cell_vec_data_now->SetItem("item 1", vec_replace_item_1);

        ReplicatableVector rep_item_1_parent(p_parent_cell_vec_data_now->GetItem("item 1"));
        TS_ASSERT_EQUALS(rep_item_1_parent[0], -8.54);

        Vec vec_item2_parent = p_parent_cell_vec_data_now->GetItem("item 2");
        ReplicatableVector rep_item_2_parent(vec_item2_parent);
        TS_ASSERT_EQUALS(rep_item_2_parent[0], 3.58);

        // We modify an entry in the Vec stored in the parent cell and we check a few lines below that it doesn't propagate to the daughter cell
        PetscVecTools::SetElement(vec_item2_parent, 0, -1.1);
        PetscVecTools::Finalise(vec_item2_parent);

        TS_ASSERT_THROWS_THIS(p_parent_cell_vec_data_now->GetItem("item 3"), "The item item 3 is not stored");

        CellPropertyCollection daughter_cell_property_collection = p_daughter_cell->rGetCellPropertyCollection().GetPropertiesType<CellVecData>();
        boost::shared_ptr<CellVecData> p_daughter_cell_vec_data = boost::static_pointer_cast<CellVecData>(daughter_cell_property_collection.GetProperty());

        ReplicatableVector rep_item_1_daughter(p_daughter_cell_vec_data->GetItem("item 1"));
        ReplicatableVector rep_item_2_daughter(p_daughter_cell_vec_data->GetItem("item 2"));
        TS_ASSERT_EQUALS(rep_item_1_daughter[0], 17.3);
        TS_ASSERT_EQUALS(rep_item_2_daughter[0], 3.58);

        // Free PETSc Vecs which were created in this test
        PetscTools::Destroy(vec_item_1);
        PetscTools::Destroy(vec_item_2);
        PetscTools::Destroy(vec_replace_item_1);
        // Vec vec_item2_parent isn't created here -- it's just pointing to an existing Vec
    }
};

#endif /*TESTCELLDATAMAPS_HPP_*/
