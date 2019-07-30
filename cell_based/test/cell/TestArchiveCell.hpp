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

#ifndef TESTARCHIVECELL_HPP_
#define TESTARCHIVECELL_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <iostream>

#include "OutputFileHandler.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "Goldbeter1991SrnModel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "CellId.hpp"
#include "CellAncestor.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "ReplicatableVector.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * This test is separate from TestCell.hpp to avoid strange errors with the
 * intel compiler - see #1569
 */
class TestArchiveCell: public AbstractCellBasedTestSuite
{
public:

    void TestArchivingOfCell()
    {
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell.arch";

        // Archive a cell
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // Create mutation state
            boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

            // Create cell-cycle model
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();

            // Create SRN
            Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();

            // Create cell property collection
            CellPropertyCollection collection;
            MAKE_PTR(CellLabel, p_label);
            collection.AddProperty(p_label);

            // Create cell
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model, p_srn_model,  false, collection));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            p_simulation_time->IncrementTimeOneStep();

            TS_ASSERT_EQUALS(p_cell->GetAge(), 0.5);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), UNSIGNED_UNSET);

            // Set ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (2u));
            p_cell->SetAncestor(p_cell_ancestor);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);

            // Set CellVecData with some actual content
            boost::shared_ptr<AbstractCellProperty> p_vec_data(CellPropertyRegistry::Instance()->Get<CellVecData>());
            p_cell->AddCellProperty(p_vec_data);
            Vec item_1 = PetscTools::CreateAndSetVec(2, -17.3); // <-17.3, -17.3>
            p_cell->GetCellVecData()->SetItem("item 1", item_1);

            // Check properties set correctly
            CellPropertyCollection& final_collection = p_cell->rGetCellPropertyCollection();
            TS_ASSERT_EQUALS(final_collection.GetSize(), 7u);
            TS_ASSERT_EQUALS(final_collection.HasProperty<WildTypeCellMutationState>(), true);
            TS_ASSERT_EQUALS(final_collection.HasProperty<ApcOneHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(final_collection.HasProperty<ApcTwoHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(final_collection.HasProperty<CellLabel>(), true);
            TS_ASSERT_EQUALS(final_collection.HasPropertyType<AbstractCellProperty>(), true);
            TS_ASSERT_EQUALS(final_collection.HasPropertyType<AbstractCellMutationState>(), true);
            TS_ASSERT_EQUALS(final_collection.HasPropertyType<CellAncestor>(), true);
            TS_ASSERT_EQUALS(final_collection.HasPropertyType<CellId>(), true);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);

            for (CellPropertyCollection::Iterator it = final_collection.Begin(); it != final_collection.End(); ++it)
            {
                TS_ASSERT_EQUALS(final_collection.HasProperty(*it), true);

                bool is_wildtype = (*it)->IsType<WildTypeCellMutationState>();
                bool is_label = (*it)->IsType<CellLabel>();
                bool is_ancestor = (*it)->IsType<CellAncestor>();
                bool is_cellid = (*it)->IsType<CellId>();
                bool is_data = (*it)->IsType<CellData>();
                bool is_vec_data = (*it)->IsType<CellVecData>();
                bool is_stem = (*it)->IsType<StemCellProliferativeType>();

                bool is_any_of_above = is_wildtype || is_label || is_ancestor || is_cellid || is_data || is_vec_data || is_stem;
                TS_ASSERT_EQUALS(is_any_of_above, true);
            }

            // Create another cell
            boost::shared_ptr<AbstractCellProperty> p_another_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            FixedG1GenerationalCellCycleModel* p_another_cell_model = new FixedG1GenerationalCellCycleModel();
            Goldbeter1991SrnModel* p_another_srn_model = new Goldbeter1991SrnModel();
            CellPtr p_another_cell(new Cell(p_another_healthy_state, p_another_cell_model, p_another_srn_model, false, collection));
            boost::shared_ptr<AbstractCellProperty> p_another_vec_data(new CellVecData);
            p_another_cell->AddCellProperty(p_another_vec_data);
            TS_ASSERT_EQUALS(p_cell->GetCellVecData()->GetNumItems(), 1u);
            TS_ASSERT_EQUALS(p_another_cell->GetCellVecData()->GetNumItems(), 0u);
            Vec another_item_1 = PetscTools::CreateAndSetVec(2, 42.0); // <42, 42>

            p_another_cell->GetCellVecData()->SetItem("item 1", another_item_1);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            CellPtr const p_const_cell = p_cell;

            // Write the cell to the archive
            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            output_arch << p_const_cell;
            // Write the second cell also
            CellPtr const p_another_const_cell = p_another_cell;
            output_arch << p_another_const_cell;

            // Tidy up
            SimulationTime::Destroy();
            PetscTools::Destroy(item_1);
            PetscTools::Destroy(another_item_1);
        }

        // Restore CellPtr
        {
            // Need to set up time to initialize a cell
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(1.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 1); // will be restored

            // Initialize a cell

            CellPtr p_cell;

            // Restore the cell
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> *p_simulation_time;
            input_arch >> p_cell;

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.5);
            TS_ASSERT_EQUALS(p_simulation_time->GetTimeStep(), 0.5);

            TS_ASSERT_EQUALS(p_cell->GetAge(), 0.5);
            TS_ASSERT_EQUALS(static_cast<FixedG1GenerationalCellCycleModel*>(p_cell->GetCellCycleModel())->GetGeneration(), 0u);
            TS_ASSERT(dynamic_cast<Goldbeter1991SrnModel*>(p_cell->GetSrnModel()));
            TS_ASSERT_EQUALS(p_cell->GetCellProliferativeType()->IsType<StemCellProliferativeType>(), true);

            AbstractCellCycleModel* p_cc_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_cc_model->GetCell(), p_cell);

            AbstractSrnModel* p_srn_model = p_cell->GetSrnModel();
            TS_ASSERT_EQUALS(p_srn_model->GetCell(), p_cell);

            CellPropertyCollection& collection = p_cell->rGetCellPropertyCollection();
            TS_ASSERT_EQUALS(collection.GetSize(), 7u);
            TS_ASSERT_EQUALS(collection.HasProperty<WildTypeCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcOneHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcTwoHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(collection.HasProperty<CellLabel>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellProperty>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<CellAncestor>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<CellId>(), true);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);

            // Check explicitly for CellVecData as it is not available by default as CellData
            TS_ASSERT_EQUALS(collection.HasPropertyType<CellVecData>(), true);

            // Check that the Vec stored in CellVecData was unarchived correctly
            boost::shared_ptr<CellVecData> p_cell_vec_data = boost::static_pointer_cast<CellVecData>(collection.GetPropertiesType<CellVecData>().GetProperty());
            PetscInt vec_size;
            VecGetSize(p_cell_vec_data->GetItem("item 1"), &vec_size);
            TS_ASSERT_EQUALS(vec_size, 2);
            ReplicatableVector rep_item_1(p_cell_vec_data->GetItem("item 1"));
            TS_ASSERT_DELTA(rep_item_1[0], -17.3, 2e-14);

            for (CellPropertyCollection::Iterator it = collection.Begin(); it != collection.End(); ++it)
            {
                TS_ASSERT_EQUALS(collection.HasProperty(*it), true);

                bool is_wildtype = (*it)->IsType<WildTypeCellMutationState>();
                bool is_label = (*it)->IsType<CellLabel>();
                bool is_ancestor = (*it)->IsType<CellAncestor>();
                bool is_cellid = (*it)->IsType<CellId>();
                bool is_data = (*it)->IsType<CellData>();
                bool is_vec_data = (*it)->IsType<CellVecData>();
                bool is_stem = (*it)->IsType<StemCellProliferativeType>();

                bool is_any_of_above = is_wildtype || is_label || is_ancestor || is_cellid || is_data || is_vec_data || is_stem;
                TS_ASSERT_EQUALS(is_any_of_above, true);
            }

            // Try another cell
            CellPtr p_another_cell;
            input_arch >> p_another_cell;
            ReplicatableVector rep_another_item_1(p_another_cell->GetCellVecData()->GetItem("item 1"));
            TS_ASSERT_DELTA(rep_another_item_1[0], 42.0, 2e-14);
        }
    }
};

#endif /*TESTARCHIVECELL_HPP_*/
