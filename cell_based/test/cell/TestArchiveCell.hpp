/*

Copyright (c) 2005-2012, University of Oxford.
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
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "SmartPointers.hpp"

/*
 * This test is seperate from TestCell.hpp to avoid strange errors with the
 * intel compiler - see #1569
 */
class TestArchiveCell: public AbstractCellBasedTestSuite
{
public:

    void TestArchivingOfCell() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell.arch";

        // Archive a cell
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // Create mutation state
            boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

            // Create cell-cycle model
            FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
            p_cell_model->SetCellProliferativeType(STEM);

            // Create cell property collection
            CellPropertyCollection collection;
            MAKE_PTR(CellLabel, p_label);
            collection.AddProperty(p_label);

            // Create cell
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model, false, collection));
            p_cell->InitialiseCellCycleModel();
            p_simulation_time->IncrementTimeOneStep();

            TS_ASSERT_EQUALS(p_cell->GetAge(), 0.5);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), UNSIGNED_UNSET);

            // Set ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (2u));
            p_cell->SetAncestor(p_cell_ancestor);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);

            // Check properties set correctly
            CellPropertyCollection& final_collection = p_cell->rGetCellPropertyCollection();
            TS_ASSERT_EQUALS(final_collection.GetSize(), 4u);
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
                TS_ASSERT((*it)->IsType<WildTypeCellMutationState>() || (*it)->IsType<CellLabel>() || (*it)->IsType<CellAncestor>() || (*it)->IsType<CellId>());
            }

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            CellPtr const p_const_cell = p_cell;

            // Write the cell to the archive
            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            output_arch << p_const_cell;

            // Tidy up
            SimulationTime::Destroy();
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
            TS_ASSERT_EQUALS(static_cast<FixedDurationGenerationBasedCellCycleModel*>(p_cell->GetCellCycleModel())->GetGeneration(), 0u);
            TS_ASSERT_EQUALS(p_cell->GetCellCycleModel()->GetCellProliferativeType(), STEM);

            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
            TS_ASSERT_EQUALS(p_model->GetCell(), p_cell);

            CellPropertyCollection& collection = p_cell->rGetCellPropertyCollection();
            TS_ASSERT_EQUALS(collection.GetSize(), 4u);
            TS_ASSERT_EQUALS(collection.HasProperty<WildTypeCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcOneHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcTwoHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(collection.HasProperty<CellLabel>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellProperty>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<CellAncestor>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<CellId>(), true);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);

            for (CellPropertyCollection::Iterator it = collection.Begin(); it != collection.End(); ++it)
            {
                TS_ASSERT_EQUALS(collection.HasProperty(*it), true);
                TS_ASSERT((*it)->IsType<WildTypeCellMutationState>() || (*it)->IsType<CellLabel>() || (*it)->IsType<CellAncestor>() || (*it)->IsType<CellId>());
            }
        }
    }

//    void TestFailingArchivingTest() throw(Exception)
//    {
//    	OutputFileHandler handler("archive", false);
//		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell.arch";
//
//		// Create and archive a cell cycle model
//		{
//			SimulationTime* p_simulation_time = SimulationTime::Instance();
//			p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
//
//			// Create cell-cycle model
//			FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
//			p_cell_model->SetCellProliferativeType(STEM);
//			p_cell_model->SetMDuration(60.0);
//
//			// Create an output archive
//			std::ofstream ofs(archive_filename.c_str());
//			boost::archive::text_oarchive output_arch(ofs);
//
//			// Write the cell cycle model to the archive
//			output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
//			output_arch << p_cell_model;
//
//			// Tidy up
//			SimulationTime::Destroy();
//		}
//
//		// load the model.
//		{
//			SimulationTime* p_simulation_time = SimulationTime::Instance();
//			p_simulation_time->SetStartTime(1.0);
//			p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 1); // will be restored
//
//			AbstractCellCycleModel* p_model;
//
//			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
//			boost::archive::text_iarchive input_arch(ifs);
//
//			input_arch >> *p_simulation_time;
//			input_arch >> p_model;
//		}
//    }
};

#endif /*TESTARCHIVECELL_HPP_*/
