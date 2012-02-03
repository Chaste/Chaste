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

            // set ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (2u));
            p_cell->SetAncestor(p_cell_ancestor);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);

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
            TS_ASSERT_EQUALS(collection.GetSize(), 3u);
            TS_ASSERT_EQUALS(collection.HasProperty<WildTypeCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcOneHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcTwoHitCellMutationState>(), false);
            TS_ASSERT_EQUALS(collection.HasProperty<CellLabel>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellProperty>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<CellAncestor>(), true);
            TS_ASSERT_EQUALS(p_cell->GetAncestor(), 2u);

            for (CellPropertyCollection::Iterator it = collection.Begin(); it != collection.End(); ++it)
            {
                TS_ASSERT_EQUALS(collection.HasProperty(*it), true);
                TS_ASSERT((*it)->IsType<WildTypeCellMutationState>() || (*it)->IsType<CellLabel>() || (*it)->IsType<CellAncestor>());
            }
        }
    }
};

#endif /*TESTARCHIVECELL_HPP_*/
