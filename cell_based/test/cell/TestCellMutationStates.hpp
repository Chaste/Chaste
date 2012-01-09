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

#ifndef TESTCELLMUTATIONSTATES_HPP_
#define TESTCELLMUTATIONSTATES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

#include "CellPropertyRegistry.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

class TestCellMutationStates : public AbstractCellBasedTestSuite
{
public:

    void TestCellMutationStateMethods() throw(Exception)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        p_state->IncrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
        p_state->DecrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_state->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");
        TS_ASSERT_EQUALS(p_state->GetColour(), 0u);

        TS_ASSERT_EQUALS(p_state->IsType<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_state->IsType<ApcOneHitCellMutationState>(), false);

        MAKE_PTR(WildTypeCellMutationState, p_wt_state);
        MAKE_PTR(ApcTwoHitCellMutationState, p_apc2_state);

        TS_ASSERT(p_wt_state->IsSame(p_state.get()));
        TS_ASSERT(p_state->IsSame(p_wt_state));
        TS_ASSERT_EQUALS(p_wt_state->IsSame(p_apc2_state.get()), false);
        TS_ASSERT_EQUALS(p_apc2_state->IsSame(p_wt_state), false);

        // Check that const-ness doesn't matter
        TS_ASSERT(p_wt_state->IsType<const WildTypeCellMutationState>());
        const WildTypeCellMutationState const_wt_state;
        TS_ASSERT(p_wt_state->IsSame(&const_wt_state));
        TS_ASSERT(const_wt_state.IsSame(p_wt_state));
        TS_ASSERT(const_wt_state.IsSame(p_wt_state.get()));
    }

    void TestArchiveCellMutationState() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation.arch";

        // Archive a mutation state
        {
            ApcOneHitCellMutationState* p_state = new ApcOneHitCellMutationState();
            p_state->IncrementCellCount();

            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_state->GetColour(), 3u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_state;
            output_arch << p_const_state;

            delete p_state;
        }

        // Restore mutation state
        {
            AbstractCellProperty* p_state;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_state;

            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);

            ApcOneHitCellMutationState* p_real_state = dynamic_cast<ApcOneHitCellMutationState*>(p_state);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 3u);

            // Tidy up
            delete p_state;
        }
    }
};

#endif /* TESTCELLMUTATIONSTATES_HPP_ */
