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

#ifndef TESTCELLMUTATIONSTATES_HPP_
#define TESTCELLMUTATIONSTATES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

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

#include "FakePetscSetup.hpp"

class TestCellMutationStates : public AbstractCellBasedTestSuite
{
public:

    void TestCellMutationStateMethods()
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

    void TestArchiveCellMutationState()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell_mutation_state.arch";

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
