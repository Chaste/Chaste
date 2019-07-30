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

#ifndef TESTCELLPROPERTYCOLLECTION_HPP_
#define TESTCELLPROPERTYCOLLECTION_HPP_

#include "CheckpointArchiveTypes.hpp"

#include <boost/shared_ptr.hpp>

#include "CellPropertyCollection.hpp"
#include "AbstractCellProperty.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

#include "OutputFileHandler.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

#define NEW_PROP(type, name) boost::shared_ptr<AbstractCellProperty> name(new type)

class TestCellPropertyCollection : public AbstractCellBasedTestSuite
{
public:
    void TestPropertyCollection()
    {
        CellPropertyCollection collection;
        TS_ASSERT_EQUALS(collection.GetSize(), 0u);

        // Test that the CellPropertyRegistry assigned to the CellPropertyCollection defaults to CellPropertyRegistry::Instance().
        TS_ASSERT_EQUALS(collection.GetCellPropertyRegistry(),CellPropertyRegistry::Instance());

        // Add some properties
        NEW_PROP(WildTypeCellMutationState, p_wt_mutation);

        TS_ASSERT_EQUALS(p_wt_mutation->GetIdentifier(), "WildTypeCellMutationState");

        collection.AddProperty(p_wt_mutation);
        NEW_PROP(ApcOneHitCellMutationState, p_apc1_mutation);
        collection.AddProperty(p_apc1_mutation);

        // Test we can't add the same *object* twice
        TS_ASSERT_THROWS_THIS(collection.AddProperty(p_wt_mutation),
                              "That property object is already in the collection.");
        NEW_PROP(WildTypeCellMutationState, p_wt_mutation_2);
        collection.AddProperty(p_wt_mutation_2);
        collection.RemoveProperty(p_wt_mutation_2);

        // Check the contents
        TS_ASSERT_EQUALS(collection.GetSize(), 2u);
        // ...by object
        TS_ASSERT_EQUALS(collection.HasProperty(p_wt_mutation), true);
        TS_ASSERT_EQUALS(collection.HasProperty(p_apc1_mutation), true);
        NEW_PROP(ApcOneHitCellMutationState, p_apc1_mutation_2);
        TS_ASSERT_EQUALS(collection.HasProperty(p_apc1_mutation_2), false);
        // ...by type
        TS_ASSERT_EQUALS(collection.HasProperty<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(collection.HasProperty<ApcOneHitCellMutationState>(), true);
        TS_ASSERT_EQUALS(collection.HasProperty<ApcTwoHitCellMutationState>(), false);
        // ..by subclass
        TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellProperty>(), true);
        TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellMutationState>(), true);
        //TS_ASSERT_EQUALS(!collection.HasProperty<AbstractCellMutationState>(), false); <-- This won't compile (yet)
        // ..by iteration
        for (CellPropertyCollection::Iterator it = collection.Begin(); it != collection.End(); ++it)
        {
            TS_ASSERT_EQUALS(collection.HasProperty(*it), true);
            TS_ASSERT((*it)->IsType<WildTypeCellMutationState>() || (*it)->IsType<ApcOneHitCellMutationState>());
            TS_ASSERT_EQUALS((*it)->IsSubType<AbstractCellMutationState>(), true);
        }

        // Remove property
        collection.RemoveProperty<WildTypeCellMutationState>();
        TS_ASSERT_EQUALS(collection.HasProperty<WildTypeCellMutationState>(), false);
        collection.RemoveProperty(p_apc1_mutation);
        TS_ASSERT_EQUALS(collection.HasProperty<ApcOneHitCellMutationState>(), false);
        TS_ASSERT_THROWS_THIS(collection.RemoveProperty<WildTypeCellMutationState>(),
                              "Collection does not contain the given property type.");
        TS_ASSERT_THROWS_THIS(collection.RemoveProperty(p_apc1_mutation),
                              "Collection does not contain the given property.");

        TS_ASSERT_EQUALS(collection.GetSize(), 0u);

        // Get matching properties
        collection.AddProperty(p_wt_mutation);
        collection.AddProperty(p_apc1_mutation);
        CellPropertyCollection mutations = collection.GetPropertiesType<AbstractCellMutationState>();
        TS_ASSERT_EQUALS(mutations.GetSize(), 2u);
        CellPropertyCollection::Iterator it = mutations.Begin();
        TS_ASSERT_EQUALS(collection.HasProperty(*it), true);
        TS_ASSERT_EQUALS((*it)->IsSubType<AbstractCellMutationState>(), true);
        TS_ASSERT( (*it)->IsType<WildTypeCellMutationState>() || (*(++it))->IsType<WildTypeCellMutationState>() );

        CellPropertyCollection wild_types = collection.GetProperties<WildTypeCellMutationState>();
        TS_ASSERT_EQUALS(wild_types.GetSize(), 1u);
        it = wild_types.Begin();
        TS_ASSERT_EQUALS((*it)->IsType<WildTypeCellMutationState>(), true);
        TS_ASSERT( *it == wild_types.GetProperty() );
        TS_ASSERT_THROWS_THIS(collection.GetProperty(),
                              "Can only call GetProperty on a collection of size 1.");
    }

    void TestArchiveCellPropertyCollection()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "properties.arch";

        // Archive a cell property collection
        {
            // Create a cell property collection
            CellPropertyCollection collection;

            NEW_PROP(BetaCateninOneHitCellMutationState, p_bcat1_mutation);
            collection.AddProperty(p_bcat1_mutation);
            NEW_PROP(ApcOneHitCellMutationState, p_apc1_mutation);
            collection.AddProperty(p_apc1_mutation);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            output_arch << static_cast<const CellPropertyCollection>(collection);
        }

        // Restore cell property collection
        {
            CellPropertyCollection collection;

            // Restore the cell property collection
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> collection;

            // Test that the collection was archived correctly
            TS_ASSERT_EQUALS(collection.GetSize(), 2u);

            TS_ASSERT_EQUALS(collection.HasProperty<BetaCateninOneHitCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcOneHitCellMutationState>(), true);
            TS_ASSERT_EQUALS(collection.HasProperty<ApcTwoHitCellMutationState>(), false);

            NEW_PROP(ApcOneHitCellMutationState, p_apc1_mutation_2);
            TS_ASSERT_EQUALS(collection.HasProperty(p_apc1_mutation_2), false);

            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellProperty>(), true);
            TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellMutationState>(), true);

            for (CellPropertyCollection::Iterator it = collection.Begin(); it != collection.End(); ++it)
            {
                TS_ASSERT_EQUALS(collection.HasProperty(*it), true);
                TS_ASSERT((*it)->IsType<BetaCateninOneHitCellMutationState>() || (*it)->IsType<ApcOneHitCellMutationState>());
                TS_ASSERT_EQUALS((*it)->IsSubType<AbstractCellMutationState>(), true);
            }
        }
    }
};

#endif /* TESTCELLPROPERTYCOLLECTION_HPP_ */
