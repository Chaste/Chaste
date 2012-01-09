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

#ifndef TESTCELLPROPERTYCOLLECTION_HPP_
#define TESTCELLPROPERTYCOLLECTION_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"

#include <boost/shared_ptr.hpp>

#include "CellPropertyCollection.hpp"
#include "AbstractCellProperty.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

#include "OutputFileHandler.hpp"

#define NEW_PROP(type, name) boost::shared_ptr<AbstractCellProperty> name(new type)

class TestCellPropertyCollection : public AbstractCellBasedTestSuite
{
public:
    void TestPropertyCollection() throw (Exception)
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

    void TestArchiveCellPropertyCollection() throw (Exception)
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
