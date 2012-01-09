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

#ifndef TESTCELLLABEL_HPP_
#define TESTCELLLABEL_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"

#include <boost/shared_ptr.hpp>

#include "CellLabel.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "OutputFileHandler.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"

class TestCellLabel : public AbstractCellBasedTestSuite
{
public:

    void TestCellLabelMethods() throw (Exception)
    {
        MAKE_PTR(CellLabel, p_label);
        TS_ASSERT_EQUALS(p_label->GetCellCount(), 0u);
        p_label->IncrementCellCount();
        TS_ASSERT_EQUALS(p_label->GetCellCount(), 1u);
        p_label->DecrementCellCount();
        TS_ASSERT_EQUALS(p_label->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_label->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");
        TS_ASSERT_EQUALS(p_label->GetColour(), 5u);

        TS_ASSERT_EQUALS(p_label->IsType<CellLabel>(), true);
        TS_ASSERT_EQUALS(p_label->IsType<ApcOneHitCellMutationState>(), false);

        MAKE_PTR(CellLabel, p_other_label);
        TS_ASSERT_EQUALS(p_other_label->IsSame(p_label.get()), true);
        TS_ASSERT_EQUALS(p_label->IsSame(p_other_label), true);

        // Check that const-ness doesn't matter
        TS_ASSERT_EQUALS(p_other_label->IsType<const CellLabel>(), true);
        const CellLabel const_label;
        TS_ASSERT_EQUALS(p_other_label->IsSame(&const_label), true);
        TS_ASSERT_EQUALS(const_label.IsSame(p_other_label), true);
        TS_ASSERT(const_label.IsSame(p_other_label.get()));
    }

    void TestCellPropertyRegistry() throw(Exception)
    {
        boost::shared_ptr<AbstractCellProperty> p_property1(
                CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_property2(
                CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        TS_ASSERT(p_property1 == p_property2);
        TS_ASSERT_EQUALS(p_property1->IsType<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_property1->IsSubType<AbstractCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_property1->IsType<CellLabel>(), false);

        std::vector<boost::shared_ptr<AbstractCellProperty> > properties =
            CellPropertyRegistry::Instance()->rGetAllCellProperties();
        TS_ASSERT_EQUALS(properties.size(), 1u);
        TS_ASSERT(properties[0] == p_property1);

        CellPropertyRegistry::Instance()->Clear();
        properties = CellPropertyRegistry::Instance()->rGetAllCellProperties();
        TS_ASSERT_EQUALS(properties.size(), 0u);

        // The taking-ownership functionality
        p_property1 = CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>();
        CellPropertyRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>();
        TS_ASSERT_EQUALS(CellPropertyRegistry::Instance()->rGetAllCellProperties().size(), 2u);
        CellPropertyRegistry* p_instance = CellPropertyRegistry::Instance();
        CellPropertyRegistry* p_registry = CellPropertyRegistry::Instance()->TakeOwnership();
        TS_ASSERT_EQUALS(p_instance, p_registry);
        TS_ASSERT_EQUALS(p_registry->rGetAllCellProperties().size(), 2u);
        TS_ASSERT_DIFFERS(CellPropertyRegistry::Instance(), p_registry);
        TS_ASSERT_EQUALS(CellPropertyRegistry::Instance()->rGetAllCellProperties().size(), 0u);
        TS_ASSERT_EQUALS(p_registry->rGetAllCellProperties().size(), 2u);

        // Tidy up
        delete p_registry;
    }

    void TestCellPropertyOrdering() throw(Exception)
    {
        CellPropertyRegistry* p_instance = CellPropertyRegistry::Instance();
        p_instance->Clear();
        p_instance->Get<ApcOneHitCellMutationState>();

        std::vector<boost::shared_ptr<AbstractCellProperty> > property_ordering;
        property_ordering.push_back(p_instance->Get<WildTypeCellMutationState>());
        property_ordering.push_back(p_instance->Get<ApcOneHitCellMutationState>());
        property_ordering.push_back(p_instance->Get<ApoptoticCellProperty>());
        property_ordering.push_back(p_instance->Get<CellLabel>());

        TS_ASSERT_EQUALS(p_instance->HasOrderingBeenSpecified(), false);
        p_instance->SpecifyOrdering(property_ordering);
        TS_ASSERT_EQUALS(p_instance->HasOrderingBeenSpecified(), true);

        TS_ASSERT_THROWS_THIS(p_instance->SpecifyOrdering(property_ordering), "An ordering has already been specified.");

        std::vector<boost::shared_ptr<AbstractCellProperty> > properties = p_instance->rGetAllCellProperties();
        TS_ASSERT_EQUALS(properties.size(), 4u);
        TS_ASSERT_EQUALS(properties[0]->IsType<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(properties[1]->IsType<ApcOneHitCellMutationState>(), true);
        TS_ASSERT_EQUALS(properties[2]->IsType<ApoptoticCellProperty>(), true);
        TS_ASSERT_EQUALS(properties[3]->IsType<CellLabel>(), true);
        p_instance->Clear();
    }

    void TestArchiveCellLabel() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "label.arch";

        // Archive a cell label
        {
            CellLabel* p_label = new CellLabel();
            p_label->IncrementCellCount();

            TS_ASSERT_EQUALS(p_label->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_label->GetColour(), 5u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_label = p_label;
            output_arch << p_const_label;

            delete p_label;
        }

        // Restore cell label
        {
            AbstractCellProperty* p_label;

            // Restore the cell label
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_label;

            CellLabel* p_real_label = dynamic_cast<CellLabel*>(p_label);
            TS_ASSERT(p_real_label != NULL);

            TS_ASSERT_EQUALS(p_real_label->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_real_label->GetColour(), 5u);

            // Tidy up
            delete p_label;
        }
    }

    void TestArchiveCellPropertyRegistry() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "property.arch";

        // Save
        {
            boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

            TS_ASSERT_EQUALS(p_label->IsType<CellLabel>(), true);
            TS_ASSERT_EQUALS(p_label->IsSubType<AbstractCellProperty>(), true);
            TS_ASSERT_EQUALS(p_label->IsSubType<AbstractCellMutationState>(), false);

            const CellPropertyRegistry* const p_registry = CellPropertyRegistry::Instance()->TakeOwnership();

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the registry to the archive
            output_arch << p_registry;

            // Tidy up
            delete p_registry;
        }

        // Restore boost::shared_ptr to cell property
        {
            CellPropertyRegistry::Instance()->Clear();

            // Initialize a registry
            CellPropertyRegistry* p_registry;

            // Restore the registry state
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_registry;

            boost::shared_ptr<AbstractCellProperty> p_label = p_registry->Get<CellLabel>();

            TS_ASSERT_EQUALS(p_label->IsType<CellLabel>(), true);
            TS_ASSERT_EQUALS(p_label->IsSubType<AbstractCellProperty>(), true);
            TS_ASSERT_EQUALS(p_label->IsSubType<AbstractCellMutationState>(), false);

            // Tidy up
            delete p_registry;
        }
    }
};

#endif /* TESTCELLLABEL_HPP_ */
