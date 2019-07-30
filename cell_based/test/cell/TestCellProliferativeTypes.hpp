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

#ifndef TESTCELLPROLIFERATIVETYPES_HPP_
#define TESTCELLPROLIFERATIVETYPES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellProliferativeType.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellPropertyRegistry.hpp"

#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCellProliferativeTypes : public AbstractCellBasedTestSuite
{
public:

    void TestCellProliferativeTypeMethods()
    {
        MAKE_PTR(StemCellProliferativeType, p_type);
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 0u);

        p_type->IncrementCellCount();
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

        p_type->DecrementCellCount();
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 0u);

        TS_ASSERT_THROWS_THIS(p_type->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");
        TS_ASSERT_EQUALS(p_type->GetColour(), 0u);

        TS_ASSERT_EQUALS(p_type->IsType<StemCellProliferativeType>(), true);
        TS_ASSERT_EQUALS(p_type->IsType<TransitCellProliferativeType>(), false);

        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        TS_ASSERT(p_stem_type->IsSame(p_type.get()));
        TS_ASSERT(p_type->IsSame(p_stem_type));

        TS_ASSERT_EQUALS(p_stem_type->IsSame(p_transit_type.get()), false);
        TS_ASSERT_EQUALS(p_transit_type->IsSame(p_stem_type), false);

        // Check that const-ness doesn't matter
        TS_ASSERT(p_stem_type->IsType<const StemCellProliferativeType>());
        const StemCellProliferativeType const_stem_type;

        TS_ASSERT(p_stem_type->IsSame(&const_stem_type));
        TS_ASSERT(const_stem_type.IsSame(p_stem_type));
        TS_ASSERT(const_stem_type.IsSame(p_stem_type.get()));
    }

    void TestArchiveStemCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "StemCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            StemCellProliferativeType* p_type = new StemCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 0u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            StemCellProliferativeType* p_real_state = dynamic_cast<StemCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 0u);

            // Tidy up
            delete p_type;
        }
    }

    void TestArchiveTransitCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "TransitCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            TransitCellProliferativeType* p_type = new TransitCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 1u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            TransitCellProliferativeType* p_real_state = dynamic_cast<TransitCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 1u);

            // Tidy up
            delete p_type;
        }
    }

    void TestArchiveDifferentiatedCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DifferentiatedCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            DifferentiatedCellProliferativeType* p_type = new DifferentiatedCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 2u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            DifferentiatedCellProliferativeType* p_real_state = dynamic_cast<DifferentiatedCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 2u);

            // Tidy up
            delete p_type;
        }
    }

    void TestArchiveDefaultCellProliferativeType()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DefaultCellProliferativeType.arch";

        // Archive a cell proliferative type
        {
            DefaultCellProliferativeType* p_type = new DefaultCellProliferativeType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 0u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell proliferative type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            DefaultCellProliferativeType* p_real_state = dynamic_cast<DefaultCellProliferativeType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 0u);

            // Tidy up
            delete p_type;
        }
    }
};

#endif /*TESTCELLPROLIFERATIVETYPES_HPP_*/
