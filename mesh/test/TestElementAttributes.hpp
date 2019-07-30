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


#ifndef _TESTELEMENTATTRIBUTES_HPP_
#define _TESTELEMENTATTRIBUTES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "OutputFileHandler.hpp"
#include "ElementAttributes.hpp"
#include "PetscTools.hpp"
#include "MutableElement.hpp"
#include "Node.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestElementAttributes : public CxxTest::TestSuite
{
public:

    void TestAttributesContainer()
    {
        ElementAttributes<2,2> element_attributes;

        element_attributes.AddAttribute(1.23);
        TS_ASSERT_EQUALS(element_attributes.rGetAttributes().size(), 1u);
        TS_ASSERT_DELTA(element_attributes.rGetAttributes()[0], 1.23, 1e-10);

        element_attributes.AddAttribute(4.56);
        TS_ASSERT_EQUALS(element_attributes.rGetAttributes().size(), 2u);
        TS_ASSERT_DELTA(element_attributes.rGetAttributes()[0], 1.23, 1e-10);
        TS_ASSERT_DELTA(element_attributes.rGetAttributes()[1], 4.56, 1e-10);

        element_attributes.SetFirstAttribute(42.0);
        TS_ASSERT_DELTA(element_attributes.GetFirstAttribute(), 42.0, 1e-10);
        TS_ASSERT_DELTA(element_attributes.rGetAttributes()[0], 42.0, 1e-10);

        // Now make a fresh attribute class with a size zero vector in it
        ElementAttributes<2,2> element_attributes2;
        // TS_ASSERT_THROWS_THIS(element_attributes2.GetFirstAttribute(), "Attempting to get element attribute when there are none defined");
        TS_ASSERT_DELTA(element_attributes2.GetFirstAttribute(), 0.0, 1e-10); ///\todo #2739 Currently defaults to zero
        element_attributes2.SetFirstAttribute(42.0);
        TS_ASSERT_DELTA(element_attributes2.GetFirstAttribute(), 42.0, 1e-10);
        TS_ASSERT_DELTA(element_attributes2.rGetAttributes()[0], 42.0, 1e-10);
    }

    void TestArchiving()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestElementAttributes", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "element_attributes.arch";

        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ElementAttributes<2,2>* p_element_attributes = new ElementAttributes<2,2>();

            p_element_attributes->AddAttribute(2.34);
            p_element_attributes->AddAttribute(5.67);

            ElementAttributes<2,2>* const p_const_element_attributes = p_element_attributes;

            output_arch << p_const_element_attributes;

            delete p_element_attributes;
        }

        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            ElementAttributes<2,2>* p_element_attributes;
            input_arch >> p_element_attributes;

            TS_ASSERT_EQUALS(p_element_attributes->rGetAttributes().size(), 2u);
            TS_ASSERT_DELTA(p_element_attributes->rGetAttributes()[0], 2.34, 1e-10);
            TS_ASSERT_DELTA(p_element_attributes->rGetAttributes()[1], 5.67, 1e-10);

            delete p_element_attributes;
        }
    }

    void TestArchivingNullPointer()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestElementAttributes", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "null_element_attributes.arch";

        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ElementAttributes<2,2>* const p_element_attributes = NULL;

            output_arch << p_element_attributes;

            delete p_element_attributes;
        }

        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            ElementAttributes<2,2>* p_element_attributes;
            input_arch >> p_element_attributes;

            TS_ASSERT(!p_element_attributes);

            delete p_element_attributes;
        }
    }

    void TestElementWithAttributes()
    {
        // Create an element
        MutableElement<2,2> this_element(0);

        // Check no attributes and exception thrown if accessed
        TS_ASSERT_EQUALS(this_element.GetNumElementAttributes(), 0u);
        TS_ASSERT_THROWS_THIS(this_element.rGetElementAttributes(), "Element has no attributes associated with it. Construct attributes first");

        // Add an attribute
        this_element.AddElementAttribute(1.23);

        // Check correct number and value
        TS_ASSERT_EQUALS(this_element.GetNumElementAttributes(), 1u);
        TS_ASSERT_DELTA(this_element.rGetElementAttributes()[0], 1.23, 1e-10);

        // Add a second attribute
        this_element.AddElementAttribute(4.56);

        // Check correct number and values
        TS_ASSERT_EQUALS(this_element.GetNumElementAttributes(), 2u);
        TS_ASSERT_DELTA(this_element.rGetElementAttributes()[0], 1.23, 1e-10);
        TS_ASSERT_DELTA(this_element.rGetElementAttributes()[1], 4.56, 1e-10);
    }
};

#endif //_TESTELEMENTATTRIBUTES_HPP_
