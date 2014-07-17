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

#ifndef TESTOUTPUTMODIFIERS_HPP_
#define TESTOUTPUTMODIFIERS_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "ArchiveLocationInfo.hpp"

#include "SingleTraceOutputModifier.hpp"
#include "ActivationOutputModifier.hpp"

class TestOutputModifiers : public CxxTest::TestSuite
{
public:

    void DontTestArchivingOfSingleTraceOutputModifier() throw(Exception)
    {
        OutputFileHandler handler("TestArchivingOfSingleTraceOutputModifier", false);
        // The next two lines ensure that different processes read/write different archive files when running in parallel
        ArchiveLocationInfo::SetArchiveDirectory(handler.FindFile(""));
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("SingleTraceOutputModifier.arch");

        // Create data structures to store variables to test for equality here

        // Save
        {
            AbstractOutputModifier* const p_abstract_class = new SingleTraceOutputModifier("SomeFileName", 123);

            // Create an output file
            std::ofstream ofs(archive_filename.c_str());
            // And create a boost output archive that goes to this file
            boost::archive::text_oarchive output_arch(ofs);

            // Record values to test into data structures
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class)
            // (if your abstract class doesn't contain the necessary variables and methods)

            output_arch << p_abstract_class;
            delete p_abstract_class;
        }

        // Load
        {
            AbstractOutputModifier* p_abstract_class_2;

            // Read from this input file
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            // And choose a boost input_archive object to translate this file
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_abstract_class_2;

            // Check things in the data structures with TS_ASSERTS here.
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class_2)
            // (if your abstract class doesn't contain the necessary variables and methods)

            TS_ASSERT_EQUALS(p_abstract_class_2->mFilename, "SomeFileName");

            delete p_abstract_class_2;

        }
    }
};

#endif // TESTOUTPUTMODIFIERS_HPP_
