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

#ifndef TESTFILECOMPARISON_HPP_
#define TESTFILECOMPARISON_HPP_

#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"

class TestFileComparison : public CxxTest::TestSuite
{
public:

    void TestBasicFunctionality() throw(Exception)
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string noised_file = "./global/test/data/same_random_data_with_1e-4_noise.txt";

        // Comparing identical files shows no difference
        FileComparison same_data(base_file, base_file);
        TS_ASSERT(same_data.CompareFiles());

        // Comparing two different files gives a failure.
        FileComparison different_data(base_file, noised_file);
        TS_ASSERT_EQUALS(different_data.CompareFiles(0,false), false);
    }

    void TestIgnoreHeader() throw(Exception)
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string changed_file = "./global/test/data/same_random_data_with_different_header.txt";

        FileComparison file_comparer(base_file, changed_file);

        // Here we examine the header lines, find a difference and get a failure
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(0,false), false);

        // Here we ignore the header line and test passes.
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(1), true);
    }
};

#endif /*TESTFILECOMPARISON_HPP_*/
