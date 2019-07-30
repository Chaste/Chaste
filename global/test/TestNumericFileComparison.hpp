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

#ifndef TESTNUMERICFILECOMPARISON_HPP_
#define TESTNUMERICFILECOMPARISON_HPP_

#include <cxxtest/TestSuite.h>
#include "NumericFileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

/*
 * NB - If copying and pasting a file comparison out of here then you should not include
 * the last argument to the constructor which suppresses its output,
 * as it is not very useful for seeing what the difference in files is.
 * But it makes the output of this test much nicer.
 */
class TestNumericFileComparison : public CxxTest::TestSuite
{
private:
    bool CalledCollectively;
    bool SuppressOutput;
    bool expected_fail_result;
public:

    void TestBasicFunctionality()
    {
        CalledCollectively = true;
        SuppressOutput = true;
        expected_fail_result = !PetscTools::AmMaster();

        std::string base_file = "./global/test/data/random_data.txt";
        std::string noised_file = "./global/test/data/same_random_data_with_1e-4_noise.txt";

        NumericFileComparison same_data(base_file, base_file, CalledCollectively, SuppressOutput);
        TS_ASSERT(same_data.CompareFiles());

        NumericFileComparison different_data(base_file, noised_file, CalledCollectively, SuppressOutput);
        TS_ASSERT(different_data.CompareFiles(1e-4));

        TS_ASSERT_EQUALS(different_data.CompareFiles(1e-9, 0, 1e-9, false), expected_fail_result);
    }

    void TestIgnoreHeader()
    {
        std::string boost_33_file = "./global/test/data/fake_archive_boost_1_33.txt";
        std::string boost_34_file = "./global/test/data/fake_archive_boost_1_34.txt";

        NumericFileComparison same_data(boost_33_file, boost_34_file, CalledCollectively, SuppressOutput);
        TS_ASSERT(same_data.CompareFiles(5.1e-6, 1));

        TS_ASSERT_EQUALS(same_data.CompareFiles(1e-9, 1, 1e-9, false), expected_fail_result);
    }

    void TestIgnoreProvenanceComment()
    {
        std::string v1_file = "./global/test/data/fake_v_day_one.txt";
        std::string v2_file = "./global/test/data/fake_v_day_two.txt";

        NumericFileComparison same_data(v1_file, v2_file, CalledCollectively, SuppressOutput);

        TS_ASSERT_EQUALS(same_data.CompareFiles(5e-4, 1, 1e-9, false), expected_fail_result); // Fails due to difference after comment
        TS_ASSERT(same_data.CompareFiles(5e-3)); // Difference after comment is below this tolerance
    }

    void TestRelativeDifference()
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string noised_file = "./global/test/data/same_random_data_with_1e-4_noise.txt";
        NumericFileComparison different_data(base_file, noised_file, CalledCollectively, SuppressOutput);

        // Lower bound on data is 1e-2 so 1e-4 absolute noise is within 1e-2 relative tolerance

        // Here we would fail absolute and pass relative -> pass.
        TS_ASSERT(different_data.CompareFiles(1e-5, 0, 1e-2, false));

        // Here we should pass absolute and fail relative -> pass
        TS_ASSERT(different_data.CompareFiles(1e-3, 0, 1e-9, false));

        // Here we should fail both -> fail
        TS_ASSERT_EQUALS(different_data.CompareFiles(1e-5, 0, 1e-3, false), expected_fail_result);

    }
};

#endif /*TESTNUMERICFILECOMPARISON_HPP_*/
