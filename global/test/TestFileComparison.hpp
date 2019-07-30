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

#ifndef TESTFILECOMPARISON_HPP_
#define TESTFILECOMPARISON_HPP_

#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

/*
 * NB - If copying and pasting a file comparison out of here then you should not include
 * the last argument to the constructor which suppresses its output,
 * and is not very useful for seeing what the difference in files is,
 * but it makes the output of this test much nicer.
 */
class TestFileComparison : public CxxTest::TestSuite
{
private:
    bool CalledCollectively;
    bool SuppressOutput;
    bool expected_fail_result;
public:

    void TestBasicFunctionality()
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string noised_file = "./global/test/data/same_random_data_with_1e-4_noise.txt";

        CalledCollectively = true;
        SuppressOutput = true;
        expected_fail_result = !PetscTools::AmMaster();

        // Comparing identical files shows no difference
        // NB Last two arguments can generally take their defaults (true and false)
        FileComparison same_data(base_file, base_file, CalledCollectively, SuppressOutput);
        TS_ASSERT(same_data.CompareFiles());

        // Comparing two different files gives a failure.
        FileComparison different_data(base_file, noised_file, CalledCollectively, SuppressOutput);

        TS_ASSERT_EQUALS(different_data.CompareFiles(0,false), expected_fail_result);
    }

    // By default the class expects to be used collectively, so only the master does file operations,
    // and there's a barrier before opening the files.  This case is "tested" in the other
    // tests - it's working if they don't deadlock!  Here we test the ability to be used non-collectively.
    void TestParallelOperation()
    {
        CalledCollectively = false;
        if (PetscTools::AmMaster())
        {
            std::string base_file = "./global/test/data/random_data.txt";
            FileComparison same_data(base_file, base_file, CalledCollectively);
            TS_ASSERT(same_data.CompareFiles());
        }
        CalledCollectively = true;
    }

    void TestFileFinderInterface()
    {
        FileFinder base_file("global/test/data/random_data.txt",
                             RelativeTo::ChasteSourceRoot);
        FileFinder noised_file("global/test/data/same_random_data_with_1e-4_noise.txt",
                               RelativeTo::ChasteSourceRoot);

        FileComparison same_data(base_file, noised_file, CalledCollectively, SuppressOutput);
        TS_ASSERT_EQUALS(same_data.CompareFiles(0,false), expected_fail_result);
    }

    void TestIgnoreHeader()
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string changed_file = "./global/test/data/same_random_data_with_different_header.txt";

        FileComparison file_comparer(base_file, changed_file, CalledCollectively, SuppressOutput);

        // Here we examine the header lines, find a difference and get a failure
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(0,false), expected_fail_result);

        // In a second file comparer make this pass by using SetIgnoreLinesBeginningWith();
        FileComparison file_comparer2(base_file, changed_file, CalledCollectively, SuppressOutput);
        file_comparer2.SetIgnoreLinesBeginningWith("header");
        TS_ASSERT_EQUALS(file_comparer2.CompareFiles(0,false), true);

        // Here we ignore the first line and test passes.
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(1), true);

        // Comment line differs so comparison fails
        file_comparer.SetIgnoreCommentLines(false);
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(1,false), expected_fail_result);
    }

    void TestIgnoreBlankLines()
    {
        std::string base_file = "./global/test/data/random_data.txt";
        std::string changed_file = "./global/test/data/same_random_data_with_blanks.txt";

        FileComparison file_comparer(base_file, changed_file, CalledCollectively, SuppressOutput);

        // A default compare fails because of the added blank lines
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(0, false), expected_fail_result);

        // Now make this pass using IgnoreBlankLines()
        file_comparer.IgnoreBlankLines();
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(0,false), true);

        // Now turn that off and make it fail again
        file_comparer.IgnoreBlankLines(false);
        TS_ASSERT_EQUALS(file_comparer.CompareFiles(0, false), expected_fail_result);
    }

    void TestBinaryFiles()
    {
        std::string base_file = "./mesh/test/data/simple_cube_binary.node";
        TS_ASSERT(FileComparison(base_file, base_file, CalledCollectively, SuppressOutput).CompareFiles());

        //A file which is the same for data purposes, but has the provenance line altered
          std::string copy_file = "./global/test/data/simple_cube_binary_copy.node";
        TS_ASSERT(FileComparison(base_file, copy_file, CalledCollectively, SuppressOutput).CompareFiles());

        //A file which has a single byte of data inserted
        std::string modified_file = "./global/test/data/simple_cube_binary_modified.node";

        TS_ASSERT_EQUALS(FileComparison(base_file, modified_file, CalledCollectively, SuppressOutput).CompareFiles(0, false), expected_fail_result);
    }
};

#endif /*TESTFILECOMPARISON_HPP_*/
