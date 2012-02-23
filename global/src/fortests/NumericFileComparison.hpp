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
#ifndef NUMERICFILECOMPARISON_HPP_
#define NUMERICFILECOMPARISON_HPP_

#include <cfloat>
#include <string>

#include "OutputFileHandler.hpp"
#include "MathsCustomFunctions.hpp"

#define A_WORD DBL_MAX
#define NOTHING_TO_READ DBL_MIN

/**
 * Compare files of numbers to see if they match to within a given tolerance.
 */
class NumericFileComparison
{
private:

    std::string mFilename1; /**< First filename */
    std::string mFilename2; /**< Second filename */

    std::ifstream* mpFile1; /**< First file */
    std::ifstream* mpFile2; /**< Second file */

public:

    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     *
     * @param fileName1  first file
     * @param fileName2  second file
     */
    NumericFileComparison(std::string fileName1, std::string fileName2):
        mFilename1(fileName1),
        mFilename2(fileName2)
    {
        mpFile1 = new std::ifstream(fileName1.c_str());

        // If it doesn't exist - throw exception
        if (!mpFile1->is_open())
        {
            delete mpFile1;
            mpFile1 = NULL;
            EXCEPTION("Couldn't open file: " + fileName1);
        }

        mpFile2 = new std::ifstream(fileName2.c_str());

        // If it doesn't exist - throw exception
        if (!mpFile2->is_open())
        {
            mpFile1->close();
            delete mpFile1;
            mpFile1 = NULL;
            delete mpFile2;
            mpFile2 = NULL;
            EXCEPTION("Couldn't open file: " + fileName2);
        }
    }

    /**
     * Close the files being compared.
     */
    ~NumericFileComparison()
    {
        if (mpFile1)
        {
            mpFile1->close();
        }
        if (mpFile2)
        {
            mpFile2->close();
        }
        delete mpFile1;
        delete mpFile2;
    }

    /**
     * Compare the files under both relative and absolute tolerances.
     * The comparison only fails if neither tolerance holds.  The
     * default settings effectively require numbers to match exactly.
     *
     * @param absTol  absolute tolerance on difference between numbers
     * @param ignoreFirstFewLines  how many lines to ignore from the comparison
     * @param relTol  relative tolerance on difference between numbers
     */
    bool CompareFiles(double absTol=DBL_EPSILON, unsigned ignoreFirstFewLines=0,
                      double relTol=DBL_EPSILON)
    {
        double data1;
        double data2;
        unsigned failures = 0;
        unsigned max_display_failures = 10;

        for (unsigned line_number=0; line_number<ignoreFirstFewLines; line_number++)
        {
            char buffer[1024];
            mpFile1->getline(buffer, 1024);
            mpFile2->getline(buffer, 1024);
            TS_ASSERT(!mpFile1->fail()); // Here we assume there are at least "ignoreFirstFewLines" lines...
            TS_ASSERT(!mpFile2->fail()); // ...and that they are lines of no more than 1024 characters
        }

        do
        {
            if (!(*mpFile1>>data1))
            {
                // Cannot read the next token from file as a number, so try a word instead
                std::string word;
                mpFile1->clear(); // reset the "failbit"
                if (*mpFile1 >> word)
                {
                    data1 = A_WORD;
                    if (word == "#" || word == "!")
                    {
                        // Ignore comment (up to 1024 characters until newline)
                        mpFile1->ignore(1024, '\n');
                    }
                }
                else
                {
                    mpFile1->clear(); // reset the "failbit"
                    data1 = NOTHING_TO_READ;
                }
            }
            if (!(*mpFile2 >> data2))
            {
                // Cannot read the next token from file as a number, so try a word instead
                std::string word;
                mpFile2->clear(); // reset the "failbit"
                if (*mpFile2 >> word)
                {
                    data2 = A_WORD;
                    if (word == "#" || word == "!")
                    {
                        // Ignore comment (up to 1024 characters until newline)
                        mpFile2->ignore(1024, '\n');
                    }
                }
                else
                {
                    mpFile2->clear(); // reset the "failbit"
                    data2 = NOTHING_TO_READ;
                }
            }

            bool ok = CompareDoubles::WithinAnyTolerance(data1, data2, relTol, absTol);
            if (!ok)
            {
                if (failures++ < max_display_failures)
                {
                    // Display error
                    CompareDoubles::WithinAnyTolerance(data1, data2, relTol, absTol, true);
                }
            }
        }
        while (data1 != NOTHING_TO_READ && data2 != NOTHING_TO_READ); // If either is a NOTHING_TO_READ, then it means that there's nothing to read from the file

        // Force CxxTest error if there were any major differences
        TS_ASSERT_EQUALS(failures, 0u);
        // If that assertion tripped...
        if (failures > 0u)
        {
#define COVERAGE_IGNORE
            // Report the paths to the files
            TS_TRACE("Files " + mFilename1 + " and " + mFilename2 + " numerically differ.");
#undef COVERAGE_IGNORE
        }
        return (failures==0);
    }
};

#endif /*NUMERICFILECOMPARISON_HPP_*/
