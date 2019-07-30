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

#ifndef ABSTRACTFILECOMPARISON_HPP_
#define ABSTRACTFILECOMPARISON_HPP_

#include <string>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"

/**
 * Abstract class for comparing two files, looking for differences in tests.
 */
class AbstractFileComparison
{
public:

    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     *
     * @param rFileFinder1  first file
     * @param rFileFinder2  second file
     * @param calledCollectively  If true there will be a barrier before opening files, and only master compares contents.
     * @param suppressOutput  If true then no errors will go to TS_TRACE(). Should only be set for the test of this class.
     */
    AbstractFileComparison(const FileFinder& rFileFinder1,
                           const FileFinder& rFileFinder2,
                           bool calledCollectively,
                           bool suppressOutput):
        mFilename1(rFileFinder1.GetAbsolutePath()),
        mFilename2(rFileFinder2.GetAbsolutePath()),
        mCalledCollectively(calledCollectively),
        mSuppressOutput(suppressOutput)
    {
        Setup();
    }

    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     *
     * @param fileName1  first file
     * @param fileName2  second file
     * @param calledCollectively  If true there will be a barrier before opening files, and only master compares contents.
     * @param suppressOutput  If true then no errors will go to TS_TRACE(). Should only be set for the test of this class.
     */
    AbstractFileComparison(std::string fileName1,
                           std::string fileName2,
                           bool calledCollectively,
                           bool suppressOutput):
        mFilename1(fileName1),
        mFilename2(fileName2),
        mCalledCollectively(calledCollectively),
        mSuppressOutput(suppressOutput)
    {
        Setup();
    }

    /**
     * Close the files being compared.
     */
    virtual ~AbstractFileComparison()
    {
        if (mpFile1)
        {
            mpFile1->close();
            delete mpFile1;
        }
        if (mpFile2)
        {
            mpFile2->close();
            delete mpFile2;
        }
    }

protected:

    std::string mFilename1; /**< First filename */
    std::string mFilename2; /**< Second filename */

    std::ifstream* mpFile1; /**< First file */
    std::ifstream* mpFile2; /**< Second file */

    unsigned mLineNum; /**< Counter for the line number we are on (in FileComparision) */

    bool mCalledCollectively; /**< If true there will be a barrier before opening files, and only master compares contents. */

    /** Whether we should suppress output from this class, just for a clean looking test */
    bool mSuppressOutput;

    /**
     * This method closes and reopens files so that another CompareFiles() command can be run on the same object.
     */
    void ResetFiles()
    {
        if (!mCalledCollectively || PetscTools::AmMaster())
        {
            // We want to reset the files to allow this method to be called again, with different tolerances for instance.
            mpFile1->close();
            mpFile2->close();
            mpFile1->open(mFilename1.c_str());
            mpFile2->open(mFilename2.c_str());
            mLineNum = 1u;
        }
    }

    /**
     * This helper method just moves forward the two file pointers to skip some header lines.
     * @param numLinesToSkip  The number of header lines to skip
     */
    void SkipHeaderLines(unsigned numLinesToSkip)
    {
        if (!mCalledCollectively || PetscTools::AmMaster())
        {
            for (unsigned line_number=0; line_number<numLinesToSkip; line_number++)
            {
                char buffer[1024];
                mpFile1->getline(buffer, 1024);
                mpFile2->getline(buffer, 1024);
                TS_ASSERT(!mpFile1->fail()); // Here we assume there are at least "ignoreFirstFewLines" lines...
                TS_ASSERT(!mpFile2->fail()); // ...and that they are lines of no more than 1024 characters
                mLineNum++;
            }
        }
    }

private:
    /**
     * Private method called only by the two constructors.
     */
    void Setup()
    {
        if (mCalledCollectively)
        {
            PetscTools::Barrier("AbstractFileComparison::Setup");
        }
        if (!mCalledCollectively || PetscTools::AmMaster())
        {
            mpFile1 = new std::ifstream(mFilename1.c_str());

            // If it doesn't exist - throw exception
            if (!mpFile1->is_open())
            {
                delete mpFile1;
                mpFile1 = NULL;
                EXCEPTION("Couldn't open file: " + mFilename1);
            }

            mpFile2 = new std::ifstream(mFilename2.c_str());

            // If it doesn't exist - throw exception
            if (!mpFile2->is_open())
            {
                mpFile1->close();
                delete mpFile1;
                mpFile1 = NULL;
                delete mpFile2;
                mpFile2 = NULL;
                EXCEPTION("Couldn't open file: " + mFilename2);
            }

            mLineNum = 1u;
        }
        else
        {
            mpFile1 = NULL;
            mpFile2 = NULL;
        }
    }
};

#endif // ABSTRACTFILECOMPARISON_HPP_
