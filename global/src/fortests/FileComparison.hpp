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
#ifndef FILECOMPARISON_HPP_
#define FILECOMPARISON_HPP_

#include "AbstractFileComparison.hpp"

/**
 * Compare files to check for any differences (in numeric and/or string values).
 *
 * By Default this class ignores all lines which (in both files) start with '#'.
 */
class FileComparison : public AbstractFileComparison
{
private:
	/** Whether or not we should ignore lines starting with '#'. True by default.*/
    bool mIgnoreCommentLines;

    /**
     * A list of strings, if found at the beginning of lines when
     * mIgnoreCommentLines is true, differences in these lines
     * are ignored.
     */
    std::vector<std::string> mCommentLineStarts;
public:

    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     *
     * @param fileName1  first file
     * @param fileName2  second file
     */
    FileComparison(std::string fileName1, std::string fileName2):
       AbstractFileComparison(fileName1,fileName2),
       mIgnoreCommentLines(true)
    {
        SetupCommentLines();
    }

    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     *
     * @param rFileName1  first file finder
     * @param rFileName2  second file finder
     */
    FileComparison(const FileFinder& rFileName1, const FileFinder& rFileName2):
       AbstractFileComparison(rFileName1,rFileName2),
       mIgnoreCommentLines(true)
    {
        SetupCommentLines();
    }

    /**
     * Set some line starts that define comments in the files.
     *
     * These are ignored by default and when mIgnoreCommentLines = true.
     */
    void SetupCommentLines()
    {
        mCommentLineStarts.push_back("#");
        mCommentLineStarts.push_back("!");
        mCommentLineStarts.push_back("Created by Chaste");
    }

    /**
     * Whether or not we should ignore lines starting with '#'
     * @param ignore  whether to ignore these lines (defaults to true)
     */
    void SetIgnoreCommentLines(bool ignore=true)
    {
        mIgnoreCommentLines = ignore;
    }

    /**
     * Set a line start which should be treated as a comment and ignored
     * (and therefore switch on mIgnoreCommentLines = true)
     *
     * @param rLineStart  The beginning of a line which should be treated as a comment
     */
    void SetIgnoreLinesBeginningWith(std::string lineStart)
    {
        mIgnoreCommentLines = true;
        mCommentLineStarts.push_back(lineStart);
    }

    /**
     * Compare the files under both relative and absolute tolerances.
     * The comparison only fails if neither tolerance holds.  The
     * default settings effectively require numbers to match exactly.
     *
     * @param ignoreFirstFewLines  how many lines to ignore from the comparison
     * @param doTsAssert  Whether to throw a TS_ASSERT internally (switched off for testing only)
     */
    bool CompareFiles(unsigned ignoreFirstFewLines=0, bool doTsAssert=true)
    {
        std::string data1;
        std::string data2;
        unsigned failures = 0;
        unsigned max_display_failures = 10;

        SkipHeaderLines(ignoreFirstFewLines);

        bool files_empty = false;
        do
        {
            std::string buffer1;
            std::string buffer2;
            getline(*mpFile1,buffer1);
            getline(*mpFile2,buffer2);

            if (mIgnoreCommentLines)
            {
                bool skip_this_line = false;
                for (unsigned i=0; i<mCommentLineStarts.size(); i++)
                {
                    // Check for lines starting with "#"
                    size_t found1 = buffer1.find(mCommentLineStarts[i]);
                    size_t found2 = buffer2.find(mCommentLineStarts[i]);
                    if (found1 == 0 && found2 == 0)
                    {
                        skip_this_line = true;
                    }
                }
                if (skip_this_line)
                {
                    continue;
                }
            }

            if (!(buffer1==buffer2) && !files_empty)
            {
                if (failures++ < max_display_failures)
                {
                    // Display error
                    std::stringstream message;
                    message << "Line " << mLineNum << " differs in files " << mFilename1 << " and " << mFilename2;
                    TS_TRACE(message.str());
                    TS_TRACE( buffer1 );
                    TS_TRACE( buffer2 );
                }
            }
            mLineNum++;
        }
        while (mpFile1->good() && mpFile2->good());
        // If either is not good(), then it means that there's nothing to read from the file, or a file input error.

        if (doTsAssert)
        {
            // Force CxxTest error if there were any major differences
            TS_ASSERT_EQUALS(failures, 0u);
            // If that assertion tripped...
            if (failures > 0u)
            {
#define COVERAGE_IGNORE
                // Report the paths to the files
                TS_TRACE("Files " + mFilename1 + " and " + mFilename2 + " differ.");
#undef COVERAGE_IGNORE
            }
        }

        ResetFiles();

        return (failures==0);
    }
};

#endif /*FILECOMPARISON_HPP_*/
