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
#ifndef FILECOMPARISON_HPP_
#define FILECOMPARISON_HPP_

#include "AbstractFileComparison.hpp"
#include <vector>
#include <boost/foreach.hpp>

/**
 * Compare files to check for any differences (in numeric and/or string values).
 *
 * By default this class ignores all lines which (in both files) start with '#' or '!'.
 */
class FileComparison : public AbstractFileComparison
{
private:
    /** Whether or not we should ignore comment lines. True by default. */
    bool mIgnoreCommentLines;

    /** Whether or not we should ignore blank lines. False by default. */
    bool mIgnoreBlankLines;

    /**
     * A list of strings, which if found at the beginning of lines when
     * #mIgnoreCommentLines is true, differences in these lines are ignored.
     */
    std::vector<std::string> mCommentLineStarts;

    /** Any lines which (in both files) contain one of these strings will be ignored. */
    std::vector<std::string> mIgnorableContent;

public:
    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     *
     * @param fileName1  first file
     * @param fileName2  second file
     * @param calledCollectively  if true there will be a barrier before opening files, and only master compares contents.
     * @param suppressOutput  if true then no errors will go to TS_TRACE(). Should only be set for the test of this class.
     */
    FileComparison(std::string fileName1, std::string fileName2, bool calledCollectively=true, bool suppressOutput = false)
        : AbstractFileComparison(fileName1, fileName2, calledCollectively, suppressOutput),
          mIgnoreCommentLines(true),
          mIgnoreBlankLines(false)
    {
        SetupCommentLines();
    }

    /**
     * Specify two files to compare, and open them for reading.
     * Actual comparison is done by calling CompareFiles.
     *
     * @param rFileName1  first file finder
     * @param rFileName2  second file finder
     * @param calledCollectively  if true there will be a barrier before opening files, and only master compares contents.
     * @param suppressOutput  if true then no errors will go to TS_TRACE(). Should only be set for the test of this class.
     */
    FileComparison(const FileFinder& rFileName1, const FileFinder& rFileName2, bool calledCollectively=true, bool suppressOutput = false)
        : AbstractFileComparison(rFileName1, rFileName2, calledCollectively, suppressOutput),
          mIgnoreCommentLines(true),
          mIgnoreBlankLines(false)
    {
        SetupCommentLines();
    }

    /**
     * Set some line starts that define comments in the files.
     *
     * These are ignored by default and when #mIgnoreCommentLines is explicitly set to true.
     */
    void SetupCommentLines()
    {
        mCommentLineStarts.push_back("#");
        mCommentLineStarts.push_back("!");
        mCommentLineStarts.push_back("Created by Chaste");
        mCommentLineStarts.push_back("<!-- Created by Chaste");
    }

    /**
     * Whether or not we should ignore lines starting with a comment symbol
     * (the default symbols are '#' and '!', set by SetupCommentLines).
     * @param ignore  whether to ignore these lines.  If set to false, existing
     *      defined comment symbols are also cleared, so that an entirely new
     *      set may be defined with SetIgnoreLinesBeginningWith.
     */
    void SetIgnoreCommentLines(bool ignore=true)
    {
        mIgnoreCommentLines = ignore;
        if (!ignore)
        {
            mCommentLineStarts.clear();
        }
    }

    /**
     * Set whether or not we should ignore blank lines.  If only one of the files being
     * compared has a blank line, if this mode is on we just skip that line and see if the
     * next content line matches the other file.
     * @param ignore  whether to ignore blank lines appearing only in one file
     */
    void IgnoreBlankLines(bool ignore=true)
    {
        mIgnoreBlankLines = ignore;
    }

    /**
     * Set an additional line start which should be treated as a comment and ignored
     * (and therefore switch on #mIgnoreCommentLines = true).
     *
     * @param lineStart  the beginning of a line which should be treated as a comment
     */
    void SetIgnoreLinesBeginningWith(std::string lineStart)
    {
        mIgnoreCommentLines = true;
        mCommentLineStarts.push_back(lineStart);
    }

    /**
     * Add the given string to #mIgnorableContent, and hence ignore differences in lines
     * which contain that text in both files.
     *
     * @param rIgnorableText  the text indicating lines to ignore
     */
    void IgnoreLinesContaining(const std::string& rIgnorableText)
    {
        mIgnorableContent.push_back(rIgnorableText);
    }

    /**
     * @return true if the files are identical, barring ignored content.
     *
     * @param ignoreFirstFewLines  how many lines to ignore from the comparison
     * @param doTsAssert  whether to throw a TS_ASSERT internally (switched off for testing only)
     */
    bool CompareFiles(unsigned ignoreFirstFewLines=0, bool doTsAssert=true)
    {
        // Usually only the master process does the checking, this can be switched off in the constructor.
        if (mCalledCollectively && !PetscTools::AmMaster())
        {
            return true;
        }

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
            getline(*mpFile1, buffer1);
            getline(*mpFile2, buffer2);

            if (mIgnoreBlankLines)
            {
                // Keep reading lines until we see non-blank, end-of-file or read error
                while (buffer1.empty() && mpFile1->good())
                {
                    getline(*mpFile1, buffer1);
                }
                while (buffer2.empty() && mpFile2->good())
                {
                    getline(*mpFile2, buffer2);
                }
            }

            if (mIgnoreCommentLines)
            {
                bool skip_this_line = false;
                for (unsigned i=0; i<mCommentLineStarts.size(); i++)
                {
                    // Check for lines starting with a comment symbol
                    size_t found1 = buffer1.find(mCommentLineStarts[i]);
                    size_t found2 = buffer2.find(mCommentLineStarts[i]);
                    if (found1 == 0 && found2 == 0)
                    {
                        skip_this_line = true;
                        break;
                    }
                }
                if (skip_this_line)
                {
                    continue;
                }
            }

            // Check for lines containing ignorable text
            if (!mIgnorableContent.empty())
            {
                bool skip_this_line = false;
                BOOST_FOREACH(const std::string& rText, mIgnorableContent)
                {
                    size_t found1 = buffer1.find(rText);
                    size_t found2 = buffer2.find(rText);
                    if (found1 != std::string::npos && found2 != std::string::npos)
                    {
                        skip_this_line = true;
                        break;
                    }
                }
                if (skip_this_line)
                {
                    continue;
                }
            }

            if (!(buffer1==buffer2) && !files_empty)
            {
                if (failures++ < max_display_failures && !mSuppressOutput)
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
            if (failures > 0u && !mSuppressOutput)
            {
                // Report the paths to the files
                TS_TRACE("Files " + mFilename1 + " and " + mFilename2 + " differ.");
            }
        }

        ResetFiles();

        return (failures==0);
    }
};

#endif /*FILECOMPARISON_HPP_*/
