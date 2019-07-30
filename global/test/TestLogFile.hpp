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

#ifndef TESTLOGFILE_HPP_
#define TESTLOGFILE_HPP_

#include <fstream>
#include <cxxtest/TestSuite.h>
#include "LogFile.hpp"
#include "Exception.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestLogFile : public CxxTest::TestSuite
{
public:

    void TestLogFileCreate()
    {
        LogFile* p_log_file = LogFile::Instance();

        // no file set yet
        TS_ASSERT_EQUALS(p_log_file->IsFileSet(), false);

        // set the file
        p_log_file->Set(1, "TestLogFile");
        TS_ASSERT_EQUALS(p_log_file->IsFileSet(), true);

        // check a new instance works correctly
        LogFile* p_same_log = LogFile::Instance();
        TS_ASSERT_EQUALS(p_same_log->IsFileSet(), true);
    }

    void TestLogStillExists()
    {
        LogFile* p_log = LogFile::Instance();
        TS_ASSERT_EQUALS(p_log->IsFileSet(), true);
    }

    void TestClose()
    {
        LogFile::Close();

        // Check file not set on a new instance
        LogFile* p_log = LogFile::Instance();
        TS_ASSERT_EQUALS(p_log->IsFileSet(), false);
    }

    void TestWritingToFile1()
    {
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(1, "TestLogFile", "log2.txt");

        (*p_log_file) << "Some stuff\n" << "Some more\n";
        (*p_log_file) << "Even more\n";
    }

    void TestWritingToFile2()
    {
        LogFile* p_log_file = LogFile::Instance();

        (*p_log_file) << ".. and another bit\n";

        p_log_file->SetPrecision(9);
        (*p_log_file) << M_PI << "\n";

        (*LogFile::Instance()) << "..and one final bit\n";
        LogFile::Close();

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";
        FileComparison(results_dir + "log2.txt", "global/test/data/good_log2.txt").CompareFiles();
    }

    void TestWritingToNoFile()
    {
        LogFile::Close();

        // Test no seg faults etc
        (*LogFile::Instance()) << "this won't be written anywhere, as no log file has been created";
    }

    void TestWritingToNewFiles()
    {
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(1, "TestLogFile");
        (*p_log_file) << "data";

        // Open a new file without closing the previous
        p_log_file->Set(1, "TestLogFile","log3.txt");
        (*p_log_file) << "data";

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";
        FileComparison(results_dir + "log.txt",  "global/test/data/good_log.txt").CompareFiles();
        FileComparison(results_dir + "log3.txt", "global/test/data/good_log.txt").CompareFiles();

        LogFile::Close();
    }

    void TestUsingMacroAndLevels()
    {
        LogFile* p_log_file = LogFile::Instance();

        // Bad level
        TS_ASSERT_THROWS_THIS(p_log_file->Set( LogFile::MaxLoggingLevel()+1, "TestLogFile"),
                "Requested level 3 should have been less than or equal to 2");

        p_log_file->Set(1, "TestLogFile", "log4.txt");

        unsigned i = 0;

        LOG(1, "Level 1 info, will be written. i = " << i);
        LOG(2, "Level 2 info, WONT be written. i = " << i);

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";

        // This will fail if optimised (and should fail) since the NDEBUG flag currently forces NO LOGGING
        FileComparison(results_dir + "log4.txt", "global/test/data/good_log4.txt").CompareFiles();

        LogFile::Close();
    }

    void TestHeaderAndElapsedTime()
    {
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(1, "TestLogFile", "log5.txt");

        p_log_file->WriteHeader("Complete human");

        // for (unsigned i=0; i<1e9; i++);
        p_log_file->WriteElapsedTime(" -> ");

        LogFile::Close();

        // The file will be different on different occasions (as the date is printed), so test by reading it in
        std::ifstream ifs((OutputFileHandler::GetChasteTestOutputDirectory()+"TestLogFile/log5.txt").c_str());
        if (ifs.is_open())
        {
            std::string line;

            // Get the second line
            getline(ifs,line);
            getline(ifs,line);

            // The date will change but the beginning of the line won't
            std::string expected_beginning_of_line = "Chaste: Complete human simulation, on";
            TS_ASSERT_EQUALS(line.substr(0,expected_beginning_of_line.size()),
                             expected_beginning_of_line);

            // Get the fourth line
            getline(ifs,line);
            getline(ifs,line);

            //Hopefully it took less than one minute(!) to do a tiny bit of writing..
            expected_beginning_of_line = " -> Elapsed time is: 0h 0m";
            TS_ASSERT_EQUALS(line.substr(0,expected_beginning_of_line.size()),
                             expected_beginning_of_line);
        }
        else
        {
            TS_FAIL("log file not written?");
        }
    }
};

#endif /*TESTLOGFILE_HPP_*/
