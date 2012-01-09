/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTLOGFILE_HPP_
#define TESTLOGFILE_HPP_

#include <cxxtest/TestSuite.h>
#include "LogFile.hpp"
#include "Exception.hpp"
#include <fstream>

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
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log2.txt  global/test/data/good_log2.txt").c_str()), 0);
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
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log.txt  global/test/data/good_log.txt").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log3.txt  global/test/data/good_log.txt").c_str()), 0);

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
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log4.txt  global/test/data/good_log4.txt").c_str()), 0);
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
