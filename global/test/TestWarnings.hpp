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

#ifndef _TESTWARNINGS_HPP_
#define _TESTWARNINGS_HPP_

#include <cxxtest/TestSuite.h>
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "FileComparison.hpp"

class TestWarnings: public CxxTest::TestSuite
{
private:
    void ThrowWarning()
    {
        WARN_ONCE_ONLY("Ozzy is Hungry.");
    }

public:

    void TestGetMessage()
    {
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        WARNING("Ozzy is near.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        WARNING("Stay alert.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Ozzy is near.");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::QuietDestroy();
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
    }

    void TestWarningsWithLogging()
    {
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(1, "TestLogFile", "log_warnings.txt");

        WARNING("This one goes into a log file"); //The line number for this WARNING is important since it gets written to the log file

        LogFile::Close();
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";
        FileComparison(results_dir + "log_warnings.txt", "global/test/data/log_warnings.txt").CompareFiles();
        Warnings::QuietDestroy();
    }

    void TestWarningOnceOnly()
    {
        for (unsigned year=1970; year<2100; year+=4)
        {
            WARN_ONCE_ONLY("Don't get your hopes up, England are not going to win the World Cup.");
        }
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        //And a streaming warning
        std::string best_hope = "Murray";
        for (unsigned year=2005; year<2020; year++)
        {
            WARN_ONCE_ONLY("Don't get your hopes up, " << best_hope <<" is not winning Wimbledon.");
        }
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        Warnings::QuietDestroy();
    }

    void TestPrintWarnings()
    {

    	WARNING("This warning is printed inside TestPrintWarnings().");
    	TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
    	Warnings::PrintWarnings();
    	TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
    	Warnings::QuietDestroy();
    }

    void TestWarningOnlyOnceReset()
    {
        for (unsigned i=0; i<10; i++)
        {
            ThrowWarning();
        }
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::QuietDestroy();

        ThrowWarning();
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::QuietDestroy();
    }

    void TestLastTestWithWarningsIsNoisy() //Needs to happen last (after any QuietDestroy()), so that a warning is printed
    {
        TS_ASSERT_THROWS_THIS(Warnings::Instance()->GetNextWarningMessage(),"There are no warnings");
        unsigned one = 1;
        WARNING("This one will get printed " << one << " time");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
    }

};

#endif //_TESTWARNINGS_HPP_
