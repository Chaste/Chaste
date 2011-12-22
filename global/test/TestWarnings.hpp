/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef _TESTWARNINGS_HPP_
#define _TESTWARNINGS_HPP_

#include <cxxtest/TestSuite.h>
#include "Warnings.hpp"
#include "LogFile.hpp"

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

        WARNING("This one goes into a log file");

        LogFile::Close();
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log_warnings.txt  global/test/data/log_warnings.txt").c_str()), 0);
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
