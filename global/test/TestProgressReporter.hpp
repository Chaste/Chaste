/*

Copyright (c) 2005-2021, University of Oxford.
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

#ifndef TESTPROGRESSREPORTER_HPP_
#define TESTPROGRESSREPORTER_HPP_

#include <cxxtest/TestSuite.h>
#include "ProgressReporter.hpp"
#include "FileComparison.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#include "Debug.hpp"

class TestProgressReporter : public CxxTest::TestSuite
{
public:

    void xTestBar()
    {
        // This is in a block to make sure the file gets closed (i.e. destructor called)
        {
            double dt = 9.0/1000.0;

            ProgressReporter progress_bar("ProgressReporter", 1.0, 10.0, dt);


            progress_bar.PrintInitialising();
            for (unsigned i=0; i<=1000; i++)
            {
                double t = 1.0 + (double)i * dt;
                progress_bar.Update(t);
            }

            progress_bar.PrintFinalising();
        }



        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "ProgressReporter/";

        PRINT_VARIABLE(results_dir);

        FileComparison(results_dir + "progress_status.txt", "global/test/data/good_progress_status.txt").CompareFiles();

    }

    void xTestBarCleansFilesUpAfterException()
    {
        double smidge = 1e-8;
        double dt = 9.0/1000.0;

        {
            ProgressReporter progress_bar("ProgressReporterException", 1.0, 10.0, dt);
            try
            {
                progress_bar.PrintInitialising();
                for (unsigned i=0; i<=900; i++)
                {
                    double t = 1.0 + (double)i * dt - smidge;
                    progress_bar.Update(t);
                }
                EXCEPTION("Throw");
            }
            catch(Exception&)
            {
            }
        }

        // File should now be closed since  progress_bar is out of scope
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "ProgressReporterException/";
        FileComparison(results_dir + "progress_status.txt", "global/test/data/bad_progress_status.txt").CompareFiles();
    }

    void TestGetTimeString()
    {
        // Invoke the default constructor
        ProgressReporter progress_reporter{};

        // Check standard calculations are being done correctly
        TS_ASSERT_EQUALS("(00:00:00:00s)", progress_reporter.GetTimeString(0.0));  // no time
        TS_ASSERT_EQUALS("(00:00:00:01s)", progress_reporter.GetTimeString(1.0));  // one second
        TS_ASSERT_EQUALS("(00:00:01:00s)", progress_reporter.GetTimeString(60.0));  // one minute
        TS_ASSERT_EQUALS("(00:01:00:00s)", progress_reporter.GetTimeString(3600.0));  // one hour
        TS_ASSERT_EQUALS("(01:00:00:00s)", progress_reporter.GetTimeString(86400.0));  // one day
        TS_ASSERT_EQUALS("(54:16:50:55s)", progress_reporter.GetTimeString(4726254.501));  // random less a bit
        TS_ASSERT_EQUALS("(54:16:50:55s)", progress_reporter.GetTimeString(4726255.499));  // random plus a bit

        // Negative duration should throw
        TS_ASSERT_THROWS_THIS(progress_reporter.GetTimeString(-1.0), "Invalid time elapsed.")

        // Obscenely long duration should throw
        auto long_time = static_cast<double>(LONG_MAX) + 10000.0;
        TS_ASSERT_THROWS_THIS(progress_reporter.GetTimeString(long_time), "Long overflow: something has gone wrong.")
    }
};

#endif /*TESTPROGRESSREPORTER_HPP_*/
