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

#ifndef TESTPROGRESSREPORTER_HPP_
#define TESTPROGRESSREPORTER_HPP_

#include <cxxtest/TestSuite.h>
#include "ProgressReporter.hpp"
#include "FileComparison.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestProgressReporter : public CxxTest::TestSuite
{
public:

    void TestBar()
    {
        // This is in a block to make sure the file gets closed (i.e. destructor called)
        {
            ProgressReporter progress_bar("ProgressReporter", 1.0, 10.0);

            progress_bar.PrintInitialising();
            for (unsigned i=0; i<=1000; i++)
            {
                double t = 1.0 + ((i+0.0)/1000)*9.0;
                progress_bar.Update(t);
            }

            progress_bar.PrintFinalising();
        }

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "ProgressReporter/";
        FileComparison(results_dir + "progress_status.txt", "global/test/data/good_progress_status.txt").CompareFiles();

    }

    void TestBarCleansFilesUpAfterException()
    {
        double smidge = 1e-8;

        {
            ProgressReporter progress_bar("ProgressReporterException", 1.0, 10.0);
            try
            {
                progress_bar.PrintInitialising();
                for (unsigned i=0; i<=900; i++)
                {
                    double t = 1.0 + ((i+0.0)/1000)*9.0 - smidge;
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
};

#endif /*TESTPROGRESSREPORTER_HPP_*/
