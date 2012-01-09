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

#ifndef TESTPROGRESSREPORTER_HPP_
#define TESTPROGRESSREPORTER_HPP_

#include <cxxtest/TestSuite.h>
#include "ProgressReporter.hpp"

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
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "progress_status.txt  global/test/data/good_progress_status.txt").c_str()), 0);
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
            catch(Exception e)
            {
            }
        }

        // File should now be closed since  progress_bar is out of scope
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "ProgressReporterException/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "progress_status.txt  global/test/data/bad_progress_status.txt").c_str()), 0);
    }
};

#endif /*TESTPROGRESSREPORTER_HPP_*/
