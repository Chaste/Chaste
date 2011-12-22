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

#ifndef TESTIONICMODELSWITHSACS_HPP_
#define TESTIONICMODELSWITHSACS_HPP_


#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "HeartConfig.hpp"
#include "SimpleStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "OutputFileHandler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "TimeStepper.hpp"
#include "NumericFileComparison.hpp"

// This is a separate class to TestIonicModels as we don't want to be calling
// RunAndCheckIonicModels, because we want to switch the SAC on at certain times..

class TestIonicModelsWithSacs : public CxxTest::TestSuite
{
private:
    //  Run normally up to stretch-time, then apply stretch and
    //  run until stretch-off-time, then return to stretch=1
    //  and run until 1000ms.
    void RunModelWithSacRecruitment(double stretch,
                                    double stretchStartTime,
                                    double stretchEndTime,
                                    std::string directory,
                                    std::string filePrefix,
                                    bool clearDir)
    {
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(
                -3 /*magnitude*/,
                3 /*duration*/,
                10.0 /*start time*/));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        double time_step = 0.01;
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, 1.0);

        CML_noble_varghese_kohl_noble_1998_basic_with_sac   n98_with_sac(p_solver, p_stimulus);

        OutputFileHandler handler(directory,clearDir);
        out_stream p_file = handler.OpenOutputFile(filePrefix+".dat");
        *p_file << 0 << " " << n98_with_sac.GetVoltage() << "\n";

        double printing_dt = 1;
        TimeStepper stepper(0, stretchStartTime, printing_dt);

        while ( !stepper.IsTimeAtEnd() )
        {
            n98_with_sac.Compute(stepper.GetTime(), stepper.GetNextTime());
            stepper.AdvanceOneTimeStep();

            *p_file << stepper.GetTime() << " " << n98_with_sac.GetVoltage() << "\n";
        }

        n98_with_sac.SetStretch(stretch);

        TimeStepper stepper2(stepper.GetTime(), stretchEndTime, printing_dt);
        while ( !stepper2.IsTimeAtEnd() )
        {
            n98_with_sac.Compute(stepper2.GetTime(), stepper2.GetNextTime());
            stepper2.AdvanceOneTimeStep();

            *p_file << stepper2.GetTime() << " " << n98_with_sac.GetVoltage() << "\n";
        }

        n98_with_sac.SetStretch(1.0);

        TimeStepper stepper3(stepper2.GetTime(), 1000, printing_dt);
        while ( !stepper3.IsTimeAtEnd() )
        {
            n98_with_sac.Compute(stepper3.GetTime(), stepper3.GetNextTime());
            stepper3.AdvanceOneTimeStep();

            *p_file << stepper3.GetTime() << " " << n98_with_sac.GetVoltage() << "\n";
        }

        p_file->close();
    }

public:

    // series of tests. one has a hardcoded test. See figure in #1345.

    void TestN98WithSacAt300msShort() throw(Exception)
    {
        double stretch = 1.1;
        double stretch_start_time = 300.0;
        double stretch_end_time = 305.0;

        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
                                   "TestN98WithSac", "sac300short", true);
    }

    void TestN98WithSacAt300msLong() throw(Exception)
    {
        double stretch = 1.1;
        double stretch_start_time = 300.0;
        double stretch_end_time = 310.0;

        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
                                   "TestN98WithSac", "sac300long", false);

        OutputFileHandler handler("TestN98WithSac",false);
        std::string file1 = handler.GetOutputDirectoryFullPath() + "sac300long.dat";
        std::string file2 = "heart/test/data/N98WithSac_sac300mslong.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles());
    }

    void TestN98WithSacAt200msShort() throw(Exception)
    {
        double stretch = 1.1;
        double stretch_start_time = 200.0;
        double stretch_end_time = 205.0;

        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
                                   "TestN98WithSac", "sac200short", false);
    }

    void TestN98WithSacAt800msLong() throw(Exception)
    {

        double stretch = 1.1;
        double stretch_start_time = 800.0;
        double stretch_end_time = 810.0;

        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
                                   "TestN98WithSac", "sac800long", false);
    }

//// other possibilities - all fine - just uncomment and run
//    void TestN98WithSacAt200msLong() throw(Exception)
//    {
//        double stretch = 1.1;
//        double stretch_start_time = 200.0;
//        double stretch_end_time = 210.0;
//
//        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
//                                   "TestN98WithSac", "sac200long", false);
//    }
//
//    void TestN98WithSacAt800msShort() throw(Exception)
//    {
//
//        double stretch = 1.1;
//        double stretch_start_time = 800.0;
//        double stretch_end_time = 805.0;
//
//        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
//                                   "TestN98WithSac", "sac800short", false);
//    }
//
//    void TestN98WithSacAt800msLong() throw(Exception)
//    {
//
//        double stretch = 1.1;
//        double stretch_start_time = 800.0;
//        double stretch_end_time = 810.0;
//
//        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
//                                   "TestN98WithSac", "sac800long", false);
//    }
//
//    void TestN98WithSacAt230msLong() throw(Exception)
//    {
//
//        double stretch = 1.1;
//        double stretch_start_time = 230.0;
//        double stretch_end_time = 240.0;
//
//        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
//                                   "TestN98WithSac", "sac230long", false);
//    }

};
#endif /*TESTIONICMODELSWITHSACS_HPP_*/
