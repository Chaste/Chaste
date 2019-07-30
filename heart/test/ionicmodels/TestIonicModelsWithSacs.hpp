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

#ifndef TESTIONICMODELSWITHSACS_HPP_
#define TESTIONICMODELSWITHSACS_HPP_


#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "HeartConfig.hpp"
#include "SimpleStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "OutputFileHandler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "TimeStepper.hpp"
#include "NumericFileComparison.hpp"

#include "FakePetscSetup.hpp"

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

    void TestN98WithSacAt300msShort()
    {
        double stretch = 1.1;
        double stretch_start_time = 300.0;
        double stretch_end_time = 305.0;

        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
                                   "TestN98WithSac", "sac300short", true);
    }

    void TestN98WithSacAt300msLong()
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

    void TestN98WithSacAt200msShort()
    {
        double stretch = 1.1;
        double stretch_start_time = 200.0;
        double stretch_end_time = 205.0;

        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
                                   "TestN98WithSac", "sac200short", false);
    }

    void TestN98WithSacAt800msLong()
    {

        double stretch = 1.1;
        double stretch_start_time = 800.0;
        double stretch_end_time = 810.0;

        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
                                   "TestN98WithSac", "sac800long", false);
    }

//// other possibilities - all fine - just uncomment and run
//    void TestN98WithSacAt200msLong()
//    {
//        double stretch = 1.1;
//        double stretch_start_time = 200.0;
//        double stretch_end_time = 210.0;
//
//        RunModelWithSacRecruitment(stretch, stretch_start_time, stretch_end_time,
//                                   "TestN98WithSac", "sac200long", false);
//    }
//
//    void TestN98WithSacAt800msShort()
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
//    void TestN98WithSacAt800msLong()
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
//    void TestN98WithSacAt230msLong()
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
