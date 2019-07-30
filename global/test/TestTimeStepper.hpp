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

#ifndef TESTTIMESTEPPER_HPP_
#define TESTTIMESTEPPER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cassert>
#include <cmath>
#include <cfloat>

#include "TimeStepper.hpp"
#include "OutputFileHandler.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"


class TestTimeStepper : public CxxTest::TestSuite
{
public:

    void TestOverflow()
    {
        TimeStepper stepper(0.0, DBL_MAX, DBL_EPSILON);
        stepper.mTotalTimeStepsTaken = (unsigned)(-1);
        TS_ASSERT(!stepper.IsTimeAtEnd());
        TS_ASSERT_THROWS_THIS(stepper.AdvanceOneTimeStep(), "Time step counter has overflowed.");
    }

    void TestAdvance()
    {
        const double smidge = 1e-10;

        double start_time = 0.0;
        double end_time = 2.0;
        double timestep = 3.7e-05;

        // This is how a time stepper is normally used
        TimeStepper my_stepper(start_time, end_time, timestep);
        while ( !my_stepper.IsTimeAtEnd() )
        {
            // do something

            my_stepper.AdvanceOneTimeStep();
        }

        // Tests
        TS_ASSERT_THROWS_THIS(TimeStepper(end_time, start_time, timestep),
                              "The simulation duration must be positive, not -2");

        // Note that the timestep in this test does not divide whole time interval nicely, therefore we can't enforce a constant timestep.
        bool enforce_constant_timestep = true;
        TS_ASSERT_THROWS_CONTAINS(TimeStepper(start_time, end_time, timestep, enforce_constant_timestep),
                                      "TimeStepper estimates non-constant timesteps will need to be used");

        TimeStepper stepper(start_time, end_time, timestep);

        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), (unsigned) floor((end_time - start_time)/timestep) );

        double real_time_step = timestep;
        unsigned time_step_number = 0;
        double current_time = start_time;

        /*
         * We'll trap for stopping times that are close to the end time
         * in order to avoid having a timestep of 1e-14 (or whatever) at
         * the end in the case of rounding errors.
         */
        double close_to_end_time = end_time - smidge*timestep;

        while (current_time < end_time)
        {
            TS_ASSERT(!stepper.IsTimeAtEnd());

            time_step_number++;

            // Determine what the value time step should really be like
            double to_time = start_time+time_step_number*timestep;

            if (to_time >= close_to_end_time)
            {
                real_time_step = end_time - current_time;
                //std::cout<<"TImes: " << timestep << " " << stepper.GetNextTimeStep() <<", difference = " << timestep - stepper.GetNextTimeStep()  << "\n";
                to_time = end_time;

                // Note that with non-constant timesteps, the final timestep can vary a lot (but not more than a normal sized timestep!).
                TS_ASSERT_DELTA(stepper.GetNextTimeStep(), timestep, timestep);
            }
            else
            {
                // Otherwise the GetNextTimeStep() returns the stored timestep,
                TS_ASSERT_EQUALS(stepper.GetNextTimeStep(),  timestep);

            }
            TS_ASSERT_DELTA(stepper.GetNextTimeStep(), real_time_step, DBL_EPSILON);
            TS_ASSERT_EQUALS(stepper.GetIdealTimeStep(), timestep);
            TS_ASSERT_EQUALS(stepper.GetTime(), current_time);
            TS_ASSERT_EQUALS(stepper.GetNextTime(), to_time);

            // Determine the new current time
            current_time = to_time;
            stepper.AdvanceOneTimeStep();

            TS_ASSERT_EQUALS(current_time, stepper.GetTime());
        }

        TS_ASSERT(stepper.IsTimeAtEnd());
        TS_ASSERT(stepper.GetTotalTimeStepsTaken()==time_step_number);

        //Stepper no longer allows increments beyond the end to be silently ignored
        TS_ASSERT_THROWS_THIS(stepper.AdvanceOneTimeStep(), "TimeStepper incremented beyond end time.");
    }

    void TestEnforceConstantTimeStep()
    {
        TimeStepper stepper(0.0, 1.0, 0.3); // timestep does not divide, but no checking

        TS_ASSERT_THROWS_THIS( TimeStepper bad_const_dt_stepper(0.0, 1.0, 0.3, true),
                "TimeStepper estimates non-constant timesteps will need to be used: "
                "check timestep divides (end_time-start_time) (or divides printing timestep). "
                "[End time=1; start=0; dt=0.3; error=0.1]");

#ifdef _MSC_VER
        _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
        TS_ASSERT_THROWS_THIS( TimeStepper bad_const_dt_stepper2(0.0, 1.0, 0.99999999, true),
                "TimeStepper estimates non-constant timesteps will need to be used: "
                "check timestep divides (end_time-start_time) (or divides printing timestep). "
                "[End time=1; start=0; dt=1; error=1e-08]");
        TimeStepper const_dt_stepper(0.0, 1.0, 0.1, true);
        unsigned counter = 0;
        while (!const_dt_stepper.IsTimeAtEnd())
        {
            counter++;
            TS_ASSERT_DELTA(const_dt_stepper.GetNextTimeStep(), 0.1, 1e-15);
            TS_ASSERT_EQUALS(const_dt_stepper.GetIdealTimeStep(), 0.1);
            const_dt_stepper.AdvanceOneTimeStep();
            if (const_dt_stepper.IsTimeAtEnd())
            {
                TS_ASSERT_EQUALS(const_dt_stepper.GetNextTimeStep(), 0.0);
            }
            else
            {
                TS_ASSERT_DELTA(const_dt_stepper.GetNextTimeStep(), 0.1, 1e-15);
            }
        }
        TS_ASSERT_EQUALS(counter,10u);
    }

    void TestAdditionalSteppingPoints()
    {
        {
            std::vector<double> additional_times_bad_order;

            additional_times_bad_order.push_back(0.7);
            additional_times_bad_order.push_back(0.2);

            TS_ASSERT_THROWS_THIS(TimeStepper stepper(0.0, 1.0, 0.1, false, additional_times_bad_order),
                                  "The additional times vector should be in ascending numerical order; "
                                  "entry 1 is less than or equal to entry 0.");
        }
        {
            std::vector<double> additional_times;
            additional_times.push_back(0.03);
            additional_times.push_back(0.25);
            additional_times.push_back(0.5);
            additional_times.push_back(0.75);

            TS_ASSERT_THROWS_THIS(TimeStepper stepper(0.0, 1.0, 0.1, false, additional_times),
                                  "Additional times are now deprecated.  Use only to check whether the given times are met: "
                                  "e.g. Electrode events should only happen on printing steps.");
      }

        std::vector<double> check_additional_times;
        check_additional_times.push_back(0.0);
        check_additional_times.push_back(0.2);
        check_additional_times.push_back(0.5);
        check_additional_times.push_back(0.7);
        TimeStepper stepper(0.0, 1.0, 0.1, false, check_additional_times);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(),10u);
        while (!stepper.IsTimeAtEnd())
        {
            stepper.AdvanceOneTimeStep();
        }

        TS_ASSERT_EQUALS(stepper.GetTotalTimeStepsTaken(),10u);
    }

    void TestResetTimeStep()
    {
        double timestep = 0.1;
        TimeStepper stepper(0.0, 0.5, timestep, false);

        TS_ASSERT_EQUALS(stepper.GetNextTimeStep(), timestep);
        TS_ASSERT_EQUALS(stepper.GetTime(), 0.0);
        TS_ASSERT_EQUALS(stepper.GetNextTime(), timestep);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), 5u);

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_EQUALS(stepper.GetTime(), 0.1);

        timestep = 0.05;
        stepper.ResetTimeStep(timestep);
        TS_ASSERT_DELTA(stepper.GetNextTimeStep(), timestep, 1e-12);
        TS_ASSERT_DELTA(stepper.GetNextTime(), 0.15, 1e-12);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), 8u);

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_DELTA(stepper.GetTime(), 0.15, 1e-12);

        timestep = 0.25;
        stepper.ResetTimeStep(timestep);
        TS_ASSERT_DELTA(stepper.GetNextTimeStep(), timestep, 1e-12);
        TS_ASSERT_DELTA(stepper.GetNextTime(), 0.4, 1e-12);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), 1u); // just an estimate

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_DELTA(stepper.GetTime(), 0.4, 1e-12);

        TS_ASSERT_DELTA(stepper.GetNextTimeStep(), 0.1, 1e-12);
        TS_ASSERT_DELTA(stepper.GetNextTime(), 0.5, 1e-12);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), 1u);

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_EQUALS(stepper.IsTimeAtEnd(), true);
    }

    void TestWithLargerStartTime()
    {
        // Abstracted from AbstractDynamicLinearPdeSolver
        TimeStepper pde_stepper(16384.1, 16384.2, 0.01, true);
        while (!pde_stepper.IsTimeAtEnd())
        {
            TS_ASSERT_DELTA(pde_stepper.GetNextTimeStep(), 0.01, DBL_EPSILON*16384.2);
            pde_stepper.AdvanceOneTimeStep();
        }
    }

    void TestWithLargeEndTime()
    {
        TimeStepper stepper(0.01*999999999, 1e7, 0.01, true);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), 1u);
        while (!stepper.IsTimeAtEnd())
        {
            TS_ASSERT_DELTA(stepper.GetNextTimeStep(), 0.01, DBL_EPSILON*1e7);
            stepper.AdvanceOneTimeStep();
        }
        TS_ASSERT_EQUALS(stepper.GetTotalTimeStepsTaken(), 1u);
    }

    void TestWithSingleTimeStep()
    {
        TimeStepper stepper(4.02, 4.04, 0.02, true);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), 1u);
        while (!stepper.IsTimeAtEnd())
        {
            TS_ASSERT_DELTA(stepper.GetNextTimeStep(), 0.02, 4.04*DBL_EPSILON);
            stepper.AdvanceOneTimeStep();
        }
        TS_ASSERT_EQUALS(stepper.GetTotalTimeStepsTaken(), 1u);
    }

    /**
     * This case used to break with the IntelProduction build, believe it or not.
     * Adding further information to the error message clearly changes the optimisations it applies,
     * and fixes the problem!
     */
    void TestIntelProductionFoolishness()
    {
        TimeStepper stepper(5.0, 10.0, 0.01, true);
        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(), 500u);
        while (!stepper.IsTimeAtEnd())
        {
            TS_ASSERT_DELTA(stepper.GetNextTimeStep(), 0.01, 10.0*DBL_EPSILON);
            TimeStepper pde_stepper(stepper.GetTime(), stepper.GetNextTime(), 0.01, true);
            TS_ASSERT_EQUALS(pde_stepper.EstimateTimeSteps(), 1u);
            while (!pde_stepper.IsTimeAtEnd())
            {
                TS_ASSERT_DELTA(pde_stepper.GetNextTimeStep(), 0.01, 10.0*DBL_EPSILON);
                pde_stepper.AdvanceOneTimeStep();
            }
            TS_ASSERT_EQUALS(pde_stepper.GetTotalTimeStepsTaken(), 1u);
            stepper.AdvanceOneTimeStep();
        }
        TS_ASSERT_EQUALS(stepper.GetTotalTimeStepsTaken(), 500u);
    }

    void TestArchiveTimeStepper()
    {
        OutputFileHandler handler("TestTimeStepper_Archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "time.arch";

        // Create and archive time-stepper
        {

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TimeStepper stepper(0.0, 1.0, 0.1, false);

            //Run for 5 steps
            for (int i=0; i<5; i++)
            {
                stepper.AdvanceOneTimeStep();
            }

            TS_ASSERT_DELTA(stepper.GetTime(), 0.5, 1e-10);
            TS_ASSERT_DELTA(stepper.GetNextTime(), 0.6, 1e-10);

            TimeStepper* const p_stepper_for_archiving = &stepper;
            output_arch << p_stepper_for_archiving;


            //Exhibit normal behaviour after archive snapshot
            stepper.AdvanceOneTimeStep();
            TS_ASSERT_DELTA(stepper.GetTime(), 0.6, 1e-10);
            TS_ASSERT_DELTA(stepper.GetNextTime(), 0.7, 1e-10);
            stepper.AdvanceOneTimeStep();
            TS_ASSERT_DELTA(stepper.GetTime(), 0.7, 1e-10);
            TS_ASSERT_DELTA(stepper.GetNextTime(), 0.8, 1e-10);

        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            TimeStepper* p_stepper;
            input_arch >> p_stepper;

            TS_ASSERT_DELTA(p_stepper->GetTime(), 0.5, 1e-10);
            TS_ASSERT_DELTA(p_stepper->GetNextTime(), 0.6, 1e-10);

            p_stepper->AdvanceOneTimeStep();
            TS_ASSERT_DELTA(p_stepper->GetTime(), 0.6, 1e-10);
            TS_ASSERT_DELTA(p_stepper->GetNextTime(), 0.7, 1e-10);

            p_stepper->AdvanceOneTimeStep();
            TS_ASSERT_DELTA(p_stepper->GetTime(), 0.7, 1e-10);
            TS_ASSERT_DELTA(p_stepper->GetNextTime(), 0.8, 1e-10);

            delete p_stepper;
        }
    }
};

#endif /*TESTTIMESTEPPER_HPP_*/
