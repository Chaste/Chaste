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
                // std::cout<<"InternalSolve "<<timestep<<" "<<real_time_step<<"\n";
                to_time = end_time;
            }

            TS_ASSERT_EQUALS(stepper.GetNextTimeStep(), real_time_step);
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

    void TestEnforceConstantTimeStep() throw(Exception)
    {
        TimeStepper stepper(0.0, 1.0, 0.3); // timestep does not divide, but no checking

        TS_ASSERT_THROWS_THIS( TimeStepper bad_const_dt_stepper(0.0, 1.0, 0.3, true),
                "TimeStepper estimates non-constant timesteps will need to be used: "
                "check timestep divides (end_time-start_time) (or divides printing timestep). "
                "[End time=1; start=0; dt=0.3; error=0.1]");

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

    void TestAdditionalSteppingPoints() throw(Exception)
    {
        std::vector<double> additional_times_bad_order;
        additional_times_bad_order.push_back(0.75);
        additional_times_bad_order.push_back(0.25);

        TS_ASSERT_THROWS_THIS(TimeStepper stepper(0.0, 1.0, 0.1, false, additional_times_bad_order),
                              "The additional times vector should be in ascending numerical order; "
                              "entry 1 is less than or equal to entry 0.");

        std::vector<double> additional_times;
        additional_times.push_back(0.03);
        additional_times.push_back(0.25);
        additional_times.push_back(0.5);
        additional_times.push_back(0.75);

        TimeStepper stepper(0.0, 1.0, 0.1, false, additional_times);

        TS_ASSERT_EQUALS(stepper.EstimateTimeSteps(),13u);

        std::vector<double> expected_times_reverse_order;
        expected_times_reverse_order.push_back(1.0);
        expected_times_reverse_order.push_back(0.9);
        expected_times_reverse_order.push_back(0.8);
        expected_times_reverse_order.push_back(0.75);
        expected_times_reverse_order.push_back(0.7);
        expected_times_reverse_order.push_back(0.6);
        expected_times_reverse_order.push_back(0.5);
        expected_times_reverse_order.push_back(0.4);
        expected_times_reverse_order.push_back(0.3);
        expected_times_reverse_order.push_back(0.25);
        expected_times_reverse_order.push_back(0.2);
        expected_times_reverse_order.push_back(0.1);
        expected_times_reverse_order.push_back(0.03);
        expected_times_reverse_order.push_back(0.0);

        std::vector<double> expected_timesteps_reverse_order;
        expected_timesteps_reverse_order.push_back(0.1);
        expected_timesteps_reverse_order.push_back(0.1);
        expected_timesteps_reverse_order.push_back(0.05);
        expected_timesteps_reverse_order.push_back(0.05);
        expected_timesteps_reverse_order.push_back(0.1);
        expected_timesteps_reverse_order.push_back(0.1);
        expected_timesteps_reverse_order.push_back(0.1);
        expected_timesteps_reverse_order.push_back(0.1);
        expected_timesteps_reverse_order.push_back(0.05);
        expected_timesteps_reverse_order.push_back(0.05);
        expected_timesteps_reverse_order.push_back(0.1);
        expected_timesteps_reverse_order.push_back(0.07);
        expected_timesteps_reverse_order.push_back(0.03);

        while (!stepper.IsTimeAtEnd())
        {
            TS_ASSERT_DELTA(stepper.GetTime(),expected_times_reverse_order.back(),1e-12);
            expected_times_reverse_order.pop_back();

            TS_ASSERT_DELTA(stepper.GetNextTimeStep(),expected_timesteps_reverse_order.back(),1e-12);
            expected_timesteps_reverse_order.pop_back();

            stepper.AdvanceOneTimeStep();
        }

        TS_ASSERT_EQUALS(stepper.GetTotalTimeStepsTaken(),13u);
    }

    void TestResetTimeStep() throw(Exception)
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

    void TestResetTimeStepWithAdditionalSteppingPoints() throw(Exception)
    {
        double timestep = 0.1;
        std::vector<double> additional_times;
        additional_times.push_back(0.33);

        TimeStepper stepper(0.0, 0.5, timestep, false, additional_times);

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_EQUALS(stepper.GetTime(), 0.1);

        timestep = 0.25;
        stepper.ResetTimeStep(timestep);
        TS_ASSERT_DELTA(stepper.GetNextTimeStep(), 0.23, 1e-12);
        TS_ASSERT_DELTA(stepper.GetNextTime(), 0.33, 1e-12);

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_EQUALS(stepper.GetTime(), 0.33);

        TS_ASSERT_DELTA(stepper.GetNextTimeStep(), 0.02, 1e-12);
        TS_ASSERT_DELTA(stepper.GetNextTime(), 0.35, 1e-12);

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_EQUALS(stepper.GetTime(), 0.35);

        TS_ASSERT_DELTA(stepper.GetNextTimeStep(), 0.15, 1e-12);
        TS_ASSERT_DELTA(stepper.GetNextTime(), 0.5, 1e-12);

        stepper.AdvanceOneTimeStep();
        TS_ASSERT_EQUALS(stepper.GetTime(), 0.5);

        TS_ASSERT_EQUALS(stepper.IsTimeAtEnd(), true);
    }

    void TestWithLargerStartTime() throw(Exception)
    {
        // Abstracted from AbstractDynamicLinearPdeSolver
        TimeStepper pde_stepper(16384.1, 16384.2, 0.01, true);
        while (!pde_stepper.IsTimeAtEnd())
        {
            TS_ASSERT_DELTA(pde_stepper.GetNextTimeStep(), 0.01, DBL_EPSILON*16384.2);
            pde_stepper.AdvanceOneTimeStep();
        }
    }

    void TestWithLargeEndTime() throw(Exception)
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

    void TestWithSingleTimeStep() throw(Exception)
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
    void TestIntelProductionFoolishness() throw(Exception)
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

            std::vector<double> additional_times;
            additional_times.push_back(0.55);
            TimeStepper stepper(0.0, 1.0, 0.1, false, additional_times);

            //Run for 5 steps
            for (int i=0; i<5; i++)
            {
                stepper.AdvanceOneTimeStep();
            }

            TS_ASSERT_DELTA(stepper.GetTime(), 0.5, 1e-10);
            TS_ASSERT_DELTA(stepper.GetNextTime(), 0.55, 1e-10);

            TimeStepper* const p_stepper_for_archiving = &stepper;
            output_arch << p_stepper_for_archiving;


            //Exhibit normal behaviour after archive snapshot
            stepper.AdvanceOneTimeStep();
            TS_ASSERT_DELTA(stepper.GetTime(), 0.55, 1e-10);
            TS_ASSERT_DELTA(stepper.GetNextTime(), 0.6, 1e-10);
            stepper.AdvanceOneTimeStep();
            TS_ASSERT_DELTA(stepper.GetTime(), 0.6, 1e-10);
            TS_ASSERT_DELTA(stepper.GetNextTime(), 0.7, 1e-10);

        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            TimeStepper* p_stepper;
            input_arch >> p_stepper;

            TS_ASSERT_DELTA(p_stepper->GetTime(), 0.5, 1e-10);
            TS_ASSERT_DELTA(p_stepper->GetNextTime(), 0.55, 1e-10);

            p_stepper->AdvanceOneTimeStep();
            TS_ASSERT_DELTA(p_stepper->GetTime(), 0.55, 1e-10);
            TS_ASSERT_DELTA(p_stepper->GetNextTime(), 0.6, 1e-10);

            p_stepper->AdvanceOneTimeStep();
            TS_ASSERT_DELTA(p_stepper->GetTime(), 0.6, 1e-10);
            TS_ASSERT_DELTA(p_stepper->GetNextTime(), 0.7, 1e-10);

            delete p_stepper;
        }
    }

};

#endif /*TESTTIMESTEPPER_HPP_*/
