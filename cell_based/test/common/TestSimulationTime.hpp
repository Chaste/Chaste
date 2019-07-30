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

#ifndef TESTSIMULATIONTIME_HPP_
#define TESTSIMULATIONTIME_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <fstream>

#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "PetscTools.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSimulationTime : public CxxTest::TestSuite
{
public:
    void TestTime()
    {
        // Create the simulation time object
        // Set the simulation length and number of time steps
        SimulationTime* p_simulation_time = SimulationTime :: Instance();

        TS_ASSERT_EQUALS(p_simulation_time->IsStartTimeSetUp(), false);

        p_simulation_time->SetStartTime(0.0);

        TS_ASSERT_EQUALS(p_simulation_time->IsStartTimeSetUp(), true);

        TS_ASSERT_EQUALS(p_simulation_time->IsEndTimeAndNumberOfTimeStepsSetUp(), false);

        //Should be able to get time before setting up the rest of the machinery
        TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.0, 1e-6);

        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 3);

        TS_ASSERT_EQUALS(p_simulation_time->IsEndTimeAndNumberOfTimeStepsSetUp(), true);

        // Get the time step
        TS_ASSERT_DELTA(p_simulation_time->GetTimeStep(), 3.33333333, 1e-6);

        // Get a second instance
        // Check that the time step is set correctly
        SimulationTime* p_simulation_time2 = SimulationTime :: Instance();
        TS_ASSERT_DELTA(p_simulation_time2->GetTimeStep(), 3.33333333, 1e-6);

        // Check that number of time steps starts at 0
        TS_ASSERT_EQUALS(p_simulation_time->GetTimeStepsElapsed(), 0u);

        // Increment the time
        p_simulation_time->IncrementTimeOneStep();

        // Check the number of time steps
        TS_ASSERT_EQUALS(p_simulation_time->GetTimeStepsElapsed(), 1u);

        // Check the simulation time from the second instance
        TS_ASSERT_DELTA(p_simulation_time2->GetTime(), 3.33333333, 1e-6);

        // Increment the time twice
        TS_ASSERT(!(p_simulation_time->IsFinished()));
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT(!(p_simulation_time->IsFinished()));
        p_simulation_time->IncrementTimeOneStep();

        // Check the simulation time from the first instance
        TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 10.0);
        TS_ASSERT(p_simulation_time->IsFinished());

        SimulationTime::Destroy();

        SimulationTime* p_simulation_time3 = SimulationTime :: Instance();
        p_simulation_time3->SetStartTime(0.0);
        p_simulation_time3->SetEndTimeAndNumberOfTimeSteps(10.0,5);
        TS_ASSERT_DELTA(p_simulation_time3->GetTimeStep(), 2.0, 1e-6);

        SimulationTime::Destroy();

        p_simulation_time3 = SimulationTime :: Instance();
        p_simulation_time3->SetStartTime(5.0);
        p_simulation_time3->SetEndTimeAndNumberOfTimeSteps(10.0,5);
        TS_ASSERT_DELTA(p_simulation_time3->GetTimeStep(), 1.0, 1e-6);

        SimulationTime::Destroy();
    }

    void TestSimulationTestDoesNotRunOver()
    {
        SimulationTime* p_simulation_time = SimulationTime :: Instance();
        p_simulation_time->SetStartTime(0.0);

        TS_ASSERT_EQUALS(p_simulation_time->IsStartTimeSetUp(), true);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 10);
        double time_step = p_simulation_time->GetTimeStep();
        TS_ASSERT_DELTA(time_step, 1.0, 1e-6);

        for (unsigned i=0; i<10; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_LESS_THAN_EQUALS(p_simulation_time->GetTime(), 10);
        }
        TS_ASSERT_THROWS_THIS(p_simulation_time->IncrementTimeOneStep(), "TimeStepper incremented beyond end time.");
        TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 10);

        SimulationTime::Destroy();
    }

    void TestResetTime()
    {
        // Create the simulation time object
        // Set the simulation length and number of time steps
        SimulationTime* p_simulation_time = SimulationTime :: Instance();
        p_simulation_time->SetStartTime(0.0);
        unsigned num_steps = 4;
        double first_end = 10.0;
        double second_end = 20.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(first_end, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            double time_should_be = i*first_end/(double)num_steps;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), time_should_be, 1e-9);
            p_simulation_time->IncrementTimeOneStep();
        }

        // Reset the end time and number of steps
        num_steps = 20;
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(second_end, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            double time_should_be = first_end + i*(second_end-first_end)/(double)num_steps;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), time_should_be, 1e-9);
            p_simulation_time->IncrementTimeOneStep();
        }

        TS_ASSERT_DELTA(p_simulation_time->GetTime(), second_end, 1e-9);

        SimulationTime::Destroy();
    }

    void TestLoopedStepping()
    {
        SimulationTime* p_simulation_time = SimulationTime :: Instance();
        double start = 0.25833333333333330373 - DBL_EPSILON;
        double end = 0.5;
        double dt = 0.0083333333333333332177;
        //As in the abstract simulation classes...
        unsigned num_time_steps = (unsigned) ((end-start)/dt+0.5);
        TS_ASSERT_EQUALS(num_time_steps, 29U);
        p_simulation_time->SetStartTime(start);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end, num_time_steps);
        TS_ASSERT_DELTA( p_simulation_time->GetTimeStep(), (1.0/120.0), DBL_EPSILON);

        unsigned step = 0;
        while (!p_simulation_time->IsFinished())
        {
            double time_should_be = start + step*dt;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), time_should_be, 1e-9);
            p_simulation_time->IncrementTimeOneStep();
            step++;
        }
        TS_ASSERT_DELTA(p_simulation_time->GetTime(), end, DBL_EPSILON);

        SimulationTime::Destroy();
    }

    void TestArchiveSimulationTime()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "time.arch";

        // Create and archive simulation time
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);
            p_simulation_time->IncrementTimeOneStep();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(2.0,6);
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.75, 1e-9);

            SerializableSingleton<SimulationTime>* const p_wrapper = p_simulation_time->GetSerializationWrapper();
            output_arch << p_wrapper;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.75, 1e-9);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 1.0, 1e-9);

            SimulationTime::Destroy();
        }

        // Restore
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(-100.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(5.0, 5);
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), -100.0, 1e-9);

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SerializableSingleton<SimulationTime>* p_wrapper;
            input_arch >> p_wrapper;

            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.75,1e-9);
            TS_ASSERT_DELTA(p_simulation_time->GetTimeStep(), 0.25, 1e-9);

            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 1.0, 1e-9);

            SimulationTime::Destroy();
        }
    }
};

#endif /*TESTSIMULATIONTIME_HPP_*/
