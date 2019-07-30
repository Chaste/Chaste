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

#include <cassert>
#include <cmath>
#include "SimulationTime.hpp"

/** Pointer to the single instance */
SimulationTime* SimulationTime::mpInstance = nullptr;

/** Shared pointer to the delegated class */
boost::shared_ptr<TimeStepper> SimulationTime::mpTimeStepper;

SimulationTime* SimulationTime::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new SimulationTime;
        mpTimeStepper.reset();
        std::atexit(Destroy);
    }
    return mpInstance;
}

SimulationTime::SimulationTime()
    :
      mStartTime(DOUBLE_UNSET)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == nullptr);
}

void SimulationTime::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

double SimulationTime::GetTimeStep() const
{
    assert(mpTimeStepper);
    return mpTimeStepper->GetIdealTimeStep();
}

void SimulationTime::IncrementTimeOneStep()
{
    assert(mpTimeStepper);
    mpTimeStepper->AdvanceOneTimeStep(); //This can now throw if the end time has been reached
}

unsigned SimulationTime::GetTimeStepsElapsed() const
{
    assert(mpTimeStepper);
    return mpTimeStepper->GetTotalTimeStepsTaken();
}

double SimulationTime::GetTime() const
{
    // NOTE: if this assertion fails, it may be because Destroy() wasn't called in the previous test
    assert(mStartTime != DOUBLE_UNSET);
    //Check if the time stepping has started
    if (mpTimeStepper)
    {
        return mpTimeStepper->GetTime();
    }
    //If time stepping hasn't started then we are still at start time
    return mStartTime;
}

void SimulationTime::SetStartTime(double startTime)
{
    assert(mStartTime == DOUBLE_UNSET);
    mStartTime = startTime;
}

void SimulationTime::SetEndTimeAndNumberOfTimeSteps(double endTime, unsigned totalTimeStepsInSimulation)
{
    // NOTE: if this assertion fails, it may be because Destroy() wasn't called in the previous test
    assert(mStartTime != DOUBLE_UNSET);
    assert(!mpTimeStepper);
    assert(endTime > mStartTime);
    assert(totalTimeStepsInSimulation != 0u);

    mpTimeStepper.reset(new TimeStepper(mStartTime, endTime, (endTime-mStartTime)/totalTimeStepsInSimulation, true));
}

void SimulationTime::ResetEndTimeAndNumberOfTimeSteps(const double& rEndTime, const unsigned& rNumberOfTimeStepsInThisRun)
{
    // NOTE: if this assertion fails, it may be because Destroy() wasn't called in the previous test
    assert(mStartTime != DOUBLE_UNSET);
    // NOTE: If this assertion fails, you should be using set rather than reset
    assert(mpTimeStepper);
    mStartTime = mpTimeStepper->GetTime();

    assert(rEndTime > mStartTime);

    // Reset the machinery that works out the time
    mpTimeStepper.reset(new TimeStepper(mStartTime, rEndTime, (rEndTime-mStartTime)/rNumberOfTimeStepsInThisRun, true));
}

bool SimulationTime::IsStartTimeSetUp() const
{
    return (mStartTime != DOUBLE_UNSET);
}

bool SimulationTime::IsEndTimeAndNumberOfTimeStepsSetUp() const
{
    if (mpTimeStepper)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool SimulationTime::IsFinished() const
{
    // assert(mpTimeStepper);
    return(mpTimeStepper->IsTimeAtEnd());
}

