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

#include <cassert>
#include <cmath>
#include "SimulationTime.hpp"

/** Pointer to the single instance */
SimulationTime* SimulationTime::mpInstance = NULL;

/** Shared pointer to the delegated class */
boost::shared_ptr<TimeStepper> SimulationTime::mpTimeStepper;

SimulationTime* SimulationTime::Instance()
{
    if (mpInstance == NULL)
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
    assert(mpInstance == NULL);
}

void SimulationTime::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
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
    if(mpTimeStepper)
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

