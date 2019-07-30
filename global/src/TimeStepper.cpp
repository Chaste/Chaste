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

#include "TimeStepper.hpp"
#include "Exception.hpp"
#include "MathsCustomFunctions.hpp"

TimeStepper::TimeStepper(double startTime, double endTime, double dt, bool enforceConstantTimeStep, std::vector<double> additionalTimes)
    : mStart(startTime),
      mEnd(endTime),
      mDt(dt),
      mTotalTimeStepsTaken(0),
      mAdditionalTimesReachedDeprecated(0),
      mTime(startTime),
      mEpsilon(DBL_EPSILON)
{
    if (startTime > endTime)
    {
        EXCEPTION("The simulation duration must be positive, not " << endTime-startTime);
    }

    // Remove any additionalTimes entries which fall too close to a time when the stepper would stop anyway
    for (unsigned i=0; i<additionalTimes.size(); i++)
    {
        if (i > 0)
        {
            if (additionalTimes[i-1] >= additionalTimes[i])
            {
                EXCEPTION("The additional times vector should be in ascending numerical order; "
                          "entry " << i << " is less than or equal to entry " << i-1 << ".");
            }
        }

        double time_interval = additionalTimes[i] - startTime;

        // When mDt divides this interval (and the interval is positive) then we are going there anyway
        if (!Divides(mDt, time_interval) && (time_interval > DBL_EPSILON))
        {
            //mAdditionalTimes.push_back(additionalTimes[i]);
            EXCEPTION("Additional times are now deprecated.  Use only to check whether the given times are met: e.g. Electrode events should only happen on printing steps.");
        }
    }

    /*
     * Note that when mEnd is large then the error of subtracting two numbers of
     * that magnitude is about DBL_EPSILON*mEnd (1e-16*mEnd). When mEnd is small
     * then the error should be around DBL_EPSILON.
     */
    if (mEnd > 1.0)
    {
        mEpsilon = DBL_EPSILON*mEnd;
    }

    // If enforceConstantTimeStep check whether the times are such that we won't have a variable dt
    if (enforceConstantTimeStep)
    {
        double expected_end_time = mStart + mDt*EstimateTimeSteps();

        if (fabs( mEnd - expected_end_time ) > mEpsilon)
        {
            EXCEPTION("TimeStepper estimates non-constant timesteps will need to be used: check timestep "
                      "divides (end_time-start_time) (or divides printing timestep). "
                      "[End time=" << mEnd << "; start=" << mStart << "; dt=" << mDt << "; error="
                      << fabs(mEnd-expected_end_time) << "]");
        }
    }

    mNextTime = CalculateNextTime();
}

double TimeStepper::CalculateNextTime()
{
    double next_time = mStart + (mTotalTimeStepsTaken + 1)*mDt;

    // Does the next time bring us very close to the end time?
    // Note that the inequality in this guard matches the inversion of the guard in the enforceConstantTimeStep
    // calculation of the constructor
    if (mEnd - next_time <= mEpsilon)
    {
        next_time = mEnd;
    }

    assert(mAdditionalTimesDeprecated.empty());

    return next_time;
}

void TimeStepper::AdvanceOneTimeStep()
{
    mTotalTimeStepsTaken++;
    if (mTotalTimeStepsTaken == 0)
    {
        EXCEPTION("Time step counter has overflowed.");
    }
    if (mTime >= mNextTime)
    {
        EXCEPTION("TimeStepper incremented beyond end time.");
    }
    mTime = mNextTime;

    mNextTime = CalculateNextTime();
}

double TimeStepper::GetTime() const
{
    return mTime;
}

double TimeStepper::GetNextTime() const
{
    return mNextTime;
}

double TimeStepper::GetNextTimeStep()
{
    double dt = mDt;

    if (mNextTime >= mEnd)
    {
        dt = mEnd - mTime;
    }
    assert(mAdditionalTimesDeprecated.empty());
    return dt;
}
double TimeStepper::GetIdealTimeStep()
{
    return mDt;
}

bool TimeStepper::IsTimeAtEnd() const
{
    return (mTime >= mEnd);
}

unsigned TimeStepper::EstimateTimeSteps() const
{
    return (unsigned) floor((mEnd - mStart)/mDt + 0.5);
}

unsigned TimeStepper::GetTotalTimeStepsTaken() const
{
    return mTotalTimeStepsTaken;
}

void TimeStepper::ResetTimeStep(double dt)
{
    assert(dt > 0);
    /*
     * The error in subtracting two numbers of the same magnitude is about
     * DBL_EPSILON times that magnitude (we use the sum of the two numbers
     * here as a conservative estimate of their maximum). When both mDt and
     * dt are small then the error should be around DBL_EPSILON.
     */
    double scale = DBL_EPSILON*(mDt + dt);
    if (mDt + dt < 1.0)
    {
        scale = DBL_EPSILON;
    }
    if (fabs(mDt-dt) > scale)
    {
        mDt = dt;
        mStart = mTime;
        mTotalTimeStepsTaken = 0;

        mNextTime = CalculateNextTime();
    }
}
