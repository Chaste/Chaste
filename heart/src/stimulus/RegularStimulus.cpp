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


#include "RegularStimulus.hpp"
#include <cmath>
#include <cassert>

#include <iostream>

RegularStimulus::RegularStimulus(double magnitudeOfStimulus, double duration, double period, double startTime, double stopTime)
{
    assert(period > 0);
    assert(period >= duration);

    mMagnitudeOfStimulus = magnitudeOfStimulus;
    mDuration = duration;
    mPeriod = period;
    mStartTime = startTime;
    mStopTime = stopTime;

    //Swell duration to avoid rounding issues
    mDuration += period*DBL_EPSILON;
}

double RegularStimulus::GetStimulus(double time)
{
    /*
     * This covers a borderline case where start time is larger than period
     * and period divides start time, e.g. start_time=2 period=1.
     *
     * fmod() will return -0 in this case and will therefore switch the stimulus on.
     */
    if (time < mStartTime)
    {
        return 0.0;
    }

    double beatTime = fmod(time-mStartTime, mPeriod);

    if (beatTime >=0 && beatTime <= mDuration && time <= mStopTime)
    {
        return mMagnitudeOfStimulus;
    }
    else
    {
        return 0.0;
    }
}

double RegularStimulus::GetPeriod()
{
    return mPeriod;
}

double RegularStimulus::GetMagnitude()
{
    return mMagnitudeOfStimulus;
}

double RegularStimulus::GetDuration()
{
    return mDuration;
}

double RegularStimulus::GetStartTime()
{
    return mStartTime;
}

void RegularStimulus::SetMagnitude(double magnitude)
{
    mMagnitudeOfStimulus = magnitude;
}

void RegularStimulus::SetPeriod(double period)
{
    mPeriod = period;
}

void RegularStimulus::SetDuration(double duration)
{
    mDuration = duration;
}

void RegularStimulus::SetStartTime(double startTime)
{
    mStartTime = startTime;
}

void RegularStimulus::SetStopTime(double stopTime)
{
    mStopTime = stopTime;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RegularStimulus)
