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


#include "RegularStimulusZeroNetCharge.hpp"
#include <cmath>

RegularStimulusZeroNetCharge::RegularStimulusZeroNetCharge(double magnitudeOfStimulus, double duration, double period, double startTime, double stopTime)
    : RegularStimulus(magnitudeOfStimulus, duration, period, startTime,stopTime)
{

}

double RegularStimulusZeroNetCharge::GetStimulus(double time)
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

    double beatTime = fmod(time-mStartTime,mPeriod);

    if (beatTime >=0 && beatTime <= mDuration/2 && time <= mStopTime)
    {
        return mMagnitudeOfStimulus;
    }
    else if (beatTime >=0 && beatTime > mDuration/2 && beatTime <= mDuration && time <= mStopTime)
    {
        return -mMagnitudeOfStimulus;
    }
    else
    {
        return 0.0;
    }
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RegularStimulusZeroNetCharge)
