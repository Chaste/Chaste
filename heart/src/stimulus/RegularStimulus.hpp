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


#ifndef _REGULARSTIMULUS_HPP_
#define _REGULARSTIMULUS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractStimulusFunction.hpp"

/**
 * Provides a periodic square-wave stimulus.
 */
class RegularStimulus : public AbstractStimulusFunction
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simple stimulus, never used directly - boost uses this.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractStimulusFunction>(*this);
        archive & mMagnitudeOfStimulus;
        archive & mDuration;
        archive & mPeriod;
        archive & mStartTime;
        archive & mStopTime;
    }

protected:

    /** The 'height' of the square wave applied */
    double mMagnitudeOfStimulus;
    /** The length of the square wave */
    double mDuration;
    /** The time between applications of the wave */
    double mPeriod;
    /** The time at which the first wave is applied */
    double mStartTime;
    /** The time at which all stimuli are removed (even if halfway through a wave)*/
    double mStopTime;

public:
    /**
     * Create a new stimulus.
     *
     * @param magnitudeOfStimulus  The size of the stimulus
     * @param duration  How long the square wave is applied for
     * @param period  The time between square waves being applied
     * @param startTime  The time at which the first wave is applied
     * @param stopTime  The time the stimulus is removed (defaults to DBL_MAX if omitted)
     */
    RegularStimulus(double magnitudeOfStimulus, double duration, double period, double startTime, double stopTime=DBL_MAX);

    /**
     * Get the magnitude of stimulus at time 'time'
     *
     * @param time  The current time
     * @return  Magnitude of stimulus
     */
    double GetStimulus(double time);

    /**
     * @return the pacing cycle length or period of the stimulus.
     */
    double GetPeriod();

    /**
     * @return the height of the stimulus square wave (magnitude of current).
     */
    double GetMagnitude();

    /**
     * @return the duration of the stimulus square wave.
     */
    double GetDuration();

    /**
     * @return the start time of the stimulus square wave.
     */
    double GetStartTime();

    /**
     * Set the magnitude of the stimulus to apply.
     * Takes units of uA/cm^2 in single cell simulations
     * or uA/cm^3 in tissue simulations
     *
     * See https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/ChasteUnits for a full discussion of this.
     *
     * @param magnitude  The magnitude of stimulus to apply.
     */
    void SetMagnitude(double magnitude);

    /**
     * set the pacing cycle length ('period') of the stimulus.
     *
     * @param period  The stimulus pacing cycle length to use.
     */
    void SetPeriod(double period);

    /**
     * set the length ('duration') of the stimulus square wave.
     *
     * @param duration  The stimulus duration to use.
     */
    void SetDuration(double duration);

    /**
     * Set the stimulus to start at a particular time.
     *
     * @param startTime the time the stimulus should begin.
     */
    void SetStartTime(double startTime);

    /**
     * Set the stop time for this stimulus. It will never be applied after this time.
     *
     * @param stopTime
     */
    void SetStopTime(double stopTime);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(RegularStimulus)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SimpleStimulus instance (using existing constructor)
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, RegularStimulus * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of SimpleStimulus.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */
     ::new(t)RegularStimulus(0.0, 0.0, 0.1, 0.0, 1.0);
}
}
} // namespace ...

#endif //_REGULARSTIMULUS_HPP_
