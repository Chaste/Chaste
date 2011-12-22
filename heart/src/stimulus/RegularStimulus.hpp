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
     * set the pacing cycle length ('period') of the stimulus.
     *
     * @param period  The stimulus pacing cycle length to use.
     */
    void SetPeriod(double period);

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
