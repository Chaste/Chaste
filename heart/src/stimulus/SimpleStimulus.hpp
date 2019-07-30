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


#ifndef _SIMPLESTIMULUS_HPP_
#define _SIMPLESTIMULUS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractStimulusFunction.hpp"

/**
 * Provides an simple stimulus of magnitude 'magnitudeOfStimulus'
 * from time 'timeOfStimulus' for duration 'duration'.
 */
class SimpleStimulus : public AbstractStimulusFunction
{
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
        archive & mTimeOfStimulus;
    }

private:
    /** The stimulus magnitude - units are:
     *    for single-cell problems - microA/cm^2
     *    for (3d) tissue simulations - microA/cm^3
     */
    double mMagnitudeOfStimulus;
    /** Duration of initial stimulus, typically in milliseconds */
    double mDuration;
    /** The time at which the stimulus starts, typically in milliseconds */
    double mTimeOfStimulus;

public:

    /**
     * Constructor.
     *
     * @param magnitudeOfStimulus  The stimulus magnitude, with units
     *    for single-cell problems - microA/cm^2
     *    for (3d) tissue simulations - microA/cm^3
     * @param duration  Duration of initial stimulus milliseconds
     * @param timeOfStimulus  The time at which the stimulus starts (defaults to 0.0) milliseconds
     */
    SimpleStimulus(double magnitudeOfStimulus, double duration, double timeOfStimulus=0.0);

    /**
     * Destructor.
     */
    virtual ~SimpleStimulus();

    /**
     * Returns the stimulus at a given time.
     *
     * @param time  time at which to return the stimulus
     */
    double GetStimulus(double time);

    /**
     * Replace the time that was specified in the constructor with a new start time.
     *
     * @param startTime
     */
    void SetStartTime(double startTime);
};


#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SimpleStimulus)

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
    Archive & ar, SimpleStimulus * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of SimpleStimulus.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */
     ::new(t)SimpleStimulus(0.0, 0.0, 0.0);
}
}
} // namespace ...

#endif //_SIMPLESTIMULUS_HPP_

