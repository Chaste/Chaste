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

