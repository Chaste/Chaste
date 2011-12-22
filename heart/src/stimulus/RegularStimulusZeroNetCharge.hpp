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


#ifndef _REGULARSTIMULUSZERONETCHARGE_HPP_
#define _REGULARSTIMULUSZERONETCHARGE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "RegularStimulus.hpp"

/**
 * Provides a periodic square-wave stimulus, where the total net charge is zero for every stimulus.
 *
 * i.e.,
 *               ------
 *           |        |
 *              |        |
 *  ---------       |     -------------
 *                  |     |
 *                  |     |
 *                  ------
 *
 *
 *           <----------->
 *              Duration
 */
class RegularStimulusZeroNetCharge : public RegularStimulus
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
        archive & boost::serialization::base_object<RegularStimulus>(*this);
    }

public:
    /**
     * Create a new stimulus.
     *
     * @param magnitudeOfStimulus  The size of the stimulus, both positive and negative
     * @param duration  How long the stimulus is applied for (time positive + time negative)
     * @param period  The time between square waves being applied
     * @param startTime  The time at which the first wave is applied
     * @param stopTime  The time the stimulus is removed (defaults to DBL_MAX if omitted)
     */
    RegularStimulusZeroNetCharge(double magnitudeOfStimulus, double duration, double period, double startTime, double stopTime=DBL_MAX);

    /**
     * Get the magnitude of stimulus at time 'time'. Re-implemented from parent class RegulsrStimulus.
     *
     * @param time  The current time
     * @return  Magnitude of stimulus at time 'time'.
     */
    double GetStimulus(double time);

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(RegularStimulusZeroNetCharge)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a RegularStimulusZeroNetCharge instance (using existing constructor)
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, RegularStimulusZeroNetCharge * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of SimpleStimulus.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */
     ::new(t)RegularStimulusZeroNetCharge(0.0, 0.0, 0.1, 0.0, 1.0);
}
}
} // namespace ...

#endif //_REGULARSTIMULUS_HPP_
