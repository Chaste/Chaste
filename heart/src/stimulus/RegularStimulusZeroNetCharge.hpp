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


#ifndef _REGULARSTIMULUSZERONETCHARGE_HPP_
#define _REGULARSTIMULUSZERONETCHARGE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "RegularStimulus.hpp"

/**
 * Provides a periodic square-wave stimulus, where the total net charge is zero for every stimulus.
 *
 * i.e.,
 *         --------
 *         |      |
 *         |      |
 *  --------      |      -------------
 *                |      |
 *                |      |
 *                --------
 *
 *
 *          <------------>
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
