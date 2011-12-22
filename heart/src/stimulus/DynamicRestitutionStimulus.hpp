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
#ifndef DYNAMICRESTITUTIONSTIMULUS_HPP_
#define DYNAMICRESTITUTIONSTIMULUS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "RegularStimulus.hpp" // For archiving this needs to be here.
#include "MultiStimulus.hpp"

/**
 * This class provides a stimulus function which
 * follows a dynamic restitution protocol. i.e.
 *
 * Run a RegularStimulus at certain frequencies for a
 * certain number of pulses. These are combined into a
 * MultiStimulus.
 *
 */
class DynamicRestitutionStimulus : public MultiStimulus
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
        archive & boost::serialization::base_object<MultiStimulus>(*this);
    }

    /**
     * Private constructor - for archiving's eyes only.
     */
    DynamicRestitutionStimulus(){};

public:

    /**
     * Constructor
     *
     * @param magnitude  The magnitude of the stimulus 'square wave'.
     * @param stimulusDuration  The duration of the stimulus 'square wave'.
     * @param startTime  The time at which to begin the S1 stimulus (this delay is automatically added to give smooth transition to S2).
     * @param pacingCycleLengths A vector containing the pacing cycle lengths (in ms) of each phase of the protocol
     * @param numberOfPulses  The number of pulses to perform at each pacing cycle length
     */
    DynamicRestitutionStimulus(double magnitude, double stimulusDuration, double startTime, std::vector<double> pacingCycleLengths, unsigned numberOfPulses);

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(DynamicRestitutionStimulus);

#endif /*DYNAMICRESTITUTIONSTIMULUS_HPP_*/
