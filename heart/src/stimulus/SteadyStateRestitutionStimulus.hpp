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
#ifndef STEADYSTATERESTITUTIONSTIMULUS_HPP_
#define STEADYSTATERESTITUTIONSTIMULUS_HPP_

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
class SteadyStateRestitutionStimulus : public MultiStimulus
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
    SteadyStateRestitutionStimulus(){};

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
    SteadyStateRestitutionStimulus(double magnitude, double stimulusDuration, double startTime, std::vector<double> pacingCycleLengths, unsigned numberOfPulses);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SteadyStateRestitutionStimulus)

#endif /*STEADYSTATERESTITUTIONSTIMULUS_HPP_*/
