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

#ifndef S1S2STIMULUS_HPP_
#define S1S2STIMULUS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "MultiStimulus.hpp"
#include "RegularStimulus.hpp" // For archiving this needs to be here.

/**
 * This class provides a stimulus function which
 * follows an S1-S2 protocol. i.e.
 *
 * Run a RegularStimulus at a certain (S1) frequency
 * Then run two pulses of a second RegularStimulus at
 * another (S2) frequency. These are combined into a
 * MultiStimulus.
 *
 * We have a number of these MultiStimulus objects, created
 * for a given vector s2Periods in the constructor. As a
 * user you set which of these S2 frequencies to run with the method
 * SetS2ExperimentPeriodIndex().
 *
 */
class S1S2Stimulus : public MultiStimulus
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
        archive & mS2Index;
        archive & mNumS2FrequencyValues;
    }

    /**
     * The Index of the S2 Frequency that we wish to run (as specified in the vector in the constructor).
     */
    unsigned mS2Index;
    /**
     * The Number of S2 Frequencies that have been set up
     */
    unsigned mNumS2FrequencyValues;

    /**
     * Private constructor - for archiving's eyes only.
     */
    S1S2Stimulus() {};

public:

    /**
     * Constructor
     *
     * @param magnitude  The magnitude of the stimulus 'square wave'.
     * @param stimulusDuration  The duration of the stimulus 'square wave'.
     * @param s1Duration  The duration of the S1 phase of the experiment, should be a multiple of the s1Period.
     * @param s1Period  The period of the S1 phase (typically ~1000ms for example).
     * @param startTime  The time at which to begin the S1 stimulus (this delay is automatically added to give smooth transition to S2).
     * @param s2Periods  A std::vector of the different S2 frequencies that should be run (for example {1000, 900, 800} (ms)).
     *
     */
    S1S2Stimulus(double magnitude, double stimulusDuration, double s1Duration, double s1Period, double startTime, std::vector<double> s2Periods);

    /**
     * Get the magnitude of the multiple stimuli at time 'time'
     *
     * @param time  time at which to return the stimulus
     * @return  Magnitude of stimulus at time 'time'.
     */
     double GetStimulus(double time);

     /**
      * Allows us to move to the 'next' S2 frequency.
      *
      * @param index  Which S2 frequency to use next.
      */
     void SetS2ExperimentPeriodIndex(unsigned index);

     /**
      * @return The number of different frequencies this stimulus can perform.
      */
     unsigned GetNumS2FrequencyValues();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(S1S2Stimulus)

#endif /*S1S2STIMULUS_HPP_*/
