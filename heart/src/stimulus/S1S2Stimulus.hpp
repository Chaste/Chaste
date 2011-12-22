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
    S1S2Stimulus(){};

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
CHASTE_CLASS_EXPORT(S1S2Stimulus);

#endif /*S1S2STIMULUS_HPP_*/
