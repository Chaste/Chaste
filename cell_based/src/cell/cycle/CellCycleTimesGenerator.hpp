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

#ifndef CELLCYCLETIMESGENERATOR_HPP_
#define CELLCYCLETIMESGENERATOR_HPP_

#include <sstream>

#include "ChasteSerialization.hpp"
#include "RandomNumberGenerator.hpp"
#include "SerializableSingleton.hpp"

/**
 * This is a helper class for FixedSequenceCellCycleModel. It comprises a singleton class
 * for producing a fixed stream of exponentially distributed G1 durations.
 */
class CellCycleTimesGenerator : public SerializableSingleton<CellCycleTimesGenerator>
{
private:
    /** The singleton definition
     */
    static CellCycleTimesGenerator* mpInstance;

    /**
     * The random seed used to generate the series of cell cycle times.
     */
    unsigned mRandomSeed;

    /**
     * The actual sequence of cell cycle times
     */
    std::vector<double> mCellCycleTimes;

    /**
     * This index is used to keep track of where in the sequence mCellCycleTimes
     * the simulation currently is. Whenever a cell divides, a new
     * cell cycle time gets drawn from mCellCycleTimes and this counter is incremented
     * by one.
     */
    unsigned mCurrentIndex;

    /**
     * The rate constant of the exponential distribution used to
     * generate mCellCycleTimes.
     */
    double mRate;

    /**
     * A flag that changes from False to True when mCellCycleTimes has been
     * generated.
     */
    bool mVectorCreated;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialization of a CellCycleTimesGenerator object must be done with care.
     * Do not serialize this singleton directly.  Instead, serialize
     * the object returned by GetSerializationWrapper.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& mRandomSeed;
        archive& mCellCycleTimes;
        archive& mCurrentIndex;
        archive& mRate;
        archive& mVectorCreated;
    }

protected:
    /**
     * Protected constructor.
     * Use Instance() to access the random number generator.
     */
    CellCycleTimesGenerator();

public:
    /**
     * @return a pointer to the CellCycleTimesGenerator object.
     * The object is created the first time this method is called.
     */
    static CellCycleTimesGenerator* Instance();

    /**
     * Destroy the current instance of the cell cycle times generator.
     * The next call to Instance will create a new instance.
     * This method *must* be called before program exit, to avoid a memory
     * leak.
     */
    static void Destroy();

    /**
     * Set the random seed for generating the sequence of G1 Durations.
     *
     * @param randomSeed the new seed
     */
    void SetRandomSeed(unsigned randomSeed);

    /**
     * Get the random seed for generating the sequence.
     *
     * @returns seed The random seed
     */
    unsigned GetRandomSeed();

    /**
     * Set the rate of the exponential distribution.
     *
     * @param rate
     */
    void SetRate(double rate);

    /**
     * Get the rate of the exponential distribution.
     *
     * @returns rate The rate of the distribution.
     */
    double GetRate();

    /**
     * Generate the cell cycle time sequence
     */
    void GenerateCellCycleTimeSequence();

    /**
     * Get the next cell cycle time
     *
     * @returns next_cell_cycle_time  the next cell cycle time in the sequence
     */
    double GetNextCellCycleTime();
};

#endif /*CELLCYCLETIMESGENERATOR_HPP_*/
