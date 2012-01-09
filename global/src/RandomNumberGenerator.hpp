/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef RANDOMNUMBERGENERATORS_HPP_
#define RANDOMNUMBERGENERATORS_HPP_

#include "ChasteSerialization.hpp"
#include "SerializableSingleton.hpp"
#include <boost/serialization/split_member.hpp>

/**
 * A special singleton class allowing one to generate different types of
 * random number in a globally consistent way.
 */
class RandomNumberGenerator : public SerializableSingleton<RandomNumberGenerator>
{
private:

    /** The random number generator seed. */
    int mSeed;

    /** The number of times the random number generator has been called. */
    long unsigned mTimesCalled;

    /** Pointer to the single instance. */
    static RandomNumberGenerator* mpInstance;

    /** Working memory for the normal random number calculations. */
    double mWorkingMem1;

    /** Working memory for the normal random number calculations. */
    double mWorkingMem2;

    /** Working memory for the normal random number calculations. */
    double mWorkingMemW;

    /** Stored normal random number (we generate two at once). */
    double mRandNum2;

    /** Whether to use the stored random number next. */
    bool mUseLastNum;

    friend class boost::serialization::access;
    /**
     * Save the RandomNumberGenerator and its member variables.
     *
     * @note do not serialize this singleton directly.  Instead, serialize
     * the object returned by GetSerializationWrapper().
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & mSeed;
        archive & mTimesCalled;
        archive & mUseLastNum;
        archive & mRandNum2;
    }

    /**
     * Load the RandomNumberGenerator and its member variables.
     *
     * @note do not serialize this singleton directly.  Instead, serialize
     * the object returned by GetSerializationWrapper().
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & mSeed;
        archive & mTimesCalled;
        archive & mUseLastNum;
        archive & mRandNum2;

        // Reset the random number generator to use the correct seed
        srandom(mSeed);

        /*
         * Call it the correct number of times to put it in the same state
         * as it used to be.
         */
        for (long unsigned i=0; i<mTimesCalled; i++)
        {
            random();
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected:

    /**
     * Protected constructor.
     * Use Instance() to access the random number generator.
     */
    RandomNumberGenerator();

public:

    /**
     * Uses a ranf() call and an approximation to the inverse PDF of the
     * normal distribution to generate an accurate estimate of a normally
     * distributed random number. Adapted from the algorithm outlined at
     * http://home.online.no/~pjacklam/notes/invnorm/#Pseudo_code_for_rational_approximation
     *
     * @return a random number from the normal distribution with mean 0
     * and standard distribution 1.
     */
    double StandardNormalRandomDeviate();

    /**
     * Generate a random number from a normal distribution with given
     * mean and standard deviation.
     *
     * @param mean the mean of the normal distribution from which the random number is drawn
     * @param sd the standard deviation of the normal distribution from which the random number is drawn
     */
    double NormalRandomDeviate(double mean, double sd);

    /**
     * Generate a uniform random number in (0,1).
     */
    double ranf();

    /**
     * Generate a random number modulo base (i.e. an integer
     * within the range 0,..,base-1),
     *
     * @param base the order of the field of positive integers from which the random number is drawn
     */
    unsigned randMod(unsigned base);

    /**
     * Produce a permutation of the integers 0,1,..,num-1, using the Knuth-algorithm
     * (also called the Fisher-Yates algorithm), a linear time unbiased method.
     * The shuffled values are returned in rValues, which doesn't need to
     * be correctly-sized when passed in.
     *
     * @param num  the number of integers to shuffle
     * @param rValues  the output permutation of shuffled values (any initial values ignored)
     */
    void Shuffle(unsigned num, std::vector<unsigned>& rValues);

    /**
     * Return a pointer to the random number generator object.
     * The object is created the first time this method is called.
     */
    static RandomNumberGenerator* Instance();

    /**
     * Destroy the current instance of the random number generator.
     * The next call to Instance will create a new instance and re-seed.
     * This method *must* be called before program exit, to avoid a memory
     * leak.
     */
    static void Destroy();

    /**
     * Reseed the random number generator.
     *
     * @param seed the new seed
     */
    void Reseed(int seed);
};

#endif /*RANDOMNUMBERGENERATORS_HPP_*/
