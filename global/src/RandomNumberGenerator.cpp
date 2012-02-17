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

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>

#include "RandomNumberGenerator.hpp"

RandomNumberGenerator* RandomNumberGenerator::mpInstance = NULL;

RandomNumberGenerator::RandomNumberGenerator()
    : mSeed(0),
      mTimesCalled(0),
      mWorkingMem1(0),
      mWorkingMem2(0),
      mWorkingMemW(0),
      mRandNum2(0),
      mUseLastNum(false)
{
    assert(mpInstance == NULL); // Ensure correct serialization
    srandom(mSeed);
}

RandomNumberGenerator* RandomNumberGenerator::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new RandomNumberGenerator();
    }
    return mpInstance;
}

void RandomNumberGenerator::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

unsigned RandomNumberGenerator::randMod(unsigned base)
{
    mTimesCalled++;
    return (random()%base);
}

double RandomNumberGenerator::ranf()
{
    mTimesCalled++;
    return (double)random() / RAND_MAX;
}

double RandomNumberGenerator::NormalRandomDeviate(double mean, double sd)
{
    return sd * StandardNormalRandomDeviate() + mean;
}

void RandomNumberGenerator::Reseed(int seed)
{
    mSeed = seed;
    srandom(mSeed);
    mTimesCalled = 0;
}

void RandomNumberGenerator::Shuffle(unsigned num, std::vector<unsigned>& rValues)
{
	rValues.resize(num);
    for (unsigned i=0; i<num; i++)
    {
        rValues[i] = i;
    }

    for (unsigned end=num-1; end>0; end--)
    {
        // Pick a random integer from {0,..,end}
        unsigned k = RandomNumberGenerator::Instance()->randMod(end+1);
        unsigned temp = rValues[end];
        rValues[end] = rValues[k];
        rValues[k] = temp;
    }
}

double  RandomNumberGenerator::StandardNormalRandomDeviate()
{
    /**
     * This little method implements the Box-Mueller/Marsaglia polar
     * method for getting normally distributed variables from uniform
     * ones.
     *
     * G.E.P. Box and Mervin E. Muller,
     * A Note on the Generation of Random Normal Deviates,
     * The Annals of Mathematical Statistics (1958),
     * Vol. 29, No. 2 pp. 610–611
     *
     * doi:10.1214/aoms/1177706645
     *
     * as given on
     * http://en.wikipedia.org/wiki/Marsaglia_polar_method
     * on 29/11/10
     */

    if (mUseLastNum) /* We use value generated at the previous call */
    {
        mUseLastNum = false;
        return mRandNum2;
    }
    else
    {
        mUseLastNum = true;
        do
        {
            mWorkingMem1 = 2.0 * this->ranf() - 1.0;
            mWorkingMem2 = 2.0 * this->ranf() - 1.0;
            mWorkingMemW = mWorkingMem1 * mWorkingMem1 + mWorkingMem2 * mWorkingMem2;
        }
        while ( mWorkingMemW >= 1.0 );

        mWorkingMemW = sqrt((-2.0*log(mWorkingMemW))/mWorkingMemW);
        mRandNum2 = mWorkingMem2*mWorkingMemW;
        return mWorkingMem1*mWorkingMemW;
    }
}
