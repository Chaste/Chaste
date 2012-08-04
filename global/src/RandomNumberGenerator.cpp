/*

Copyright (c) 2005-2012, University of Oxford.
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
