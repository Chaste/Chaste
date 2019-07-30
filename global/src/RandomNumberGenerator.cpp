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

#include "RandomNumberGenerator.hpp"
#include "Exception.hpp"

#if BOOST_VERSION < 106400 // #2893
// Forward compatibility with Boost 1.64 onwards
#include "Boost165ExponentialDistribution.hpp"
#include "Boost165GammaDistribution.hpp"
#endif

RandomNumberGenerator* RandomNumberGenerator::mpInstance = nullptr;

RandomNumberGenerator::RandomNumberGenerator()
        : mMersenneTwisterGenerator(0u),
          mGenerateUnitReal(mMersenneTwisterGenerator, boost::uniform_real<>()),
#if BOOST_VERSION < 106400 // #2585 and #2893
          mGenerateStandardNormal(mMersenneTwisterGenerator, boost::random::normal_distribution_v165<>(0.0, 1.0))
#else
          mGenerateStandardNormal(mMersenneTwisterGenerator, boost::normal_distribution<>(0.0, 1.0))
#endif
{
    assert(mpInstance == nullptr); // Ensure correct serialization
}

RandomNumberGenerator* RandomNumberGenerator::Instance()
{
    if (mpInstance == nullptr)
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
        mpInstance = nullptr;
    }
}

unsigned RandomNumberGenerator::randMod(unsigned base)
{
    assert(base > 0u);
    /*
     * The contents of this method are copied out of
     * boost/include/boost/random/uniform_smallint.hpp lines 235 - 255
     * as of v 1.48 (preserved at least as far as 1.51).
     * to make sure we get the same
     * result on earlier versions of boost.
     *
     * It was then simplified as we know '_min' is zero, '_max' is 'base-1u'
     * and all the types are unsigneds.
     */

    // equivalent to (eng() - eng.min()) % (_max - _min + 1) + _min,
    // but guarantees no overflow.
    unsigned base_range = boost::random::detail::subtract<unsigned>()((mMersenneTwisterGenerator.max)(), (mMersenneTwisterGenerator.min)());
    unsigned val = boost::random::detail::subtract<unsigned>()(mMersenneTwisterGenerator(), (mMersenneTwisterGenerator.min)());

    if (base - 1u >= base_range)
    {
        // This was in the original boost file for when '_min' is large, but here it is zero so
        // we shouldn't ever reach this.
        NEVER_REACHED;
        //return val;
    }
    else
    {
        return (val % base);
    }
}

double RandomNumberGenerator::ranf()
{
    return mGenerateUnitReal();
}

double RandomNumberGenerator::StandardNormalRandomDeviate()
{
    return mGenerateStandardNormal();
}

double RandomNumberGenerator::NormalRandomDeviate(double mean, double stdDev)
{
    return stdDev * StandardNormalRandomDeviate() + mean;
}

double RandomNumberGenerator::GammaRandomDeviate(double shape, double scale)
{
// make a gamma distribution and `merge' this distribution with our random number generator

#if BOOST_VERSION < 106400
    boost::random::gamma_distribution_v165<> gd(shape);
    boost::variate_generator<boost::mt19937&, boost::random::gamma_distribution_v165<> > var_gamma(mMersenneTwisterGenerator, gd);
#else
    boost::gamma_distribution<> gd(shape);
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > var_gamma(mMersenneTwisterGenerator, gd);
#endif
    return scale * var_gamma();
}

double RandomNumberGenerator::ExponentialRandomDeviate(double scale)
{
// make an exponential distribution and `merge' this distribution with our random number generator

#if BOOST_VERSION < 106400
    boost::random::exponential_distribution_v165<> ed(scale);
    boost::variate_generator<boost::mt19937&, boost::random::exponential_distribution_v165<> > var_exponential(mMersenneTwisterGenerator, ed);
#else
    boost::exponential_distribution<> ed(scale);
    boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> > var_exponential(mMersenneTwisterGenerator, ed);
#endif
    // return the random number
    return var_exponential();
}

void RandomNumberGenerator::Reseed(unsigned seed)
{
    mMersenneTwisterGenerator.seed(seed);

    // Because this does some Box-Muller type thing it remembers if you don't reset it - see #2633
    mGenerateStandardNormal.distribution().reset();

    // Probably don't need to do this, but it probably is good practice!
    mGenerateUnitReal.distribution().reset();
}

void RandomNumberGenerator::Shuffle(unsigned num, std::vector<unsigned>& rValues)
{
    rValues.resize(num);
    for (unsigned i = 0; i < num; i++)
    {
        rValues[i] = i;
    }

    for (unsigned end = num - 1; end > 0; end--)
    {
        // Pick a random integer from {0,..,end}
        unsigned k = RandomNumberGenerator::Instance()->randMod(end + 1);
        unsigned temp = rValues[end];
        rValues[end] = rValues[k];
        rValues[k] = temp;
    }
}
