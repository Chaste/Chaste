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

#include "Exception.hpp"
#include "CellCycleTimesGenerator.hpp"

CellCycleTimesGenerator* CellCycleTimesGenerator::mpInstance = NULL;

CellCycleTimesGenerator::CellCycleTimesGenerator()
    : mRandomSeed(0u),
      mCurrentIndex(0u),
      mRate(1.0/2.0),
      mVectorCreated(false)
{
    assert(mpInstance == NULL); // Ensure correct serialization
}

CellCycleTimesGenerator* CellCycleTimesGenerator::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CellCycleTimesGenerator();
    }
    return mpInstance;
}

void CellCycleTimesGenerator::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

void CellCycleTimesGenerator::SetRandomSeed(unsigned randomSeed)
{
    mRandomSeed = randomSeed;
}

unsigned CellCycleTimesGenerator::GetRandomSeed()
{
    return mRandomSeed;
}

void CellCycleTimesGenerator::GenerateCellCycleTimeSequence()
{
    if (mVectorCreated)
    {
        EXCEPTION("Trying to generate the cell cycle times twice. Need to call CellCycleTimesGenerator::Destroy() first.");
    }
    else
    {
        unsigned number_stored_times = 15000u;
        mCellCycleTimes.reserve(number_stored_times);

        RandomNumberGenerator* p_random_number_generator = RandomNumberGenerator::Instance();
        p_random_number_generator->Reseed(mRandomSeed);

        for (unsigned index = 0; index < 15000u; index++)
        {
            mCellCycleTimes.push_back( p_random_number_generator->ExponentialRandomDeviate(mRate) );
        }

        mCurrentIndex = 0;
        mVectorCreated = true;
    }
}

void CellCycleTimesGenerator::SetRate(double rate)
{
    if( mVectorCreated)
    {
        EXCEPTION("You cannot reset the rate after cell cycle times are created.");
    }
    else
    {
        mRate = rate;
    }
}

double CellCycleTimesGenerator::GetRate()
{
    return mRate;
}

double CellCycleTimesGenerator::GetNextCellCycleTime()
{
    if(!mVectorCreated)
    {
        EXCEPTION("When using FixedSequenceCellCycleModel one must call CellCycleTimesGenerator::Instance()->GenerateCellCycleTimeSequence()"
                " before the start of the simulation.");
    }
    else
    {
        double new_cell_cycle_time = mCellCycleTimes[mCurrentIndex];
        mCurrentIndex += 1u;
        return new_cell_cycle_time;
    }
}

