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

#include "StochasticOxygenBasedCellCycleModel.hpp"

StochasticOxygenBasedCellCycleModel::StochasticOxygenBasedCellCycleModel()
    : SimpleOxygenBasedCellCycleModel()
{
}

StochasticOxygenBasedCellCycleModel::StochasticOxygenBasedCellCycleModel(const StochasticOxygenBasedCellCycleModel& rModel)
   : SimpleOxygenBasedCellCycleModel(rModel),
     mStochasticG2Duration(rModel.mStochasticG2Duration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentHypoxicDuration, mCurrentHypoxiaOnsetTime,
     * mHypoxicConcentration, mQuiescentConcentration and mCriticalHypoxicDuration
     * are initialized in the SimpleOxygenBasedCellCycleModel constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration and mStochasticG2Duration are (re)set as
     * soon as InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

void StochasticOxygenBasedCellCycleModel::GenerateStochasticG2Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double mean = AbstractPhaseBasedCellCycleModel::GetG2Duration();
    double standard_deviation = 1.0;

    mStochasticG2Duration = p_gen->NormalRandomDeviate(mean, standard_deviation);

    // Check that the normal random deviate has not returned a small or negative G2 duration
    if (mStochasticG2Duration < mMinimumGapDuration)
    {
        mStochasticG2Duration = mMinimumGapDuration;
    }
}

void StochasticOxygenBasedCellCycleModel::InitialiseDaughterCell()
{
    SimpleOxygenBasedCellCycleModel::InitialiseDaughterCell();
    GenerateStochasticG2Duration();
}

void StochasticOxygenBasedCellCycleModel::Initialise()
{
    AbstractSimplePhaseBasedCellCycleModel::Initialise();
    GenerateStochasticG2Duration();
}

void StochasticOxygenBasedCellCycleModel::ResetForDivision()
{
    SimpleOxygenBasedCellCycleModel::ResetForDivision();
    GenerateStochasticG2Duration();
}

double StochasticOxygenBasedCellCycleModel::GetG2Duration() const
{
    return mStochasticG2Duration;
}

AbstractCellCycleModel* StochasticOxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new StochasticOxygenBasedCellCycleModel(*this);
}

void StochasticOxygenBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    SimpleOxygenBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticOxygenBasedCellCycleModel)
