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

#include "AbstractPhaseBasedCellCycleModel.hpp"

AbstractPhaseBasedCellCycleModel::AbstractPhaseBasedCellCycleModel()
    :   AbstractCellCycleModel(),
        mCurrentCellCyclePhase(M_PHASE),
        mG1Duration(DOUBLE_UNSET),
        mMinimumGapDuration(0.01), // an educated guess
        // Default parameter values all have units of hours.
        mStemCellG1Duration(14.0),
        mTransitCellG1Duration(2.0),
        mSDuration(5.0),               // apparently between 5-6 hours normally
        mG2Duration(4.0),              // apparently 3-4 hours normally
        mMDuration(1.0)               // taken from Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
{
}

AbstractPhaseBasedCellCycleModel::~AbstractPhaseBasedCellCycleModel()
{
}

AbstractPhaseBasedCellCycleModel::AbstractPhaseBasedCellCycleModel(const AbstractPhaseBasedCellCycleModel& rModel)
    : AbstractCellCycleModel(rModel),
      mCurrentCellCyclePhase(rModel.mCurrentCellCyclePhase),
      mG1Duration(rModel.mG1Duration),
      mMinimumGapDuration(rModel.mMinimumGapDuration),
      mStemCellG1Duration(rModel.mStemCellG1Duration),
      mTransitCellG1Duration(rModel.mTransitCellG1Duration),
      mSDuration(rModel.mSDuration),
      mG2Duration(rModel.mG2Duration),
      mMDuration(rModel.mMDuration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

void AbstractPhaseBasedCellCycleModel::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    mCurrentCellCyclePhase = M_PHASE;
}

bool AbstractPhaseBasedCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ((mCurrentCellCyclePhase != G_ZERO_PHASE) &&
            (GetAge() >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()))
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}

double AbstractPhaseBasedCellCycleModel::GetG1Duration() const
{
    return mG1Duration;
}

double AbstractPhaseBasedCellCycleModel::GetStemCellG1Duration() const
{
    return mStemCellG1Duration;
}

double AbstractPhaseBasedCellCycleModel::GetTransitCellG1Duration() const
{
    return mTransitCellG1Duration;
}

double AbstractPhaseBasedCellCycleModel::GetSG2MDuration() const
{
    return mSDuration + mG2Duration + mMDuration;
}

double AbstractPhaseBasedCellCycleModel::GetSDuration() const
{
    return mSDuration;
}

double AbstractPhaseBasedCellCycleModel::GetG2Duration() const
{
    return mG2Duration;
}

double AbstractPhaseBasedCellCycleModel::GetMDuration() const
{
    return mMDuration;
}

CellCyclePhase AbstractPhaseBasedCellCycleModel::GetCurrentCellCyclePhase() const
{
    return mCurrentCellCyclePhase;
}

double AbstractPhaseBasedCellCycleModel::GetAverageTransitCellCycleTime()
{
    return mTransitCellG1Duration + GetSG2MDuration();
}

double AbstractPhaseBasedCellCycleModel::GetAverageStemCellCycleTime()
{
    return mStemCellG1Duration + GetSG2MDuration();
}

double AbstractPhaseBasedCellCycleModel::GetMinimumGapDuration() const
{
    return mMinimumGapDuration;
}

void AbstractPhaseBasedCellCycleModel::SetStemCellG1Duration(double stemCellG1Duration)
{
    assert(stemCellG1Duration >= 0.0);
    mStemCellG1Duration = stemCellG1Duration;
}

void AbstractPhaseBasedCellCycleModel::SetTransitCellG1Duration(double transitCellG1Duration)
{
    assert(transitCellG1Duration >= 0.0);
    mTransitCellG1Duration = transitCellG1Duration;
}

void AbstractPhaseBasedCellCycleModel::SetSDuration(double SDuration)
{
    assert(SDuration >= 0.0);
    mSDuration = SDuration;
}

void AbstractPhaseBasedCellCycleModel::SetG2Duration(double G2Duration)
{
    assert(G2Duration >= 0.0);
    mG2Duration = G2Duration;
}

void AbstractPhaseBasedCellCycleModel::SetMDuration(double MDuration)
{
    assert(MDuration >= 0.0);
    mMDuration = MDuration;
}

void AbstractPhaseBasedCellCycleModel::SetMinimumGapDuration(double minimumGapDuration)
{
    assert(minimumGapDuration > 0.0);
    mMinimumGapDuration = minimumGapDuration;
}

void AbstractPhaseBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<StemCellG1Duration>" << mStemCellG1Duration << "</StemCellG1Duration>\n";
    *rParamsFile << "\t\t\t<TransitCellG1Duration>" << mTransitCellG1Duration << "</TransitCellG1Duration>\n";
    *rParamsFile << "\t\t\t<SDuration>" << mSDuration << "</SDuration>\n";
    *rParamsFile << "\t\t\t<G2Duration>" << mG2Duration << "</G2Duration>\n";
    *rParamsFile << "\t\t\t<MDuration>" << mMDuration << "</MDuration>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
