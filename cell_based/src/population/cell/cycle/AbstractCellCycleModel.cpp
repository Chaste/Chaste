/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "AbstractCellCycleModel.hpp"

AbstractCellCycleModel::AbstractCellCycleModel()
    : mBirthTime(SimulationTime::Instance()->GetTime()),
      mCurrentCellCyclePhase(M_PHASE),
      mG1Duration(DOUBLE_UNSET),
      mReadyToDivide(false),
      mDimension(UNSIGNED_UNSET),
      mMinimumGapDuration(0.01), // an educated guess
      // Default parameter values all have units of hours.
      mStemCellG1Duration(14.0),
      mTransitCellG1Duration(2.0),
      mSDuration(5.0),               // apparently between 5-6 hours normally
      mG2Duration(4.0),              // apparently 3-4 hours normally
      mMDuration(1.0)               // taken from Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
{
}

AbstractCellCycleModel::~AbstractCellCycleModel()
{
}

void AbstractCellCycleModel::Initialise()
{
}

void AbstractCellCycleModel::InitialiseDaughterCell()
{
}

void AbstractCellCycleModel::SetCell(CellPtr pCell)
{
    mpCell = pCell;
}

CellPtr AbstractCellCycleModel::GetCell()
{
    assert(mpCell != NULL);
    return mpCell;
}

void AbstractCellCycleModel::SetBirthTime(double birthTime)
{
    mBirthTime = birthTime;
}

double AbstractCellCycleModel::GetBirthTime() const
{
    return mBirthTime;
}

double AbstractCellCycleModel::GetAge()
{
    return SimulationTime::Instance()->GetTime() - mBirthTime;
}

CellCyclePhase AbstractCellCycleModel::GetCurrentCellCyclePhase()
{
    return mCurrentCellCyclePhase;
}

void AbstractCellCycleModel::ResetForDivision()
{
    assert(mReadyToDivide);
    mCurrentCellCyclePhase = M_PHASE;
    mReadyToDivide = false;
}

double AbstractCellCycleModel::GetG1Duration()
{
    return mG1Duration;
}

///////////////////////////////////////////////////////////////////////
// Getter methods
///////////////////////////////////////////////////////////////////////

double AbstractCellCycleModel::GetStemCellG1Duration()
{
    return mStemCellG1Duration;
}

double AbstractCellCycleModel::GetTransitCellG1Duration()
{
    return mTransitCellG1Duration;
}

double AbstractCellCycleModel::GetSG2MDuration()
{
    return mSDuration + mG2Duration + mMDuration;
}

double AbstractCellCycleModel::GetSDuration()
{
    return mSDuration;
}

double AbstractCellCycleModel::GetG2Duration()
{
    return mG2Duration;
}

double AbstractCellCycleModel::GetMDuration()
{
    return mMDuration;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

void AbstractCellCycleModel::SetStemCellG1Duration(double stemCellG1Duration)
{
    assert(stemCellG1Duration > 0.0);
    mStemCellG1Duration = stemCellG1Duration;
}
void AbstractCellCycleModel::SetTransitCellG1Duration(double transitCellG1Duration)
{
    assert(transitCellG1Duration > 0.0);
    mTransitCellG1Duration = transitCellG1Duration;
}
void AbstractCellCycleModel::SetSDuration(double SDuration)
{
    assert(SDuration > 0.0);
    mSDuration = SDuration;
}
void AbstractCellCycleModel::SetG2Duration(double G2Duration)
{
    assert(G2Duration > 0.0);
    mG2Duration = G2Duration;
}
void AbstractCellCycleModel::SetMDuration(double MDuration)
{
    assert(MDuration > 0.0);
    mMDuration = MDuration;
}

bool AbstractCellCycleModel::ReadyToDivide()
{
    assert(mpCell != NULL);

    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ( (mCurrentCellCyclePhase != G_ZERO_PHASE) &&
             (GetAge() >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()) )
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}

void AbstractCellCycleModel::SetDimension(unsigned dimension)
{
    if (dimension != 1 && dimension !=2 && dimension != 3)
    {
        EXCEPTION("Dimension must be 1, 2 or 3");
    }
    mDimension = dimension;
}

unsigned AbstractCellCycleModel::GetDimension()
{
    return mDimension;
}

double AbstractCellCycleModel::GetAverageTransitCellCycleTime()
{
    return mTransitCellG1Duration + GetSG2MDuration();
}

double AbstractCellCycleModel::GetAverageStemCellCycleTime()
{
    return mStemCellG1Duration + GetSG2MDuration();
}

bool AbstractCellCycleModel::CanCellTerminallyDifferentiate()
{
    return true;
}

void AbstractCellCycleModel::SetMinimumGapDuration(double minimumGapDuration)
{
    assert(minimumGapDuration > 0.0);
    mMinimumGapDuration = minimumGapDuration;
}

double AbstractCellCycleModel::GetMinimumGapDuration()
{
    return mMinimumGapDuration;
}

void AbstractCellCycleModel::OutputCellCycleModelInfo(out_stream& rParamsFile)
{
    std::string cell_cycle_model_type = GetIdentifier();

    *rParamsFile << "\t\t<" << cell_cycle_model_type << ">\n";
    OutputCellCycleModelParameters(rParamsFile);
    *rParamsFile << "\t\t</" << cell_cycle_model_type << ">\n";
}

void AbstractCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<StemCellG1Duration>" << mStemCellG1Duration << "</StemCellG1Duration>\n";
    *rParamsFile << "\t\t\t<TransitCellG1Duration>" << mTransitCellG1Duration << "</TransitCellG1Duration>\n";
    *rParamsFile << "\t\t\t<SDuration>" << mSDuration << "</SDuration>\n";
    *rParamsFile << "\t\t\t<G2Duration>" << mG2Duration << "</G2Duration>\n";
    *rParamsFile << "\t\t\t<MDuration>" << mMDuration << "</MDuration>\n";
}
