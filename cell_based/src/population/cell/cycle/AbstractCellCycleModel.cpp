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

void AbstractCellCycleModel::SetCellProliferativeType(CellProliferativeType cellType)
{
    mCellProliferativeType = cellType;
}

CellProliferativeType AbstractCellCycleModel::GetCellProliferativeType() const
{
    return mCellProliferativeType;
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
