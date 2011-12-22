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

#include "AbstractSimpleCellCycleModel.hpp"
#include "Exception.hpp"

AbstractSimpleCellCycleModel::AbstractSimpleCellCycleModel()
{
}

AbstractSimpleCellCycleModel::~AbstractSimpleCellCycleModel()
{
}

void AbstractSimpleCellCycleModel::Initialise()
{
    SetG1Duration();
}

void AbstractSimpleCellCycleModel::InitialiseDaughterCell()
{
    AbstractCellCycleModel::InitialiseDaughterCell();
    SetG1Duration();
}

void AbstractSimpleCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    switch (mCellProliferativeType)
    {
        case STEM:
            mG1Duration = GetStemCellG1Duration();
            break;
        case TRANSIT:
            mG1Duration = GetTransitCellG1Duration();
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }
}

void AbstractSimpleCellCycleModel::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    mBirthTime = SimulationTime::Instance()->GetTime();
    SetG1Duration();
}

void AbstractSimpleCellCycleModel::UpdateCellCyclePhase()
{
    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    if (mCellProliferativeType == DIFFERENTIATED)
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (time_since_birth < GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
}

void AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
