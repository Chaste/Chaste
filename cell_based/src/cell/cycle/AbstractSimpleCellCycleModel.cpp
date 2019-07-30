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

#include "AbstractSimpleCellCycleModel.hpp"

AbstractSimpleCellCycleModel::AbstractSimpleCellCycleModel()
    : AbstractCellCycleModel(),
      mCellCycleDuration(DOUBLE_UNSET)
{
}

AbstractSimpleCellCycleModel::~AbstractSimpleCellCycleModel()
{
}

AbstractSimpleCellCycleModel::AbstractSimpleCellCycleModel(const AbstractSimpleCellCycleModel& rModel)
    : AbstractCellCycleModel(rModel),
      mCellCycleDuration(rModel.mCellCycleDuration)
{
    /*
     * Initialize only the member variable defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

void AbstractSimpleCellCycleModel::Initialise()
{
    SetCellCycleDuration();
}

void AbstractSimpleCellCycleModel::InitialiseDaughterCell()
{
    AbstractCellCycleModel::InitialiseDaughterCell();
    SetCellCycleDuration();
}

bool AbstractSimpleCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        if (GetAge() >= mCellCycleDuration )
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}

void AbstractSimpleCellCycleModel::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    SetCellCycleDuration();
}

double AbstractSimpleCellCycleModel::GetCellCycleDuration() const
{
    return mCellCycleDuration;
}

void AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
