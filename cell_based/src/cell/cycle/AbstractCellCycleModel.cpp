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

#include "AbstractCellCycleModel.hpp"

AbstractCellCycleModel::AbstractCellCycleModel()
    : mBirthTime(SimulationTime::Instance()->GetTime()),
      mReadyToDivide(false),
      mDimension(UNSIGNED_UNSET)
{
}

AbstractCellCycleModel::~AbstractCellCycleModel()
{
}

AbstractCellCycleModel::AbstractCellCycleModel(const AbstractCellCycleModel& rModel)
    : mBirthTime(rModel.mBirthTime),
      mReadyToDivide(rModel.mReadyToDivide),
      mDimension(rModel.mDimension)
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
    assert(mpCell != nullptr);
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

void AbstractCellCycleModel::ResetForDivision()
{
    assert(mReadyToDivide);
    mReadyToDivide = false;
    mBirthTime = SimulationTime::Instance()->GetTime();
}

void AbstractCellCycleModel::SetDimension(unsigned dimension)
{
    if (dimension != 1 && dimension !=2 && dimension != 3 && dimension != UNSIGNED_UNSET)
    {
        EXCEPTION("Dimension must be 1, 2, 3 or UNSIGNED_UNSET");
    }
    mDimension = dimension;
}

unsigned AbstractCellCycleModel::GetDimension() const
{
    return mDimension;
}

bool AbstractCellCycleModel::CanCellTerminallyDifferentiate()
{
    return true;
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
 // None to output
}
