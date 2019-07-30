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

#include "AbstractSrnModel.hpp"

AbstractSrnModel::AbstractSrnModel()
    : mSimulatedToTime(SimulationTime::Instance()->GetTime())
{
}

AbstractSrnModel::~AbstractSrnModel()
{
}

void AbstractSrnModel::Initialise()
{
}

void AbstractSrnModel::InitialiseDaughterCell()
{
}

AbstractSrnModel::AbstractSrnModel(const AbstractSrnModel& rModel)
    : mSimulatedToTime(rModel.GetSimulatedToTime())
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
}

void AbstractSrnModel::ResetForDivision()
{
    // Make sure we're at the current time
    SimulateToCurrentTime();
    assert(mSimulatedToTime == SimulationTime::Instance()->GetTime());
}

void AbstractSrnModel::SetCell(CellPtr pCell)
{
    mpCell = pCell;
}

CellPtr AbstractSrnModel::GetCell()
{
    assert(mpCell != nullptr);
    return mpCell;
}

void AbstractSrnModel::SetSimulatedToTime(double simulatedToTime)
{
    mSimulatedToTime = simulatedToTime;
}

double AbstractSrnModel::GetSimulatedToTime() const
{
    return mSimulatedToTime;
}

/*
///\todo - make abstract? (#2752)
void AbstractSrnModel::SimulateToCurrentTime()
{
    // This should be overridden
    SetSimulatedToTime(SimulationTime::Instance()->GetTime());
}
*/

void AbstractSrnModel::OutputSrnModelInfo(out_stream& rParamsFile)
{
    std::string srn_model_type = GetIdentifier();

    *rParamsFile << "\t\t<" << srn_model_type << ">\n";
    OutputSrnModelParameters(rParamsFile);
    *rParamsFile << "\t\t</" << srn_model_type << ">\n";
}

void AbstractSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
}
