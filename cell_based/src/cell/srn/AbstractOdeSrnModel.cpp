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

#include "AbstractOdeSrnModel.hpp"

AbstractOdeSrnModel::AbstractOdeSrnModel(unsigned stateSize, boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractSrnModel(),
      CellCycleModelOdeHandler(SimulationTime::Instance()->GetTime(), pOdeSolver),
      mStateSize(stateSize)
{
}

AbstractOdeSrnModel::AbstractOdeSrnModel(const AbstractOdeSrnModel& rModel)
    : AbstractSrnModel(rModel),
      CellCycleModelOdeHandler(rModel),
      mInitialConditions(rModel.mInitialConditions),
      mStateSize(rModel.mStateSize)
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

AbstractOdeSrnModel::~AbstractOdeSrnModel()
{
}

void AbstractOdeSrnModel::SimulateToCurrentTime()
{
    assert(mpOdeSystem != nullptr);
    assert(SimulationTime::Instance()->IsStartTimeSetUp());

    double current_time = SimulationTime::Instance()->GetTime();

    // Run ODEs if needed
    if (current_time > mLastTime)
    {
        if (!this->mFinishedRunningOdes)
        {
            // Update whether a stopping event has occurred
            this->mFinishedRunningOdes = SolveOdeToTime(current_time);
        }
        else
        {
            // ODE model finished, just increasing time...
        }
    }

    // Update the SimulatedToTime value
    mLastTime = current_time;
    SetSimulatedToTime(current_time);
}

void AbstractOdeSrnModel::Initialise(AbstractOdeSystem* pOdeSystem)
{
    assert(mpOdeSystem == nullptr);
    assert(mpCell != nullptr);

    mpOdeSystem = pOdeSystem;
    if (mInitialConditions == std::vector<double>())
    {
        mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    }
    else
    {
        mpOdeSystem->SetStateVariables(mInitialConditions);
    }

    SetLastTime(mSimulatedToTime);
}

void AbstractOdeSrnModel::ResetForDivision()
{
    AbstractSrnModel::ResetForDivision();
    assert(mLastTime == mSimulatedToTime);
    mFinishedRunningOdes = false;
}

void AbstractOdeSrnModel::SetInitialConditions(std::vector<double> initialConditions)
{
    assert(initialConditions.size() == mStateSize);
    mInitialConditions = initialConditions;
}

void AbstractOdeSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AbstractOdeSrnModel)
