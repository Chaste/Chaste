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

#include "CellCycleModelOdeHandler.hpp"

CellCycleModelOdeHandler::CellCycleModelOdeHandler(double lastTime,
                                                   boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : mDt(DOUBLE_UNSET),
      mpOdeSystem(nullptr),
      mpOdeSolver(pOdeSolver),
      mLastTime(lastTime),
      mFinishedRunningOdes(false)
{
}

CellCycleModelOdeHandler::~CellCycleModelOdeHandler()
{
    if (mpOdeSystem != nullptr)
    {
        delete mpOdeSystem;
    }
}

CellCycleModelOdeHandler::CellCycleModelOdeHandler(const CellCycleModelOdeHandler& rHandler)
    : mDt(rHandler.mDt),
      mpOdeSystem(rHandler.mpOdeSystem),
      mpOdeSolver(rHandler.mpOdeSolver),
      mLastTime(rHandler.mLastTime),
      mFinishedRunningOdes(rHandler.mFinishedRunningOdes)
{
}

void CellCycleModelOdeHandler::SetOdeSystem(AbstractOdeSystem* pOdeSystem)
{
    mpOdeSystem = pOdeSystem;
}

AbstractOdeSystem* CellCycleModelOdeHandler::GetOdeSystem() const
{
    return mpOdeSystem;
}

const boost::shared_ptr<AbstractCellCycleModelOdeSolver> CellCycleModelOdeHandler::GetOdeSolver() const
{
    return mpOdeSolver;
}

void CellCycleModelOdeHandler::SetDt(double timeStep)
{
    mDt = timeStep;
}

double CellCycleModelOdeHandler::GetDt()
{
    if (mDt == DOUBLE_UNSET)
    {
        if (mpOdeSolver->IsAdaptive())
        {
            mDt = SimulationTime::Instance()->GetTimeStep();
        }
        else
        {
            mDt = 0.0001; // Some models need this, so let's pick a safe default
        }
    }
    return mDt;
}

bool CellCycleModelOdeHandler::SolveOdeToTime(double currentTime)
{
    bool stopping_event_occurred = false;
    if (mLastTime < currentTime)
    {
        AdjustOdeParameters(currentTime);

        mpOdeSolver->SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, GetDt());

        stopping_event_occurred = mpOdeSolver->StoppingEventOccurred();
        if (stopping_event_occurred)
        {
            mLastTime = mpOdeSolver->GetStoppingTime();
        }
        else
        {
            mLastTime = currentTime;
        }
    }
    return stopping_event_occurred;
}

void CellCycleModelOdeHandler::AdjustOdeParameters(double currentTime)
{
}

void CellCycleModelOdeHandler::SetLastTime(double lastTime)
{
    mLastTime = lastTime;
}

void CellCycleModelOdeHandler::SetStateVariables(const std::vector<double>& rStateVariables)
{
    assert(mpOdeSystem);
    mpOdeSystem->SetStateVariables(rStateVariables);
}

std::vector<double> CellCycleModelOdeHandler::GetProteinConcentrations() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->rGetStateVariables();
}

void CellCycleModelOdeHandler::SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations)
{
    assert(mpOdeSystem != nullptr);
    assert(proteinConcentrations.size()==mpOdeSystem->rGetStateVariables().size());
    mLastTime = lastTime;
    mpOdeSystem->SetStateVariables(proteinConcentrations);
}
