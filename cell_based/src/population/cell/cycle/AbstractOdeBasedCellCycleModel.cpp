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

#include "AbstractOdeBasedCellCycleModel.hpp"
#include <iostream>
#include <cassert>
#include "Exception.hpp"

AbstractOdeBasedCellCycleModel::AbstractOdeBasedCellCycleModel(double lastTime,
                                                               boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : CellCycleModelOdeHandler(lastTime, pOdeSolver),
      mDivideTime(lastTime),
      mFinishedRunningOdes(false),
      mG2PhaseStartTime(DBL_MAX)
{
    AbstractCellCycleModel::SetBirthTime(lastTime);
}

AbstractOdeBasedCellCycleModel::~AbstractOdeBasedCellCycleModel()
{
}

void AbstractOdeBasedCellCycleModel::SetBirthTime(double birthTime)
{
    AbstractCellCycleModel::SetBirthTime(birthTime);
    mLastTime = birthTime;
    mDivideTime = birthTime;
}

std::vector<double> AbstractOdeBasedCellCycleModel::GetProteinConcentrations() const
{
    assert(mpOdeSystem != NULL);
    return mpOdeSystem->rGetStateVariables();
}

void AbstractOdeBasedCellCycleModel::SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations)
{
    assert(mpOdeSystem != NULL);
    assert(proteinConcentrations.size()==mpOdeSystem->rGetStateVariables().size());
    mLastTime = lastTime;
    mpOdeSystem->SetStateVariables(proteinConcentrations);
}

void AbstractOdeBasedCellCycleModel::UpdateCellCyclePhase()
{
    assert(mpOdeSystem != NULL);

    double current_time = SimulationTime::Instance()->GetTime();

    // Update the phase from M to G1 when necessary
    if (mCurrentCellCyclePhase == M_PHASE)
    {
        double m_duration = GetMDuration();
        if (GetAge() >= m_duration)
        {
            mCurrentCellCyclePhase = G_ONE_PHASE;
            mLastTime = m_duration + mBirthTime;
        }
        else
        {
            // Still dividing; don't run ODEs
            return;
        }
    }

    if (current_time > mLastTime)
    {
        if (!mFinishedRunningOdes)
        {
            // Update whether a stopping event has occurred
            mFinishedRunningOdes = SolveOdeToTime(current_time);

            // Check no concentrations have gone negative
            for (unsigned i=0; i<mpOdeSystem->GetNumberOfStateVariables(); i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i] < -DBL_EPSILON)
                {
                    #define COVERAGE_IGNORE
                    EXCEPTION("A protein concentration " << i << " has gone negative (" <<
                              mpOdeSystem->rGetStateVariables()[i] << ")\n"
                              << "Chaste predicts that the CellCycleModel numerical method is probably unstable.");
                    #undef COVERAGE_IGNORE
                }
            }

            if (mFinishedRunningOdes)
            {
                // Update durations of each phase
                mG1Duration = GetOdeStopTime() - mBirthTime - GetMDuration();
                mG2PhaseStartTime = GetOdeStopTime() + GetSDuration();
                mDivideTime = mG2PhaseStartTime + GetG2Duration();

                // Update phase
                if (current_time >= mG2PhaseStartTime)
                {
                    mCurrentCellCyclePhase = G_TWO_PHASE;
                }
                else
                {
                    mCurrentCellCyclePhase = S_PHASE;
                }
            }
        }
        else
        {
            // ODE model finished, just increasing time until division...
            if (current_time >= mG2PhaseStartTime)
            {
                mCurrentCellCyclePhase = G_TWO_PHASE;
            }
        }
    }
}

void AbstractOdeBasedCellCycleModel::ResetForDivision()
{
    assert(mFinishedRunningOdes);
    AbstractCellCycleModel::ResetForDivision();
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
    mFinishedRunningOdes = false;
    mG1Duration = DBL_MAX;
    mDivideTime = DBL_MAX;
}

double AbstractOdeBasedCellCycleModel::GetOdeStopTime()
{
    double stop_time = DOUBLE_UNSET;
    if (mpOdeSolver->StoppingEventOccurred())
    {
        stop_time = mpOdeSolver->GetStoppingTime();
    }
    return stop_time;
}

void AbstractOdeBasedCellCycleModel::SetFinishedRunningOdes(bool finishedRunningOdes)
{
    mFinishedRunningOdes = finishedRunningOdes;
}

void AbstractOdeBasedCellCycleModel::SetDivideTime(double divideTime)
{
    mDivideTime = divideTime;
}

void AbstractOdeBasedCellCycleModel::SetG2PhaseStartTime(double g2PhaseStartTime)
{
    mG2PhaseStartTime = g2PhaseStartTime;
}

void AbstractOdeBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
