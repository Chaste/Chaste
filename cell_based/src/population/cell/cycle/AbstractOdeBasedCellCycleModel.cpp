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
