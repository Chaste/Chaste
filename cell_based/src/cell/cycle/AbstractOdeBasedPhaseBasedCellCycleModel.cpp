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

#include "AbstractOdeBasedPhaseBasedCellCycleModel.hpp"

AbstractOdeBasedPhaseBasedCellCycleModel::AbstractOdeBasedPhaseBasedCellCycleModel(double lastTime,
                                                               boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : CellCycleModelOdeHandler(lastTime, pOdeSolver),
      mDivideTime(lastTime),
      mG2PhaseStartTime(DBL_MAX)
{
    AbstractPhaseBasedCellCycleModel::SetBirthTime(lastTime);
}

AbstractOdeBasedPhaseBasedCellCycleModel::~AbstractOdeBasedPhaseBasedCellCycleModel()
{
}

AbstractOdeBasedPhaseBasedCellCycleModel::AbstractOdeBasedPhaseBasedCellCycleModel(const AbstractOdeBasedPhaseBasedCellCycleModel& rModel)
    : AbstractPhaseBasedCellCycleModel(rModel),
      CellCycleModelOdeHandler(rModel),
      mDivideTime(rModel.mDivideTime),
      mG2PhaseStartTime(rModel.mG2PhaseStartTime)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

void AbstractOdeBasedPhaseBasedCellCycleModel::SetBirthTime(double birthTime)
{
    AbstractPhaseBasedCellCycleModel::SetBirthTime(birthTime);
    mLastTime = birthTime;
    mDivideTime = birthTime;
}

void AbstractOdeBasedPhaseBasedCellCycleModel::UpdateCellCyclePhase()
{
    assert(mpOdeSystem != nullptr);

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
        if (!this->mFinishedRunningOdes)
        {
            // Update whether a stopping event has occurred
            this->mFinishedRunningOdes = SolveOdeToTime(current_time);

            // Check no concentrations have gone negative
            for (unsigned i=0; i<mpOdeSystem->GetNumberOfStateVariables(); i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i] < -DBL_EPSILON)
                {
                    // LCOV_EXCL_START
                    EXCEPTION("A protein concentration " << i << " has gone negative (" <<
                              mpOdeSystem->rGetStateVariables()[i] << ")\n"
                              << "Chaste predicts that the CellCycleModel numerical method is probably unstable.");
                    // LCOV_EXCL_STOP
                }
            }

            if (this->mFinishedRunningOdes)
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

void AbstractOdeBasedPhaseBasedCellCycleModel::ResetForDivision()
{
    assert(this->mFinishedRunningOdes);
    AbstractPhaseBasedCellCycleModel::ResetForDivision();
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
    this->mFinishedRunningOdes = false;
    mG1Duration = DBL_MAX;
    mDivideTime = DBL_MAX;
}

double AbstractOdeBasedPhaseBasedCellCycleModel::GetOdeStopTime()
{
    double stop_time = DOUBLE_UNSET;
    if (mpOdeSolver->StoppingEventOccurred())
    {
        stop_time = mpOdeSolver->GetStoppingTime();
    }
    return stop_time;
}

void AbstractOdeBasedPhaseBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractPhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
