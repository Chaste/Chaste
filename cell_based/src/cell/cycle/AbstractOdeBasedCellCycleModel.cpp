/*

Copyright (c) 2005-2016, University of Oxford.
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
      mFinishedRunningOdes(false)
{
    AbstractCellCycleModel::SetBirthTime(lastTime);
}

AbstractOdeBasedCellCycleModel::~AbstractOdeBasedCellCycleModel()
{
}

AbstractOdeBasedCellCycleModel::AbstractOdeBasedCellCycleModel(const AbstractOdeBasedCellCycleModel& rModel)
    : AbstractCellCycleModel(rModel),
      CellCycleModelOdeHandler(rModel),
      mDivideTime(rModel.mDivideTime),
      mFinishedRunningOdes(false)
{
    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables will already
     * have been correctly initialized in its constructor or parent classes.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     *
     */
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

bool AbstractOdeBasedCellCycleModel::ReadyToDivide()
{
    assert(mpCell != NULL);

    if (!mReadyToDivide)
    {
        assert(mpOdeSystem != NULL);

        double current_time = SimulationTime::Instance()->GetTime();

        if (current_time > mLastTime)
        {
            // Update whether a stopping event has occurred
            mReadyToDivide = SolveOdeToTime(current_time);

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
            if (mReadyToDivide)
            {
                mDivideTime = GetOdeStopTime();
            }
        }
    }
    return mReadyToDivide;
}

void AbstractOdeBasedCellCycleModel::ResetForDivision()
{
    assert(mReadyToDivide);
    AbstractCellCycleModel::ResetForDivision();
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
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

void AbstractOdeBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
