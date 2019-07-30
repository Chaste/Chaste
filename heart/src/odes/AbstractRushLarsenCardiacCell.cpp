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

#include "AbstractRushLarsenCardiacCell.hpp"

#include <cassert>
#include <cmath>

#include "Exception.hpp"
#include "OdeSolution.hpp"
#include "TimeStepper.hpp"

AbstractRushLarsenCardiacCell::AbstractRushLarsenCardiacCell(unsigned numberOfStateVariables,
                                                             unsigned voltageIndex,
                                                             boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver>(),
                          numberOfStateVariables,
                          voltageIndex,
                          pIntracellularStimulus)
{}

AbstractRushLarsenCardiacCell::~AbstractRushLarsenCardiacCell()
{}

OdeSolution AbstractRushLarsenCardiacCell::Compute(double tStart, double tEnd, double tSamp)
{
    // In this method, we iterate over timesteps, doing the following for each:
    //   - update V using a forward Euler step
    //   - do as in ComputeExceptVoltage(t) to update the remaining state variables
    //     using Rush Larsen method or forward Euler as appropriate

    // Check length of time interval
    if (tSamp < mDt)
    {
        tSamp = mDt;
    }
    const unsigned n_steps = (unsigned) floor((tEnd - tStart)/tSamp + 0.5);
    assert(fabs(tStart+n_steps*tSamp - tEnd) < 1e-12);
    const unsigned n_small_steps = (unsigned) floor(tSamp/mDt+0.5);
    assert(fabs(mDt*n_small_steps - tSamp) < 1e-12);

    // Initialise solution store
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(n_steps);
    solutions.rGetSolutions().push_back(rGetStateVariables());
    solutions.rGetTimes().push_back(tStart);
    solutions.SetOdeSystemInformation(this->mpSystemInfo);

    std::vector<double> dy(mNumberOfStateVariables, 0);
    std::vector<double> alpha(mNumberOfStateVariables, 0);
    std::vector<double> beta(mNumberOfStateVariables, 0);

    // Loop over time
    for (unsigned i=0; i<n_steps; i++)
    {
        double curr_time = tStart;
        for (unsigned j=0; j<n_small_steps; j++)
        {
            curr_time = tStart + i*tSamp + j*mDt;
            EvaluateEquations(curr_time, dy, alpha, beta);
            UpdateTransmembranePotential(dy);
            ComputeOneStepExceptVoltage(dy, alpha, beta);
            VerifyStateVariables();
        }

        // Update solutions
        solutions.rGetSolutions().push_back(rGetStateVariables());
        solutions.rGetTimes().push_back(curr_time+mDt);
    }

    return solutions;
}

void AbstractRushLarsenCardiacCell::ComputeExceptVoltage(double tStart, double tEnd)
{
    mSetVoltageDerivativeToZero = true;
    TimeStepper stepper(tStart, tEnd, mDt);

    std::vector<double> dy(mNumberOfStateVariables, 0);
    std::vector<double> alpha(mNumberOfStateVariables, 0);
    std::vector<double> beta(mNumberOfStateVariables, 0);

    while (!stepper.IsTimeAtEnd())
    {
        EvaluateEquations(stepper.GetTime(), dy, alpha, beta);
        ComputeOneStepExceptVoltage(dy, alpha, beta);

#ifndef NDEBUG
        // Check gating variables are still in range
        VerifyStateVariables();
#endif // NDEBUG

        stepper.AdvanceOneTimeStep();
    }
    mSetVoltageDerivativeToZero = false;
}

void AbstractRushLarsenCardiacCell::SolveAndUpdateState(double tStart, double tEnd)
{
    TimeStepper stepper(tStart, tEnd, mDt);

    std::vector<double> dy(mNumberOfStateVariables, 0);
    std::vector<double> alpha(mNumberOfStateVariables, 0);
    std::vector<double> beta(mNumberOfStateVariables, 0);

    while (!stepper.IsTimeAtEnd())
    {
        EvaluateEquations(stepper.GetTime(), dy, alpha, beta);
        UpdateTransmembranePotential(dy);
        ComputeOneStepExceptVoltage(dy, alpha, beta);
        VerifyStateVariables();

        stepper.AdvanceOneTimeStep();
    }
}

void AbstractRushLarsenCardiacCell::UpdateTransmembranePotential(const std::vector<double> &rDY)
{
    unsigned v_index = GetVoltageIndex();
    rGetStateVariables()[v_index] += mDt*rDY[v_index];
}
