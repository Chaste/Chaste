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
