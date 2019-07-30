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


#include "BackwardEulerIvpOdeSolver.hpp"
#include <cmath>

void BackwardEulerIvpOdeSolver::ComputeResidual(AbstractOdeSystem* pAbstractOdeSystem,
                                                double timeStep,
                                                double time,
                                                std::vector<double>& rCurrentYValues,
                                                std::vector<double>& rCurrentGuess)
{
    std::vector<double> dy(mSizeOfOdeSystem); //For JC to optimize
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, rCurrentGuess, dy);
    for (unsigned i=0; i<mSizeOfOdeSystem; i++)
    {
        mResidual[i] = rCurrentGuess[i] - timeStep * dy[i] - rCurrentYValues[i];
    }
}

void BackwardEulerIvpOdeSolver::ComputeJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                                                double timeStep,
                                                double time,
                                                std::vector<double>& rCurrentYValues,
                                                std::vector<double>& rCurrentGuess)
{
    for (unsigned i=0; i<mSizeOfOdeSystem; i++)
    {
        for (unsigned j=0; j<mSizeOfOdeSystem; j++)
        {
            mJacobian[i][j] = 0.0;
        }
    }

    if (pAbstractOdeSystem->GetUseAnalyticJacobian() && !mForceUseOfNumericalJacobian)
    {
        // The ODE system has an analytic jacobian, so use that
        AbstractOdeSystemWithAnalyticJacobian* p_ode_system
            = static_cast<AbstractOdeSystemWithAnalyticJacobian*>(pAbstractOdeSystem);
        p_ode_system->AnalyticJacobian(rCurrentGuess, mJacobian, time, timeStep);
    }
    else
    {
        ComputeNumericalJacobian(pAbstractOdeSystem,
                                 timeStep,
                                 time,
                                 rCurrentYValues,
                                 rCurrentGuess);
    }
}

void BackwardEulerIvpOdeSolver::SolveLinearSystem()
{
    double fact;
    for (unsigned i=0; i<mSizeOfOdeSystem; i++)
    {
        for (unsigned ii=i+1; ii<mSizeOfOdeSystem; ii++)
        {
            fact = mJacobian[ii][i]/mJacobian[i][i];
            for (unsigned j=i; j<mSizeOfOdeSystem; j++)
            {
                mJacobian[ii][j] -= fact*mJacobian[i][j];
            }
            mResidual[ii] -= fact*mResidual[i];
        }
    }
    // This needs to int, since a downloop in unsigned won't terminate properly
    for (int i=mSizeOfOdeSystem-1; i>=0; i--)
    {
        mUpdate[i] = mResidual[i];
        for (unsigned j=i+1; j<mSizeOfOdeSystem; j++)
        {
            mUpdate[i] -= mJacobian[i][j]*mUpdate[j];
        }
        mUpdate[i] /= mJacobian[i][i];
    }
}

double BackwardEulerIvpOdeSolver::ComputeNorm(double* pVector)
{
    double norm = 0.0;
    for (unsigned i=0; i<mSizeOfOdeSystem; i++)
    {
        if (fabs(pVector[i]) > norm)
        {
            norm = fabs(pVector[i]);
        }
    }
    return norm;
}

void BackwardEulerIvpOdeSolver::ComputeNumericalJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                                                         double timeStep,
                                                         double time,
                                                         std::vector<double>& rCurrentYValues,
                                                         std::vector<double>& rCurrentGuess)
{
    std::vector<double> residual(mSizeOfOdeSystem);
    std::vector<double> residual_perturbed(mSizeOfOdeSystem);
    std::vector<double> guess_perturbed(mSizeOfOdeSystem);

    double epsilon = mNumericalJacobianEpsilon;

    ComputeResidual(pAbstractOdeSystem, timeStep, time, rCurrentYValues, rCurrentGuess);
    for (unsigned i=0; i<mSizeOfOdeSystem; i++)
    {
        residual[i] = mResidual[i];
    }

    for (unsigned global_column=0; global_column<mSizeOfOdeSystem; global_column++)
    {
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            guess_perturbed[i] = rCurrentGuess[i];
        }

        guess_perturbed[global_column] += epsilon;

        ComputeResidual(pAbstractOdeSystem, timeStep, time, rCurrentYValues, guess_perturbed);
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            residual_perturbed[i] = mResidual[i];
        }

        // Compute residual_perturbed - residual
        double one_over_eps = 1.0/epsilon;
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            mJacobian[i][global_column] = one_over_eps*(residual_perturbed[i] - residual[i]);
        }
    }
}

void BackwardEulerIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                    double timeStep,
                                                    double time,
                                                    std::vector<double>& rCurrentYValues,
                                                    std::vector<double>& rNextYValues)
{
    // Check the size of the ODE system matches the solvers expected
    assert(mSizeOfOdeSystem == pAbstractOdeSystem->GetNumberOfStateVariables());

    unsigned counter = 0;
//        const double eps = 1e-6 * rCurrentGuess[0]; // Our tolerance (should use min(guess) perhaps?)
    const double eps = 1e-6; // JonW tolerance
    double norm = 2*eps;

    std::vector<double> current_guess(mSizeOfOdeSystem);
    current_guess.assign(rCurrentYValues.begin(), rCurrentYValues.end());

    while (norm > eps)
    {
        // Calculate Jacobian and mResidual for current guess
        ComputeResidual(pAbstractOdeSystem, timeStep, time, rCurrentYValues, current_guess);
        ComputeJacobian(pAbstractOdeSystem, timeStep, time, rCurrentYValues, current_guess);
//            // Update norm (our style)
//            norm = ComputeNorm(mResidual);

        // Solve Newton linear system
        SolveLinearSystem();

        // Update norm (JonW style)
        norm = ComputeNorm(mUpdate);

        // Update current guess
        for (unsigned i=0; i<mSizeOfOdeSystem; i++)
        {
            current_guess[i] -= mUpdate[i];
        }

        counter++;
        assert(counter < 20); // avoid infinite loops
    }
    rNextYValues.assign(current_guess.begin(), current_guess.end());
}

BackwardEulerIvpOdeSolver::BackwardEulerIvpOdeSolver(unsigned sizeOfOdeSystem)
{
    mSizeOfOdeSystem = sizeOfOdeSystem;

    // default epsilon
    mNumericalJacobianEpsilon = 1e-6;
    mForceUseOfNumericalJacobian = false;

    // allocate memory
    mResidual = new double[mSizeOfOdeSystem];
    mUpdate = new double[mSizeOfOdeSystem];

    mJacobian = new double*[mSizeOfOdeSystem];
    for (unsigned i=0; i<mSizeOfOdeSystem; i++)
    {
        mJacobian[i] = new double[mSizeOfOdeSystem];
    }
}

BackwardEulerIvpOdeSolver::~BackwardEulerIvpOdeSolver()
{
    // Delete pointers
    delete[] mResidual;
    delete[] mUpdate;

    for (unsigned i=0; i<mSizeOfOdeSystem; i++)
    {
        delete[] mJacobian[i];
    }
    delete[] mJacobian;
}

void BackwardEulerIvpOdeSolver::SetEpsilonForNumericalJacobian(double epsilon)
{
    assert(epsilon > 0);
    mNumericalJacobianEpsilon = epsilon;
}

void BackwardEulerIvpOdeSolver::ForceUseOfNumericalJacobian()
{
    mForceUseOfNumericalJacobian = true;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BackwardEulerIvpOdeSolver)
