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

#include "NhsModelWithBackwardSolver.hpp"
#include <iostream>
#include <cmath>
#include "LogFile.hpp"
#include "Exception.hpp"
#include "TimeStepper.hpp"

const double NhsModelWithBackwardSolver::mTolerance = 1e-10;

double NhsModelWithBackwardSolver::ImplicitSolveForQ()
{
    mTemporaryStateVariables[2] = (mTemporaryStateVariables[2] + mDt*mA1*mDLambdaDt)/(1 + mAlpha1*mDt);
    mTemporaryStateVariables[3] = (mTemporaryStateVariables[3] + mDt*mA2*mDLambdaDt)/(1 + mAlpha2*mDt);
    mTemporaryStateVariables[4] = (mTemporaryStateVariables[4] + mDt*mA3*mDLambdaDt)/(1 + mAlpha3*mDt);

    return mTemporaryStateVariables[2] + mTemporaryStateVariables[3] + mTemporaryStateVariables[4];
}

void NhsModelWithBackwardSolver::CalculateCaTropAndZDerivatives(double calciumTroponin, double z, double Q,
                                                                double& dCaTrop, double& dz)
{
//As in straight Nhs, we don't cover the exception code
#define COVERAGE_IGNORE
    if(calciumTroponin < 0)
    {
        EXCEPTION("CalciumTrop concentration went negative");
    }
    if(z<0)
    {
        EXCEPTION("z went negative");
    }
    if(z>1)
    {
        EXCEPTION("z became greater than 1");
    }
#undef COVERAGE_IGNORE

    double T0 = CalculateT0(z);

    double Ta;
    if(Q>0)
    {
        Ta = T0*(1+(2+mA)*Q)/(1+Q);
    }
    else
    {
        Ta = T0*(1+mA*Q)/(1-Q);
    }

    dCaTrop =   mKon * mCalciumI * ( mCalciumTroponinMax - calciumTroponin)
             - mKrefoff * (1-Ta/(mGamma*mTref)) * calciumTroponin;

    double ca_trop_to_ca_trop50_ratio_to_n = pow(calciumTroponin/mCalciumTrop50, mN);

    dz =   mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n * (1-z)
         - mAlphaR1 * z
         - mAlphaR2 * pow(z,mNr) / (pow(z,mNr) + pow(mKZ,mNr));
}



void NhsModelWithBackwardSolver::CalculateBackwardEulerResidual(double calciumTroponin, double z, double Q,
                                                                double& residualComponent1, double& residualComponent2)
{
    double dcatrop;
    double dz;
    CalculateCaTropAndZDerivatives(calciumTroponin,z,Q,dcatrop,dz);

    residualComponent1 = calciumTroponin - mDt*dcatrop - mTemporaryStateVariables[0];
    residualComponent2 = z - mDt*dz - mTemporaryStateVariables[1];
}



NhsModelWithBackwardSolver::NhsModelWithBackwardSolver()
{
    mTemporaryStateVariables.resize(5);
}



void NhsModelWithBackwardSolver::RunDoNotUpdate(double startTime, double endTime, double timestep)
{
    assert(startTime < endTime);

    mDt = timestep;

    mTemporaryStateVariables = mStateVariables;

    // loop in time
    TimeStepper stepper(startTime, endTime, timestep);

    while ( !stepper.IsTimeAtEnd() )
    {
        /////////////////////////////////////////////////////////
        // Q1,Q2,Q3 using backward euler can solved straightaway
        /////////////////////////////////////////////////////////
        double new_Q = ImplicitSolveForQ();

        ////////////////////////////////////////////////////////////////////
        // Solve the 2D nonlinear problem for Backward Euler Ca_trop and z
        ////////////////////////////////////////////////////////////////////

        // see what the residual is
        double catrop_guess = mTemporaryStateVariables[0];
        double z_guess = mTemporaryStateVariables[1];
        double f1,f2; // f=[f1,f2]=residual

        CalculateBackwardEulerResidual(catrop_guess, z_guess, new_Q, f1, f2);
        double norm_resid = sqrt(f1*f1+f2*f2);

        // solve using Newton's method, no damping. Stop if num iterations
        // reaches 15 (very conservative)
        unsigned counter = 0;
        while ((norm_resid>mTolerance) && (counter++<15))
        {
            // numerically approximate the jacobian J
            double j11,j12,j21,j22; // J = [j11, j12; j21 j22]
            double temp1,temp2;

            double h = std::max(fabs(catrop_guess/100),1e-8);
            CalculateBackwardEulerResidual(catrop_guess+h, z_guess, new_Q, temp1, temp2);
            j11 = (temp1-f1)/h;
            j21 = (temp2-f2)/h;

            h = std::max(fabs(z_guess/100),1e-8);
            CalculateBackwardEulerResidual(catrop_guess, z_guess+h, new_Q, temp1, temp2);
            j12 = (temp1-f1)/h;
            j22 = (temp2-f2)/h;

            // compute u = J^{-1} f (exactly, as a 2D problem)
            double one_over_det = 1.0/(j11*j22-j12*j21);
            double u1 = one_over_det*(j22*f1  - j12*f2);
            double u2 = one_over_det*(-j21*f1 + j11*f2);

            catrop_guess -= u1;
            z_guess -= u2;

            CalculateBackwardEulerResidual(catrop_guess, z_guess, new_Q, f1, f2);
            norm_resid = sqrt(f1*f1+f2*f2);
        }
        assert(counter<15); // if this fails, see corresponding code in old NhsModelWithImplicitSolver

        mTemporaryStateVariables[0] = catrop_guess;
        mTemporaryStateVariables[1] = z_guess;

        stepper.AdvanceOneTimeStep();
    }
}

double NhsModelWithBackwardSolver::GetNextActiveTension()
{
    double T0 = CalculateT0(mTemporaryStateVariables[1]);
    double Q = mTemporaryStateVariables[2]+mTemporaryStateVariables[3]+mTemporaryStateVariables[4];

    if(Q>0)
    {
        return T0*(1+(2+mA)*Q)/(1+Q);
    }
    else
    {
        return T0*(1+mA*Q)/(1-Q);
    }
}

void NhsModelWithBackwardSolver::RunAndUpdate(double startTime, double endTime, double timestep)
{
    RunDoNotUpdate(startTime, endTime, timestep);
    UpdateStateVariables();
}

