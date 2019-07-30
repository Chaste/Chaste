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
#include "NhsContractionModel.hpp"
#include "OdeSystemInformation.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "UblasCustomFunctions.hpp"

#include <cmath>

//
// Model-scope constant parameters
//
const double NhsContractionModel::mKon = 100;
const double NhsContractionModel::mKrefoff = 0.2;
const double NhsContractionModel::mGamma = 2;
const double NhsContractionModel::mCalciumTroponinMax = 0.07;
const double NhsContractionModel::mAlphaR1 = 0.002;
const double NhsContractionModel::mAlphaR2 = 0.0017;
const double NhsContractionModel::mKZ = 0.15;
const unsigned NhsContractionModel::mNr = 3u;
const double NhsContractionModel::mBeta1 = -4;
const double NhsContractionModel::mAlpha0 = 0.008;
const unsigned NhsContractionModel::mN = 3u;
const double NhsContractionModel::mZp = 0.85;
const double NhsContractionModel::mCalcium50ref = 0.00105;
const double NhsContractionModel::mTref = 56.2;
const double NhsContractionModel::mBeta0 = 4.9;
const double NhsContractionModel::mA = 0.35;
const double NhsContractionModel::mA1 = -29;
const double NhsContractionModel::mA2 = 138;
const double NhsContractionModel::mA3 = 129;
const double NhsContractionModel::mAlpha1 = 0.03;
const double NhsContractionModel::mAlpha2 = 0.130;
const double NhsContractionModel::mAlpha3 = 0.625;


/*
 * ============================== PRIVATE FUNCTIONS =====================================
 */
void NhsContractionModel::CalculateCalciumTrop50()
{
    double Ca50ref_times_one_plus_beta1_times_lam_minus_one = mCalcium50ref * (1 + mBeta1*(mLambda-1));
    double one_plus_beta0_times_lam_minus_one_over_two_gamma = (1 + mBeta0*(mLambda-1))/(2*mGamma);

    mCalciumTrop50 = mCalciumTroponinMax * Ca50ref_times_one_plus_beta1_times_lam_minus_one;
    mCalciumTrop50 /= (Ca50ref_times_one_plus_beta1_times_lam_minus_one + (1-one_plus_beta0_times_lam_minus_one_over_two_gamma)*mKrefoff/mKon);
}


double NhsContractionModel::CalculateT0(double z)
{
    double calcium_ratio_to_n = SmallPow(mCalciumTrop50/mCalciumTroponinMax, mN);

    double z_max = mAlpha0 - mK2*calcium_ratio_to_n;
    z_max /= mAlpha0 + (mAlphaR1 + mK1)*calcium_ratio_to_n;

    return z * mTref * (1+mBeta0*(mLambda-1)) / z_max;
}

/*
 * ============================== PUBLIC FUNCTIONS =====================================
 */

NhsContractionModel::NhsContractionModel()
    :   AbstractOdeBasedContractionModel(5) // five state variables
{
    mpSystemInfo = OdeSystemInformation<NhsContractionModel>::Instance();
    ResetToInitialConditions();

    mLambda = 1.0;
    mDLambdaDt = 0.0;
    mCalciumI = 0.0;

    // Initialise mCalciumTrop50!!
    CalculateCalciumTrop50();

    double zp_to_n_plus_K_to_n = SmallPow(mZp,mNr) + SmallPow(mKZ,mNr);

    mK1 = mAlphaR2 * SmallPow(mZp,mNr-1) * mNr * SmallPow(mKZ,mNr);
    mK1 /= zp_to_n_plus_K_to_n * zp_to_n_plus_K_to_n;

    mK2 = mAlphaR2 * SmallPow(mZp,mNr)/zp_to_n_plus_K_to_n;
    mK2 *= 1 - mNr*SmallPow(mKZ,mNr)/zp_to_n_plus_K_to_n;
}

void NhsContractionModel::SetStretchAndStretchRate(double lambda, double dlambdaDt)
{
    assert(lambda>0.0);
    mLambda = lambda;
    mDLambdaDt = dlambdaDt;
    // lambda changed so update mCalciumTrop50!!
    CalculateCalciumTrop50();
}

void NhsContractionModel::SetInputParameters(ContractionModelInputParameters& rInputParameters)
{
    assert(rInputParameters.intracellularCalciumConcentration != DOUBLE_UNSET);
    assert(rInputParameters.intracellularCalciumConcentration > 0.0);
    mCalciumI = rInputParameters.intracellularCalciumConcentration;
}

void NhsContractionModel::SetIntracellularCalciumConcentration(double calciumConcentration)
{
    assert(calciumConcentration > 0.0);
    mCalciumI = calciumConcentration;
}

double NhsContractionModel::GetCalciumTroponinValue()
{
    return mStateVariables[0];
}

void NhsContractionModel::EvaluateYDerivatives(double time,
                                               const std::vector<double> &rY,
                                               std::vector<double> &rDY)
{
    //// if making changes here, see also NhsModelWithBackwardSolver::CalculateCaTropAndZDerivatives()

    const double& calcium_troponin = rY[0];
    const double& z = rY[1];
    const double& Q1 = rY[2];
    const double& Q2 = rY[3];
    const double& Q3 = rY[4];

    // check the state vars are in the expected range
    // LCOV_EXCL_START
    if (calcium_troponin < 0)
    {
        EXCEPTION("CalciumTrop concentration went negative");
    }
    if (z<0)
    {
        EXCEPTION("z went negative");
    }
    if (z>1)
    {
        EXCEPTION("z became greater than 1");
    }
    // LCOV_EXCL_STOP


    double Q = Q1 + Q2 + Q3;
    double T0 = CalculateT0(z);

    double Ta;
    if (Q>0)
    {
        Ta = T0*(1+(2+mA)*Q)/(1+Q);
    }
    else
    {
        Ta = T0*(1+mA*Q)/(1-Q);
    }

    rDY[0] =   mKon * mCalciumI * ( mCalciumTroponinMax - calcium_troponin)
             - mKrefoff * (1-Ta/(mGamma*mTref)) * calcium_troponin;

    double ca_trop_to_ca_trop50_ratio_to_n = SmallPow(calcium_troponin/mCalciumTrop50, mN);

    rDY[1] =   mAlpha0 * ca_trop_to_ca_trop50_ratio_to_n * (1-z)
             - mAlphaR1 * z
             - mAlphaR2 * SmallPow(z,mNr) / (SmallPow(z,mNr) + SmallPow(mKZ,mNr));


    rDY[2] = mA1 * mDLambdaDt - mAlpha1 * Q1;
    rDY[3] = mA2 * mDLambdaDt - mAlpha2 * Q2;
    rDY[4] = mA3 * mDLambdaDt - mAlpha3 * Q3;
}


double NhsContractionModel::GetActiveTension()
{
    double T0 = CalculateT0(mStateVariables[1]);
    double Q = mStateVariables[2]+mStateVariables[3]+mStateVariables[4];

    if (Q>0)
    {
        return T0*(1+(2+mA)*Q)/(1+Q);
    }
    else
    {
        return T0*(1+mA*Q)/(1-Q);
    }
}

template<>
void OdeSystemInformation<NhsContractionModel>::Initialise(void)
{
    this->mVariableNames.push_back("CalciumTroponin");
    this->mVariableUnits.push_back("microMols");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("z");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Q1");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Q2");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("Q3");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0);

    this->mInitialised = true;
}
