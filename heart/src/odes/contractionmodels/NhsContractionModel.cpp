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
    #define COVERAGE_IGNORE
    if(calcium_troponin < 0)
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


    double Q = Q1 + Q2 + Q3;
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

    if(Q>0)
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
