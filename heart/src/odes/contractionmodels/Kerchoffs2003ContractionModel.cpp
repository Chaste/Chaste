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


#include "Kerchoffs2003ContractionModel.hpp"
#include "Exception.hpp"
#include "TimeStepper.hpp"
#include <iostream>

const double Kerchoffs2003ContractionModel::a6 = 2.0; // 1/um
const double Kerchoffs2003ContractionModel::a7 = 1.5; // um
const double Kerchoffs2003ContractionModel::T0 = 180; // kPa
const double Kerchoffs2003ContractionModel::Ea = 20;  // 1/um
const double Kerchoffs2003ContractionModel::v0 = 0.0075; // um/ms
const double Kerchoffs2003ContractionModel::ls0 = 1.9; // um
const double Kerchoffs2003ContractionModel::ld = -0.4; // um
const double Kerchoffs2003ContractionModel::mActivationVoltage = 0.0;
const double Kerchoffs2003ContractionModel::mDeactivationVoltage = -70.0;

Kerchoffs2003ContractionModel::Kerchoffs2003ContractionModel()
    : AbstractOdeBasedContractionModel(1)
{
    mpSystemInfo = OdeSystemInformation<Kerchoffs2003ContractionModel>::Instance();

    mSarcomereLength = ls0;

    mStateVariables.push_back(mSarcomereLength-1.0/Ea); //steady state

    mIsActivated = false;
    mElectricallyUnactivated = true;
    mActivationTime = 0.0;
    mTime = 0.0;

    this->mParameters.resize(3);
    SetParameter("tr", 75.0); // ms
    SetParameter("td", 75.0); // ms
    SetParameter("b", 150.0); // ms/um
}


void Kerchoffs2003ContractionModel::EvaluateYDerivatives(double time,
                                                         const std::vector<double>& rY,
                                                         std::vector<double>& rDY)
{
    double lc = rY[0];
    rDY[0]=( Ea*(mSarcomereLength-lc) - 1 )*v0;
}


void Kerchoffs2003ContractionModel::SetInputParameters(ContractionModelInputParameters& rInputParameters)
{
    assert(rInputParameters.voltage != DOUBLE_UNSET);

    if (mIsActivated && (rInputParameters.voltage < mDeactivationVoltage))
    {
        // inactive (resting) - note don't set mIsActivated=false yet
        // as the cell may yet be producing force, and the code is such
        // that if mIsActivated=false, Ta=0
        mElectricallyUnactivated = true;
    }

    if (!mIsActivated && (rInputParameters.voltage > mActivationVoltage))
    {
        // activated
        mIsActivated = true;
        mActivationTime = mTime;
        mElectricallyUnactivated = false;
    }
}

void Kerchoffs2003ContractionModel::SetStretchAndStretchRate(double stretch, double stretchRate)
{
    mSarcomereLength = stretch*ls0;
}


double Kerchoffs2003ContractionModel::GetActiveTension(double lc)
{
    double f_iso = 0;
    if (lc > a7)
    {
        f_iso = T0 * pow(tanh(a6*(lc-a7)),2);
    }

    double f_twitch = 0;
    double b = this->GetParameter("b");
    double t_max = b*(mSarcomereLength - ld);
    if (mIsActivated)
    {
        double t_a = mTime - mActivationTime;

        if (t_a < t_max)
        {
            double tr = this->GetParameter("tr");
            double td = this->GetParameter("td");
            f_twitch = pow( tanh(t_a/tr)*tanh((t_max-t_a)/td), 2);
        }
        else if (mElectricallyUnactivated)
        {
            // t_a < t_ma => f_twitch=0 => Ta=0
            // In this case, if electrically unactivated as well,
            // set the state to be unactivated.
            mIsActivated = false;
        }
    }

    // expl is unstable for dt = 0.01, 0.001, impl is fine
    return (mSarcomereLength/ls0)*f_iso*f_twitch*(mSarcomereLength-lc)*Ea;
}


double Kerchoffs2003ContractionModel::GetActiveTension()
{
    return GetActiveTension(mStateVariables[0]);
}

double Kerchoffs2003ContractionModel::GetNextActiveTension()
{
    return GetActiveTension(mTemporaryStateVariables[0]);
}

template<>
void OdeSystemInformation<Kerchoffs2003ContractionModel>::Initialise()
{
    this->mVariableNames.push_back("lc");
    this->mVariableUnits.push_back("um");

    this->mParameterNames.push_back("tr");
    this->mParameterUnits.push_back("ms");
    this->mParameterNames.push_back("td");
    this->mParameterUnits.push_back("ms");
    this->mParameterNames.push_back("b");
    this->mParameterUnits.push_back("ms/um");

    this->mInitialised = true;
}
