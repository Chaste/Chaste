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
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include <cmath>

//
// Model-scope constant parameters
//
const double FitzHughNagumo1961OdeSystem::mAlpha = -0.08;
const double FitzHughNagumo1961OdeSystem::mGamma = 3.00;
const double FitzHughNagumo1961OdeSystem::mEpsilon = 0.005;


FitzHughNagumo1961OdeSystem::FitzHughNagumo1961OdeSystem(
        boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCell(pOdeSolver, 2, 0, pIntracellularStimulus)
{
    mpSystemInfo = OdeSystemInformation<FitzHughNagumo1961OdeSystem>::Instance();

    Init();
}

FitzHughNagumo1961OdeSystem::~FitzHughNagumo1961OdeSystem(void)
{
}

void FitzHughNagumo1961OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
{
    double membrane_V = rY[0]; // v
    double recovery_variable = rY[1]; // w

    double i_stim = GetIntracellularAreaStimulus(time);

    // dV/dt
    double membrane_V_prime = 0;
    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
    if (!mSetVoltageDerivativeToZero)
    {
        membrane_V_prime = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)-recovery_variable+i_stim;
    }

    // dw/dt
    double recovery_variable_prime = mEpsilon*(membrane_V-mGamma*recovery_variable);

    rDY[0] = membrane_V_prime;
    rDY[1] = recovery_variable_prime;
}

double FitzHughNagumo1961OdeSystem::GetIIonic(const std::vector<double>* pStateVariables)
{
    if (!pStateVariables) pStateVariables = &mStateVariables;
    double membrane_V = (*pStateVariables)[mVoltageIndex];
    double recovery_variable = (*pStateVariables)[1];
    double fake_ionic_current = membrane_V*(membrane_V-mAlpha)*(1-membrane_V)-recovery_variable;
    return fake_ionic_current;
}

template<>
void OdeSystemInformation<FitzHughNagumo1961OdeSystem>::Initialise(void)
{
    /*
     * State variables
     */
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("w");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}
