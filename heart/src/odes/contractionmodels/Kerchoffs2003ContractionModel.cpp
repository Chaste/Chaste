/*

Copyright (C) University of Oxford, 2005-2012

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
const double Kerchoffs2003ContractionModel::tr = 75; // ms
const double Kerchoffs2003ContractionModel::td = 75; // ms
const double Kerchoffs2003ContractionModel::b = 150; // ms/um
const double Kerchoffs2003ContractionModel::ld = -0.4; // um

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
    if(lc > a7)
    {
        f_iso = T0 * pow(tanh(a6*(lc-a7)),2);
    }

    double f_twitch = 0;
    double t_max = b*(mSarcomereLength - ld);
    if(mIsActivated)
    {
        double t_a = mTime - mActivationTime;

        if(t_a < t_max)
        {
            f_twitch = pow( tanh(t_a/tr)*tanh((t_max-t_a)/td), 2);
        }
        else if(mElectricallyUnactivated)
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
    this->mInitialised = true;
}
