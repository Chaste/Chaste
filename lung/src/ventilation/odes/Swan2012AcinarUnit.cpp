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

#include "Swan2012AcinarUnit.hpp"
#include "MathsCustomFunctions.hpp"
#include <cmath>
#include <iostream>

Swan2012AcinarUnit::Swan2012AcinarUnit() : mQ(0.0),
                                           mLambda(1.26),
                                           mDt(-1),
                                           mPaw(0.0),
                                           mPpl(0.0),
                                           mRaw(0.0),
                                           mA(0.433),
                                           mB(-0.611),
                                           mXi(2500),
                                           mV0(1.0)
{

}


Swan2012AcinarUnit::~Swan2012AcinarUnit()
{
}


void Swan2012AcinarUnit::SetTimestep(double dt)
{
}


void Swan2012AcinarUnit::SolveAndUpdateState(double tStart, double tEnd)
{

}

void Swan2012AcinarUnit::ComputeExceptFlow(double tStart, double tEnd)
{
    double Pe = CalculateStaticRecoilPressure(mLambda);

    mPaw = Pe + mPpl;
}

void Swan2012AcinarUnit::UpdateFlow(double tStart, double tEnd)
{
    double dt = tEnd - tStart;
    double compliance_factor = CalculateDerivativeStaticRecoilPressureByStrain();

    mQ = (1 - std::exp(-dt*compliance_factor/mRaw))*mQ*mRaw/compliance_factor/dt;

    double V = GetVolume();
    V += dt*mQ;
    mLambda = std::pow(V/mV0, 1.0/3.0);
}

void Swan2012AcinarUnit::SetFlow(double flow)
{
    mQ = flow;
}

double Swan2012AcinarUnit::GetFlow()
{
    return mQ;
}

void Swan2012AcinarUnit::SetAirwayPressure(double pressure)
{
    mPaw = pressure;
}

double Swan2012AcinarUnit::GetAirwayPressure()
{
    return mPaw;
}

void Swan2012AcinarUnit::SetPleuralPressure(double pressure)
{
    mPpl = pressure;
}

void Swan2012AcinarUnit::SetTerminalBronchioleResistance(double raw)
{
    mRaw = raw;
}

double Swan2012AcinarUnit::GetStretchRatio()
{
    return mLambda;
}

void Swan2012AcinarUnit::SetStretchRatio(double lambda)
{
    mLambda = lambda;
}


double Swan2012AcinarUnit::GetVolume()
{
    return mLambda*mLambda*mLambda*mV0;
}

void Swan2012AcinarUnit::SetUndeformedVolume(double v0)
{
    mV0 = v0;
}

double Swan2012AcinarUnit::CalculateDerivativeVolumeByStrain()
{
    return 3*mLambda*mLambda*mV0;
}

double Swan2012AcinarUnit::CalculateDerivativeStaticRecoilPressureByStrain()
{
    double gamma = CalculateGamma(mLambda);
    return ((3.0*mXi/2.0)*(3*mA + mB)*(3*mA + mB)*(mLambda*mLambda - 1)*(mLambda*mLambda - 1)*std::exp(gamma) +
            (mXi/2.0)*(3*mA + mB)*(mLambda*mLambda + 1)*std::exp(gamma)/(mLambda*mLambda))/CalculateDerivativeVolumeByStrain();
}

double Swan2012AcinarUnit::CalculateGamma(double lambda)
{
    return (3.0/4.0)*(3*mA + mB)*(lambda*lambda - 1)*(lambda*lambda - 1);
}

double Swan2012AcinarUnit::CalculateStaticRecoilPressure(double lambda)
{
    double gamma = CalculateGamma(lambda);
    return mXi*std::exp(gamma)/(2.0*lambda)*(3*mA + mB)*(lambda*lambda -1);
}
