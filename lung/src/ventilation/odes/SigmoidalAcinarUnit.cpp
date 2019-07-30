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

#include "SigmoidalAcinarUnit.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "MathsCustomFunctions.hpp"


SigmoidalAcinarUnit::SigmoidalAcinarUnit() : mQ(0.0),
                                             mDt(-1),
                                             mPaw(0.0),
                                             mPpl(0.0),
                                             mRaw(0.0),
                                             mV(0.0)

{
    mA = 2/1e3; //m^3 (RV 2L)
    mB = (6 - 2)/1e3; //m^3 (TLC 6L, RV 2L)
    mC = 667; //Pa (6.8 cmH2O)
    mD = 300; //Pa (3.8 cmH2O)
}


SigmoidalAcinarUnit::~SigmoidalAcinarUnit()
{
}


void SigmoidalAcinarUnit::SetTimestep(double dt)
{
}


void SigmoidalAcinarUnit::SolveAndUpdateState(double tStart, double tEnd)
{
    double dt = tEnd - tStart;

    double compliance = CalculateCurrentCompliance();

    mQ = (1 - std::exp(-dt/(mRaw*compliance)))*mQ*mRaw*compliance/dt;
    mV += dt*mQ;
    mPaw = (mV - mA)/compliance + mPpl;
}

void SigmoidalAcinarUnit::ComputeExceptFlow(double tStart, double tEnd)
{
    //double compliance = CalculateCurrentCompliance();
    //mPaw = (mV - mA)/compliance + mPpl-mC;
    mPaw = -mD*log(mB/(mV - mA) -1) + mPpl + mC;
}


void SigmoidalAcinarUnit::SetFlow(double flow)
{
    mQ = flow;
}

double SigmoidalAcinarUnit::GetFlow()
{
    return mQ;
}

void SigmoidalAcinarUnit::UpdateFlow(double tStart, double tEnd)
{
    double dt = tEnd - tStart;
    double compliance_factor = CalculateCurrentDerivativeVolumeOverCompliance();
    //double compliance = CalculateCurrentCompliance();

    mQ = (1 - std::exp(-dt*compliance_factor/mRaw))*mQ*mRaw/compliance_factor/dt;
    mV += dt*mQ;
}

void SigmoidalAcinarUnit::SetAirwayPressure(double pressure)
{
    mPaw = pressure;
}

void SigmoidalAcinarUnit::SetPleuralPressure(double pressure)
{
    mPpl = pressure;
}

double SigmoidalAcinarUnit::GetAirwayPressure()
{
    return mPaw;
}

void SigmoidalAcinarUnit::SetTerminalBronchioleResistance(double raw)
{
    mRaw = raw;
}

double SigmoidalAcinarUnit::GetVolume()
{
    return mV;
}

void SigmoidalAcinarUnit::SetUndeformedVolume(double v)
{
    mV = v;
}

double SigmoidalAcinarUnit::GetStretchRatio()
{
    return 0.0;
}

void SigmoidalAcinarUnit::SetStretchRatio(double lambda)
{
}

double SigmoidalAcinarUnit::CalculateCurrentCompliance()
{
    assert(mV > mA);
    assert(mV < (mA+mB));
    return -(mA-mV)*(mA+mB-mV)/(mB*mD);
}

double SigmoidalAcinarUnit::CalculateCurrentDerivativeVolumeOverCompliance()
{
    assert(mV > mA);
    assert(mV < (mA+mB));
    return -(mB*mD)/((mA-mV)*(mA+mB-mV));
}

void SigmoidalAcinarUnit::SetA(double a)
{
    mA = a;
}

void SigmoidalAcinarUnit::SetB(double b)
{
    mB = b;
}

void SigmoidalAcinarUnit::SetC(double c)
{
    mC = c;
}

void SigmoidalAcinarUnit::SetD(double d)
{
    mD = d;
}
