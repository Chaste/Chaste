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

#include "SimpleBalloonExplicitAcinarUnit.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>


SimpleBalloonExplicitAcinarUnit::SimpleBalloonExplicitAcinarUnit() : mQ(0.0),
                                                     mDt(-1),
                                                     mPaw(0.0),
                                                     mPpl(0.0),
                                                     mRaw(0.0),
                                                     mCa(0.0),
                                                     mV0(0.0)

{

}


SimpleBalloonExplicitAcinarUnit::~SimpleBalloonExplicitAcinarUnit()
{
}


void SimpleBalloonExplicitAcinarUnit::SetTimestep(double dt)
{
}


void SimpleBalloonExplicitAcinarUnit::SolveAndUpdateState(double tStart, double tEnd)
{
    double dt = tEnd - tStart;

    mV0 += dt*mQ;
    mPaw = mV0/mCa + mPpl;
}

void SimpleBalloonExplicitAcinarUnit::ComputeExceptFlow(double tStart, double tEnd)
{
    mPaw = mV0/mCa + mPpl;
}


void SimpleBalloonExplicitAcinarUnit::SetFlow(double flow)
{
    mQ = flow;
}

double SimpleBalloonExplicitAcinarUnit::GetFlow()
{
    return mQ;
}

void SimpleBalloonExplicitAcinarUnit::UpdateFlow(double tStart, double tEnd)
{
    double dt = tEnd - tStart;
    mV0 += dt*mQ;
}

void SimpleBalloonExplicitAcinarUnit::SetAirwayPressure(double pressure)
{
    mPaw = pressure;
}

void SimpleBalloonExplicitAcinarUnit::SetPleuralPressure(double pressure)
{
    mPpl = pressure;
}

double SimpleBalloonExplicitAcinarUnit::GetAirwayPressure()
{
    return mPaw;
}

void SimpleBalloonExplicitAcinarUnit::SetTerminalBronchioleResistance(double raw)
{
    mRaw = raw;
}

double SimpleBalloonExplicitAcinarUnit::GetVolume()
{
    return mV0;
}

void SimpleBalloonExplicitAcinarUnit::SetUndeformedVolume(double v0)
{
    mV0 = v0;
}

double SimpleBalloonExplicitAcinarUnit::GetStretchRatio()
{
    return 0.0;
}

void SimpleBalloonExplicitAcinarUnit::SetStretchRatio(double lambda)
{
}

void SimpleBalloonExplicitAcinarUnit::SetCompliance(double ca)
{
    mCa = ca;
}
