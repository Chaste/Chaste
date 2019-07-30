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

#include <cmath>
#include <iostream>
#include <boost/math/tools/roots.hpp>
#include <boost/bind.hpp>

#include "BoostTolerance.hpp"
#include "LaPradAirwayWall.hpp"
#include "MathsCustomFunctions.hpp"
#include "Exception.hpp"

LaPradAirwayWall::LaPradAirwayWall() : mTargetPressure(0),
                                       mRIn(0),
                                       mROut(0),
                                       mk1(0),
                                       mk2(0),
                                       mk3(0)
{
}

LaPradAirwayWall::~LaPradAirwayWall() {}

void LaPradAirwayWall::SetTimestep(double dt) {}

double LaPradAirwayWall::CalculatePressureRadiusResidual(double radius)
{

    mTargetPressure = mAirwayPressure - mPleuralPressure;

    double rin = radius;

    double areaOfAirwayWall = M_PI*(mROut*mROut - mRIn*mRIn);
    double rout = sqrt(rin*rin + areaOfAirwayWall/M_PI);
    long double functionValues[10000];
    double rValues[10000];
    double pressure;

    for (int i = 0; i < 10000; i++)
    {

        double RVal = mRIn + (double)i*(mROut - mRIn)/(10000. - 1.);
        rValues[i] = rin + (double)i*(rout - rin)/(10000. - 1.);
        functionValues[i] = ((rValues[i]/RVal)*(rValues[i]/RVal) - (RVal/rValues[i])*(RVal/rValues[i]))*(mk1*sqrt(1. + (RVal/rValues[i])*(RVal/rValues[i]) + (rValues[i]/RVal)*(rValues[i]/RVal) - 3.) + mk2*sqrt(1 + (RVal/rValues[i])*(RVal/rValues[i]) + (rValues[i]/RVal)*(rValues[i]/RVal) - 3.)*exp(mk3*(1. + (RVal/rValues[i])*(RVal/rValues[i]) + (rValues[i]/RVal)*(rValues[i]/RVal) - 3.)*(1 + (RVal/rValues[i])*(RVal/rValues[i]) + (rValues[i]/RVal)*(rValues[i]/RVal) - 3.)))/rValues[i];
    }


    pressure = (0.5*(functionValues[0] + functionValues[10000 - 1]));
    for (int i = 1; i < (10000 - 1); i++)
    {
        pressure = pressure + functionValues[i];
    }
    pressure = pressure*(rValues[1] - rValues[0]);

    double residual = mTargetPressure - pressure;

    return residual;

}

void LaPradAirwayWall::SolveAndUpdateState(double tStart, double tEnd)
{

    double guess = (mRIn + mROut)/2.;
    double factor = 2.;

    Tolerance tol = 0.000001;
    boost::uintmax_t maxIterations = 500u;

    std::pair<double, double> found = boost::math::tools::bracket_and_solve_root(boost::bind(&LaPradAirwayWall::CalculatePressureRadiusResidual, this, _1), guess, factor, false, tol, maxIterations);
    mDeformedAirwayRadius = found.first;
}

void LaPradAirwayWall::SetRIn(double RIn)
{
    assert (RIn >= 0.0);
    mRIn = RIn;
}

void LaPradAirwayWall::SetROut(double ROut)
{
    assert (ROut >= 0.0);
    mROut = ROut;
}

void LaPradAirwayWall::Setk1(double k1)
{
    assert (k1 >= 0.0);
    mk1 = k1;
}

void LaPradAirwayWall::Setk2(double k2)
{
    assert (k2 >= 0.);
    mk2 = k2;
}

void LaPradAirwayWall::Setk3(double k3)
{
    assert (k3 >= 0.);
    mk3 = k3;
}

double LaPradAirwayWall::GetLumenRadius()
{
    return mDeformedAirwayRadius;
}

void LaPradAirwayWall::SetAirwayPressure(double pressure)
{
    mAirwayPressure = pressure;
}

void LaPradAirwayWall::SetPleuralPressure(double pressure)
{
    mPleuralPressure = pressure;
}
