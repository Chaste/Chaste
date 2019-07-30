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
#include "HiornsAirwayWall.hpp"
#include "MathsCustomFunctions.hpp"
#include "Exception.hpp"


HiornsAirwayWall::HiornsAirwayWall() : mTargetPressure(0),
                                       mRIn(0),
                                       mROut(0),
                                       mmu(0),
                                       mphi1(0),
                                       mphi2(0),
                                       mC1(0),
                                       mC2(0),
                                       mA(0)
{
}

HiornsAirwayWall::~HiornsAirwayWall() {}

void HiornsAirwayWall::SetTimestep(double dt) {}

double HiornsAirwayWall::CalculatePressureRadiusResidual(double radius)
{

    mTargetPressure = mAirwayPressure - mPleuralPressure;

    double rin = radius;

    double areaOfAirwayWall = M_PI*(mROut*mROut - mRIn*mRIn);
    double rout = sqrt(rin*rin + areaOfAirwayWall/M_PI);
    double pressure;

    pressure = mmu*log((rin*mROut)/(rout*mRIn)) + mmu*((rin*rin - mRIn*mRIn)/2.)*((1./rin*rin) - (1./rout*rout)) + 2*cos(mphi2)*cos(mphi2)*mA*log(mROut/mRIn);

    if (rin - mRIn > 0.)
        {

        //  THE COMMENTED OUT CODE CORRESPONDS TO THE PRELIMINARY VERSION OF THE HIORNS QUASI-STATIC
        //  AIRWAY WALL MODEL (I.E. EQUATION 1.1 IN THE HIORNS LOOKUP TABLES NOTES). THE VERSION THAT
        //  WE USE (FOLLOWING THE COMMENTED SECTION) INSTEAD IS MODIFIED SO THAT THE CONTRIBUTION OF
        //  THE COLLAGEN TO THE STRAIN-ENERGY LAW SATISFIES THE MODEL OF HOLZAPFEL ET AL.

        //  const int numberTrap = 10000;
        //  double integrandVals[numberTrap];
        //  double tVals[numberTrap];
        //  double integralStuff = 0.;

        //  double upperlimit = -1.;
        //  double lowerlimit = -1.;

        //  lowerlimit = ((sqrt(mC2)*(rout*rout - mROut*mROut)*cos(mphi1)*cos(mphi1))/((mROut)*(mROut)));
        //  upperlimit = ((sqrt(mC2)*(rin*rin - mRIn*mRIn)*cos(mphi1)*cos(mphi1))/((mRIn)*(mRIn)));

        //  for (int it = 0; it < numberTrap; it++)
        //  {

            //  double tVal = lowerlimit + ((double)it/((double)numberTrap - 1.))*(upperlimit - lowerlimit);
            //  tVals[it] = tVal;
            //  integrandVals[it] = (2.*exp(tVal*tVal)/sqrt(M_PI));

        //  }

        //  integralStuff = (0.5*(integrandVals[0] + integrandVals[numberTrap - 1]));
        //  integralStuff = integralStuff*(tVals[1] - tVals[0]);
        //  for (int i = 1; i < (numberTrap - 1); i++)
        //  {
            //  integralStuff = integralStuff + integrandVals[i]*(tVals[1] - tVals[0]);
        //  }

        //  pressure = pressure + mC1*sqrt(M_PI/mC2)*cos(mphi1)*cos(mphi1)*integralStuff;

            double integralStuff = exp(mC2*((rin*rin - mRIn*mRIn)/(mRIn*mRIn))*((rin*rin - mRIn*mRIn)/(mRIn*mRIn))*cos(mphi1)*cos(mphi1)*cos(mphi1)*cos(mphi1)) - exp(mC2*((rout*rout - mROut*mROut)/(mROut*mROut))*((rout*rout - mROut*mROut)/(mROut*mROut))*cos(mphi1)*cos(mphi1)*cos(mphi1)*cos(mphi1));
            pressure = pressure + (mC1/mC2)*cos(mphi1)*cos(mphi1)*integralStuff;

        }

    double residual = mTargetPressure - pressure;

    return residual;

}

void HiornsAirwayWall::SolveAndUpdateState(double tStart, double tEnd)
{

    double guess = (mRIn + mROut)/2.;
    double factor = 2.;

    Tolerance tol = 0.000001;
    boost::uintmax_t maxIterations = 500u;

    std::pair<double, double> found = boost::math::tools::bracket_and_solve_root(boost::bind(&HiornsAirwayWall::CalculatePressureRadiusResidual, this, _1), guess, factor, false, tol, maxIterations);
    mDeformedAirwayRadius = found.first;
}

void HiornsAirwayWall::SetRIn(double RIn)
{
    assert (RIn >= 0.0);
    mRIn = RIn;
}

void HiornsAirwayWall::SetROut(double ROut)
{
    assert (ROut >= 0.0);
    mROut = ROut;
}

void HiornsAirwayWall::Setmu(double mu)
{
    assert (mu >= 0.0);
    mmu = mu;
}

void HiornsAirwayWall::Setphi1(double phi1)
{
    mphi1 = phi1;
}

void HiornsAirwayWall::Setphi2(double phi2)
{
    mphi2 = phi2;
}

void HiornsAirwayWall::SetC1(double C1)
{
    mC1 = C1;
}

void HiornsAirwayWall::SetC2(double C2)
{
    mC2 = C2;
}

void HiornsAirwayWall::SetA(double A)
{
    mA = A;
}

double HiornsAirwayWall::GetLumenRadius()
{
    return mDeformedAirwayRadius;
}

void HiornsAirwayWall::SetAirwayPressure(double pressure)
{
    mAirwayPressure = pressure;
}

void HiornsAirwayWall::SetPleuralPressure(double pressure)
{
    mPleuralPressure = pressure;
}
