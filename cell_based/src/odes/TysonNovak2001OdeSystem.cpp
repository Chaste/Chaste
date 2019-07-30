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

#include "TysonNovak2001OdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "MathsCustomFunctions.hpp"

TysonNovak2001OdeSystem::TysonNovak2001OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystemWithAnalyticJacobian(6)
{
    mpSystemInfo = OdeSystemInformation<TysonNovak2001OdeSystem>::Instance();

    Init();

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

TysonNovak2001OdeSystem::~TysonNovak2001OdeSystem()
{
    // Do nothing
}

void TysonNovak2001OdeSystem::Init()
{
    // Initialise model parameter values
    mK1 = 0.04;
    mK2d = 0.04;
    mK2dd = 1.0;
    mK2ddd = 1.0;
    mCycB_threshold = 0.1;
    mK3d = 1.0;
    mK3dd = 10.0;
    mK4d = 2.0;
    mK4 = 35;
    mJ3 = 0.04;
    mJ4 = 0.04;
    mK5d = 0.005;
    mK5dd = 0.2;
    mK6 = 0.1;
    mJ5 = 0.3;
    mN = 4u;
    mK7 = 1.0;
    mK8 = 0.5;
    mJ7 = 1e-3;
    mJ8 = 1e-3;
    mMad = 1.0;
    mK9 = 0.1;
    mK10 = 0.02;
    mK11 = 1.0;
    mK12d = 0.2;
    mK12dd = 50.0;
    mK12ddd = 100.0;
    mKeq = 1e3;
    mK13 = 1.0;
    mK14 = 1.0;
    mK15d = 1.5;
    mK15dd = 0.05;
    mK16d = 1.0;
    mK16dd = 3.0;
    mJ15 = 0.01;
    mJ16 = 0.01;
    mMu = 0.01;
    mMstar = 10.0;
}

void TysonNovak2001OdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double x1 = rY[0];
    double x2 = rY[1];
    double x3 = rY[2];
    double x4 = rY[3];
    double x5 = rY[4];
    double x6 = rY[5];

    double dx1 = 0.0;
    double dx2 = 0.0;
    double dx3 = 0.0;
    double dx4 = 0.0;
    double dx5 = 0.0;
    double dx6 = 0.0;

    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;

    /**
     * 1. [CycB]
     * 2. [Cdh1]
     * 3. [Cdc20T]
     * 4. [Cdc20A]
     * 5. [IEP]
     * 6. m - mass of the cell
     */

    dx1 = mK1-(mK2d+mK2dd*x2)*x1;

    // The commented line below models the start transition, no cycling, without Cdc20A
//    temp1 = ((mK3d)*(1.0-x2))/(mJ3+1.0-x2);

    temp1 = ((mK3d+mK3dd*x4)*(1.0-x2))/(mJ3+1.0-x2);
    temp2 = (mK4*x6*x1*x2)/(mJ4+x2);
    dx2 = temp1-temp2;

    temp1 = mK5dd*(SmallPow(x1*x6/mJ5,mN)/(1+SmallPow(x1*x6/mJ5,mN)));
    temp2 = mK6*x3;
    dx3 = mK5d + temp1 - temp2;

    temp1 = (mK7*x5*(x3-x4))/(mJ7+x3-x4);
    temp2 = (mK8*mMad*x4)/(mJ8+x4);
    temp3 = mK6*x4;
    dx4 = temp1 - temp2 - temp3;

    dx5 = mK9*x6*x1*(1.0-x5) - mK10*x5;

    dx6 = mMu*x6*(1.0-x6/mMstar);

    // Multiply by 60 beacuase the Tyson and Novak 2001 paper has time in minutes, not hours
    rDY[0] = dx1*60.0;
    rDY[1] = dx2*60.0;
    rDY[2] = dx3*60.0;
    rDY[3] = dx4*60.0;
    rDY[4] = dx5*60.0;
    rDY[5] = dx6*60.0;
}

void TysonNovak2001OdeSystem::AnalyticJacobian(const std::vector<double>& rSolutionGuess, double** jacobian, double time, double timeStep)
{
    timeStep *= 60.0; // to scale Jacobian so in hours not minutes
    double x1 = rSolutionGuess[0];
    double x2 = rSolutionGuess[1];
    double x3 = rSolutionGuess[2];
    double x4 = rSolutionGuess[3];
    double x5 = rSolutionGuess[4];
    double x6 = rSolutionGuess[5];

    // f1
    double df1_dx1 = -mK2d - mK2dd*x2;
    double df1_dx2 = -mK2dd*x1;

    jacobian[0][0] =  1-timeStep*df1_dx1;
    jacobian[0][1] = -timeStep*df1_dx2;

    // f2
    double df2_dx1 = -mK4*x6*x2/(mJ4+x2);
    double df2_dx2 = -mJ3*(mK3d + mK3dd*x4)/(SmallPow((mJ3 + 1 - x2),2))
                     -mJ4*mK4*x6*x1/(SmallPow((mJ4+x2),2));
    double df2_dx4 =  mK3dd*(1-x2)/(mJ3+1-x2);
    double df2_dx6 = -mK4*x1*x2/(mJ4+x2);

    jacobian[1][0] = -timeStep*df2_dx1;
    jacobian[1][1] =  1-timeStep*df2_dx2;
    jacobian[1][3] = -timeStep*df2_dx4;
    jacobian[1][5] = -timeStep*df2_dx6;

    //f3
    double z = x1*x6/mJ5;
    double df3_dx1 = (mK5dd*x6/mJ5)*mN*SmallPow(z,mN-1)/(SmallPow((1-SmallPow(z,mN)),2));
    double df3_dx3 = -mK6;
    double df3_dx6 = (mK5dd*x1/mJ5)*mN*SmallPow(z,mN-1)/(SmallPow((1-SmallPow(z,mN)),2));

    jacobian[2][0] = -timeStep*df3_dx1;
    jacobian[2][2] = 1-timeStep*df3_dx3;
    jacobian[2][5] = -timeStep*df3_dx6;

    // f4
    double df4_dx3 =  mJ7*mK7*x5/(SmallPow(mJ7+x3-x4,2));
    double df4_dx4 = -mJ7*mK7*x5/(SmallPow(mJ7+x3-x4,2)) - mK6 - mJ8*mK8*mMad/(SmallPow(mJ8+x4,2));
    double df4_dx5 =  mK7*(x3-x4)/(mJ7+x3-x4);

    jacobian[3][2] = -timeStep*df4_dx3;
    jacobian[3][3] = 1-timeStep*df4_dx4;
    jacobian[3][4] = -timeStep*df4_dx5;

    // f5
    double df5_dx1 =  mK9*x6*(1-x5);
    double df5_dx5 = -mK10 - mK9*x6*x1;
    double df5_dx6 =  mK9*x1*(1-x5);

    jacobian[4][0] = -timeStep*df5_dx1;
    jacobian[4][4] = 1-timeStep*df5_dx5;
    jacobian[4][5] = -timeStep*df5_dx6;

    // f6
    double df6_dx6 = mMu - 2*mMu*x6/mMstar;

    jacobian[5][5] = 1-timeStep*df6_dx6;
}

bool TysonNovak2001OdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    // Only call this a stopping condition if the mass of the cell is over 0.6
    // (normally cycles from 0.5-1.0 ish!)
    return ( (rY[5] > 0.6 ) && (rY[0] < mCycB_threshold) && dy[0] < 0.0 );
}

double TysonNovak2001OdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    // Only call this a stopping condition if the mass of the cell is over 0.6
    // (normally cycles from 0.5-1.0 ish!)
    if (rY[5]<0.6)
    {
        return 1.0;
    }

    if (dy[0] >= 0.0)
    {
        return 1.0;
    }
    return rY[0]-mCycB_threshold;
}

template<>
void OdeSystemInformation<TysonNovak2001OdeSystem>::Initialise()
{
    /*
     * Initialise state variables.
     *
     * These initial conditions are the approximate steady state
     * solution values while the commented out conditions are taken
     * from the Tyson and Novak 2001 paper.
     */
    this->mVariableNames.push_back("CycB");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(0.1);
    this->mInitialConditions.push_back(0.099999999999977);

    this->mVariableNames.push_back("Cdh1");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(9.8770e-01);
    this->mInitialConditions.push_back(0.989026454281841);

    this->mVariableNames.push_back("Cdc20T");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(1.5011e+00);
    this->mInitialConditions.push_back(1.547942029285891);

    this->mVariableNames.push_back("Cdc20A");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(1.2924e+00);
    this->mInitialConditions.push_back(1.421110920155839);

    this->mVariableNames.push_back("IEP");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(6.5405e-01);
    this->mInitialConditions.push_back(0.672838844290094);

    this->mVariableNames.push_back("mass");
    this->mVariableUnits.push_back("");
//    this->mInitialConditions.push_back(4.7039e-01);
    this->mInitialConditions.push_back(0.970831277863956 / 2);

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001OdeSystem)
