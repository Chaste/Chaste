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

#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

Alarcon2004OxygenBasedCellCycleOdeSystem::Alarcon2004OxygenBasedCellCycleOdeSystem(double oxygenConcentration,
                                                                                   bool isLabelled,
                                                                                   std::vector<double> stateVariables)
    : AbstractOdeSystem(6),
      mOxygenConcentration(oxygenConcentration),
      mIsLabelled(isLabelled)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Alarcon2004OxygenBasedCellCycleOdeSystem>);

    /**
     * State variables
     *
     * 0. x = Cdh1-APC complexes
     * 1. y = cyclin-CDK
     * 2. z = p27
     * 3. m = mass
     * 4. u = RBNP
     * 5. oxygenConcentration
     */
    Init(); // set up parameters

    // Parameter values are taken from the Alarcon et al. (2004) paper
    if (mIsLabelled) // labelled "cancer" cells
    {
        ma1 = 0.04;
        mc1 = 0.007;
        mxThreshold = 0.004;
        myThreshold = 0.05;
    }
    else // normal cells
    {
        ma1 = 0.05;
        mc1 = 0.1;
        mxThreshold = 0.004;
        myThreshold = 0.2;
    }

    // Cell-specific initial conditions
    SetDefaultInitialCondition(3, 0.5*mMstar);
    SetDefaultInitialCondition(5, oxygenConcentration);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Alarcon2004OxygenBasedCellCycleOdeSystem::~Alarcon2004OxygenBasedCellCycleOdeSystem()
{
    // Do nothing
}

void Alarcon2004OxygenBasedCellCycleOdeSystem::Init()
{
    // Parameter values are taken from the Alarcon et al. (2004) paper
    ma2 = 1.0;
    ma3 = 0.25;
    ma4 = 0.04;
    mb3 = 10.0;
    mb4 = 5.5;
    mc2 = 0.01;
    md1 = 0.01;
    md2 = 0.1;
    mJ3 = 0.04;
    mJ4 = 0.04;
    mEta = 0.01;
    mMstar = 10.0;
    mB = 0.01;
}

void Alarcon2004OxygenBasedCellCycleOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double x = rY[0];
    double y = rY[1];
    double z = rY[2];
    double mass = rY[3];
    double u = rY[4];
    double oxygen_concentration = rY[5];

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double dmass = 0.0;
    double du = 0.0;

    /*
     % The variables are
     % 1. x = Cdh1-APC complexes
     % 2. y = cyclin-CDK
     % 3. z = p27
     % 4. m = mass
     % 5. u = RBNP
    */

    dx = ((1 + mb3*u)*(1-x))/(mJ3 + 1 - x) - (mb4*mass*x*y)/(mJ4 + x);
    dy = ma4 -(ma1 + ma2*x + ma3*z)*y;

    // Parameter values are taken from the Alarcon et al. (2004) paper
    if (mIsLabelled) // labelled "cancer" cells
    {
        dz = mc1 - mc2*oxygen_concentration*z/(mB + oxygen_concentration);
    }
    else // normal cells
    {
        dz = mc1*(1 - mass/mMstar) - mc2*oxygen_concentration*z/(mB + oxygen_concentration);
    }

    dmass = mEta*mass*(1 - mass/mMstar);
    du = md1 - (md2 + md1*y)*u;

    // Rescale time to be in hours
    rDY[0] = 60.0*dx;
    rDY[1] = 60.0*dy;
    rDY[2] = 60.0*dz;
    rDY[3] = 60.0*dmass;
    rDY[4] = 60.0*du;
    rDY[5] = 0.0; // do not change the oxygen concentration
}

bool Alarcon2004OxygenBasedCellCycleOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    return (rY[0] < mxThreshold && rY[1] > myThreshold);
}

template<>
void CellwiseOdeSystemInformation<Alarcon2004OxygenBasedCellCycleOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Cdh1_APC_complexes");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(0.9);

    this->mVariableNames.push_back("cyclin_CDK");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("p27");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("mass");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("RBNP");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("O2");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mInitialised = true;
}

void Alarcon2004OxygenBasedCellCycleOdeSystem::SetIsLabelled(bool isLabelled)
{
    mIsLabelled = isLabelled;
}

bool Alarcon2004OxygenBasedCellCycleOdeSystem::IsLabelled() const
{
    return mIsLabelled;
}

double Alarcon2004OxygenBasedCellCycleOdeSystem::GetOxygenConcentration() const
{
    return mOxygenConcentration;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleOdeSystem)
