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

#include "WntCellCycleOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"
#include "IsNan.hpp"

// These #includes are needed for the constructor and EvaluateYDerivatives()
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

WntCellCycleOdeSystem::WntCellCycleOdeSystem(double wntLevel,
                                             boost::shared_ptr<AbstractCellMutationState> pMutationState,
                                             std::vector<double> stateVariables)
    : AbstractOdeSystem(9),
      mpMutationState(pMutationState),
      mWntLevel(wntLevel)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<WntCellCycleOdeSystem>);

    /**
     * State variables.
     *
     * 0. r = pRb
     * 1. e = E2F1 (This is the S-phase indicator)
     * 2. i = CycD (inactive)
     * 3. j = CycD (active)
     * 4. p = pRb-p
     * 5. c = APC (Active)
     * 6. b1 = Beta-Catenin (1st allele's copy)
     * 7. b2 = Beta-Catenin (2nd allele's copy)
     * 8. wntLevel
     */

    Init(); // set up parameter values

    // Set up a Wnt signalling pathway in a steady state
    double destruction_level = ma5d/(ma4d*wntLevel+ma5d);
    double beta_cat_level_1 = -1.0;
    double beta_cat_level_2 = -1.0;

    if (!mpMutationState)
    {
        // No mutations specified
    }
    else if (mpMutationState && mpMutationState->IsType<ApcOneHitCellMutationState>())
    {
        // APC +/- : only half are active
        beta_cat_level_1 = 0.5*ma2d/(ma2d+0.5*ma3d*destruction_level);
        beta_cat_level_2 = 0.5*ma2d/(ma2d+0.5*ma3d*destruction_level);
    }
    else if (mpMutationState && mpMutationState->IsType<ApcTwoHitCellMutationState>())
    {
        // APC -/-
        destruction_level = 0.0; // no active destruction complex
        beta_cat_level_1 = 0.5; // fully active beta-catenin
        beta_cat_level_2 = 0.5; // fully active beta-catenin
    }
    else if (mpMutationState && mpMutationState->IsType<BetaCateninOneHitCellMutationState>())
    {
        // Beta-cat delta 45
        beta_cat_level_1 = 0.5*ma2d/(ma2d+ma3d*destruction_level);
        beta_cat_level_2 = 0.5;
    }
    else
    {
        // healthy cells
        beta_cat_level_1 = 0.5*ma2d/(ma2d+ma3d*destruction_level);
        beta_cat_level_2 = 0.5*ma2d/(ma2d+ma3d*destruction_level);
    }

    // Cell-specific initial conditions
    SetDefaultInitialCondition(5, destruction_level);
    SetDefaultInitialCondition(6, beta_cat_level_1);
    SetDefaultInitialCondition(7, beta_cat_level_2);
    SetDefaultInitialCondition(8, wntLevel);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

void WntCellCycleOdeSystem::SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState)
{
    mpMutationState = pMutationState;
}

WntCellCycleOdeSystem::~WntCellCycleOdeSystem()
{
    // Do nothing
}

void WntCellCycleOdeSystem::Init()
{
    // Initialise model parameter values
    // Swat (2004) Parameters
    double k1 = 1.0;
    double k2 = 1.6;
    double k3 = 0.05;
    double k16 = 0.4;
    double k34 = 0.04;
    double k43 = 0.01;
    double k61 = 0.3;
    double k23 = 0.3;
    double a = 0.04;
    double J11 = 0.5;
    double J12 = 5.0;
    double J61 = 5.0;
    double J62 = 8.0;
    double J13 = 0.002;
    double J63 = 2.0;
    double Km1 = 0.5;
    double Km2 = 4.0;
    double Km4 = 0.3;
    double kp = 0.05;
    double phi_pRb = 0.005;
    double phi_E2F1 = 0.1;
    double phi_CycDi = 0.023;
    double phi_CycDa = 0.03;
    double phi_pRbp = 0.06;

    // Mirams et al. parameter values
    double a1 = 0.423;
    double a2 = 2.57e-4;
    double a3 = 1.72;
    double a4 = 10.0;
    double a5 = 0.5;
    double WntMax = 10.0;
    double mitogenic_factorF = 6.0e-4;
    double APC_Total = 0.02;

    // Non-dimensionalise...
    mk2d = k2/(Km2*phi_E2F1);
    mk3d = k3*a1*mitogenic_factorF/(Km4*phi_E2F1*a2);
    mk34d = k34/phi_E2F1;
    mk43d = k43/phi_E2F1;
    mk23d = k23*Km2/(Km4*phi_E2F1);
    mad = a/Km2;
    mJ11d = J11*phi_E2F1/k1;
    mJ12d = J12*phi_E2F1/k1;
    mJ13d = J13*phi_E2F1/k1;
    mJ61d = J61*phi_E2F1/k1;
    mJ62d = J62*phi_E2F1/k1;
    mJ63d = J63*phi_E2F1/k1;
    mKm1d = Km1/Km2;
    mkpd = kp/(Km2*phi_E2F1);
    mphi_r = phi_pRb/phi_E2F1;
    mphi_i = phi_CycDi/phi_E2F1;
    mphi_j = phi_CycDa/phi_E2F1;
    mphi_p = phi_pRbp/phi_E2F1;
    ma2d = a2/phi_E2F1;
    ma3d = a3*APC_Total/phi_E2F1;
    ma4d = a4*WntMax/phi_E2F1;
    ma5d = a5/phi_E2F1;
    mk16d = k16*Km4/phi_E2F1;
    mk61d = k61/phi_E2F1;
    mPhiE2F1 = phi_E2F1;
}

void WntCellCycleOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double r = rY[0];
    double e = rY[1];
    double i = rY[2];
    double j = rY[3];
    double p = rY[4];
    double c = rY[5];
    double b1 = rY[6];
    double b2 = rY[7];
    double wnt_level = rY[8];

    double dx1 = 0.0;
    double dx2 = 0.0;
    double dx3 = 0.0;
    double dx4 = 0.0;
    double dx5 = 0.0;
    double dx6 = 0.0;
    double dx7 = 0.0;
    double dx8 = 0.0;

    /*
     * The variables are
     * 1. r = pRb
     * 2. e = E2F1
     * 3. i = CycD (inactive)
     * 4. j = CycD (active)
     * 5. p = pRb-p
     * 6. c = APC (Active)
     * 7. b = Beta-Catenin
    */

    // Bit back-to-front, but work out the Wnt section first...

    // Mutations take effect by altering the level of beta-catenin
    if (!mpMutationState)
    {
        // No mutations specified
    }
    else if (mpMutationState->IsType<ApcOneHitCellMutationState>()) // APC +/-
    {
        dx6 = ma5d*(1.0-c) - ma4d*wnt_level*c;
        dx7 = ma2d*(0.5-b1) - 0.5*ma3d*b1*c;
        dx8 = ma2d*(0.5-b2) - 0.5*ma3d*b2*c;
    }
    else if (mpMutationState->IsType<ApcTwoHitCellMutationState>()) // APC -/-
    {
        dx6 = 0.0;
        dx7 = ma2d*(0.5-b1);
        dx8 = ma2d*(0.5-b2);
    }
    else if (mpMutationState->IsType<BetaCateninOneHitCellMutationState>()) // Beta-Cat D45
    {
        dx6 = ma5d*(1.0-c) - ma4d*wnt_level*c;
        dx7 = ma2d*(0.5-b1) - ma3d*b1*c;
        dx8 = ma2d*(0.5-b2);
    }
    else
    {
        // da
        dx6 = ma5d*(1.0-c) - ma4d*wnt_level*c;
        // db
        dx7 = ma2d*(0.5-b1) - ma3d*b1*c;
        dx8 = ma2d*(0.5-b2) - ma3d*b2*c;
    }

    // Now the cell cycle stuff...

    // dr
    dx1 = e/(mKm1d+e)*mJ11d/(mJ11d+r)*mJ61d/(mJ61d+p) - mk16d*r*j+mk61d*p-mphi_r*r;
    // de
    dx2 = mkpd+mk2d*(mad*mad+e*e)/(1+e*e)*mJ12d/(mJ12d+r)*mJ62d/(mJ62d+p) - e;
    // di
    dx3 = mk3d*(b1+b2) + mk23d*e*mJ13d/(mJ13d+r)*mJ63d/(mJ63d+p) + mk43d*j - mk34d*i*j/(1+j) - mphi_i*i;
    // dj
    dx4 = mk34d*i*j/(1+j) - (mk43d+mphi_j)*j;
    // dp
    dx5 = mk16d*r*j - mk61d*p - mphi_p*p;

    double factor = mPhiE2F1*60.0;  // convert non-dimensional d/dt s to d/dt in hours

    rDY[0] = dx1*factor;
    rDY[1] = dx2*factor;
    rDY[2] = dx3*factor;
    rDY[3] = dx4*factor;
    rDY[4] = dx5*factor;
    rDY[5] = dx6*factor;
    rDY[6] = dx7*factor; // beta-cat allele 1
    rDY[7] = dx8*factor; // beta-cat allele 2
    rDY[8] = 0.0; // do not change the Wnt level
}

const boost::shared_ptr<AbstractCellMutationState> WntCellCycleOdeSystem::GetMutationState() const
{
    return mpMutationState;
}

bool WntCellCycleOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    double r = rY[0];
    double e = rY[1];
    double p = rY[4];
    double dY1 = mkpd+mk2d*(mad*mad+e*e)/(1+e*e)*mJ12d/(mJ12d+r)*mJ62d/(mJ62d+p) - e;
    double factor = mPhiE2F1*60.0;  // Convert non-dimensional d/dt s to d/dt in hours.
    dY1 = dY1*factor;

    assert(!std::isnan(rY[1]));
    assert(!std::isnan(dY1));
    return (rY[1] > 1.0 && dY1 > 0.0);
}

double WntCellCycleOdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    return rY[1] - 1.0;
}

template<>
void CellwiseOdeSystemInformation<WntCellCycleOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("pRb");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(7.357000000000000e-01);

    this->mVariableNames.push_back("E2F1");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(1.713000000000000e-01);

    this->mVariableNames.push_back("CycD_i");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(6.900000000000001e-02);

    this->mVariableNames.push_back("CycD_a");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(3.333333333333334e-03);

    this->mVariableNames.push_back("pRb_p");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(1.000000000000000e-04);

    this->mVariableNames.push_back("APC");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Beta_Cat1");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Beta_Cat2");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Wnt");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mInitialised = true;
}

double WntCellCycleOdeSystem::GetWntLevel() const
{
    return mWntLevel;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(WntCellCycleOdeSystem)
