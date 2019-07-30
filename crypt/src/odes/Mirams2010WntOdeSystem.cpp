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
#include "Mirams2010WntOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

// These #includes are needed for the constructor and EvaluateYDerivatives()
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"

Mirams2010WntOdeSystem::Mirams2010WntOdeSystem(double wntLevel,
                                               boost::shared_ptr<AbstractCellMutationState> pMutationState,
                                               std::vector<double> stateVariables)
    : AbstractOdeSystem(3),
      mpMutationState(pMutationState),
      mWntLevel(wntLevel)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Mirams2010WntOdeSystem>);

    /**
     * State variables.
     *
     * 0. b1 = Beta-Catenin (1st allele's copy)
     * 1. b2 = Beta-Catenin (2nd allele's copy)
     * 2. wntLevel
     */

    Init(); // set up parameter values

    // Set up rough guesses for the initial steady states in this Wnt conc.
    double b1 = 0;
    double b2 = 0;
    b1 = (mA/2.0) / (((wntLevel + mB)/(mC*wntLevel + mD)) + mF);
    if (!mpMutationState)
    {
        // No mutations specified
    }
    else if (mpMutationState->IsType<BetaCateninOneHitCellMutationState>())
    {
        b2 = (mA/2.0)/mF;
    }
    else
    {
        b2 = (mA/2.0) / (((wntLevel + mB)/(mC*wntLevel + mD)) + mF);
    }

    SetDefaultInitialCondition(0, b1);
    SetDefaultInitialCondition(1, b2);
    SetDefaultInitialCondition(2, wntLevel);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

void Mirams2010WntOdeSystem::SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState)
{
    mpMutationState = pMutationState;
}

Mirams2010WntOdeSystem::~Mirams2010WntOdeSystem()
{
    // Do nothing
}

void Mirams2010WntOdeSystem::Init()
{
    // Initialise model parameter values
    mA = 25.38;     // nM hr^{-1}
    mB = 0.1;       // dimensionless
    mC = 6.386;     // hr
    mD = 9.818e-2;  // hr
    mE = 1.2e3;     // nM
    mF = 1.54e-2;   // hr^{-1}
}

void Mirams2010WntOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double x1 = rY[0];
    double x2 = rY[1];
    double wnt_level = rY[2];

    double dx1 = 0.0;
    double dx2 = 0.0;
    /*
     * The variables are
     * 1. b = Beta-Catenin1
     * 2. b = Beta-Catenin2
    */

    double c = mC;
    double d = mD;
    // Mutations take effect by altering the parameters
    if (mpMutationState->IsType<ApcOneHitCellMutationState>()) // APC +/-
    {
        c = 31.87;
        d = 0.490;
    }
    else if (mpMutationState->IsType<ApcTwoHitCellMutationState>()) // APC -/-
    {
        c = 71.21;
        d = 1.095;
    }

    // da
    dx1 = mA/2.0 - (((wnt_level + mB)/(c*wnt_level + d))*(mE/(mE+x1)) + mF)*x1;
    // db
    if (mpMutationState->IsType<BetaCateninOneHitCellMutationState>())
    {
        dx2 = mA/2.0 - mF*x2;
    }
    else
    {
        dx2 = mA/2.0 - (((wnt_level + mB)/(c*wnt_level + d))*(mE/(mE+x2)) + mF)*x2;
    }

    rDY[0] = dx1;
    rDY[1] = dx2;
    rDY[2] = 0.0; // do not change the Wnt level
}

const boost::shared_ptr<AbstractCellMutationState> Mirams2010WntOdeSystem::GetMutationState() const
{
    return mpMutationState;
}

template<>
void CellwiseOdeSystemInformation<Mirams2010WntOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Beta_Cat_Allele1");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(50.0); // will be filled in later

    this->mVariableNames.push_back("Beta_Cat_Allele2");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(50.0); // will be filled in later

    this->mVariableNames.push_back("Wnt Level");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    this->mInitialised = true;
}

double Mirams2010WntOdeSystem::GetWntLevel() const
{
    return mWntLevel;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Mirams2010WntOdeSystem)
