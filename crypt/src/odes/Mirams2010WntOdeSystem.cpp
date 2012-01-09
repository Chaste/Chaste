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
#include "Mirams2010WntOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

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
