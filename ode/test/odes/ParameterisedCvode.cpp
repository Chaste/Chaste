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

#include "ParameterisedCvode.hpp"
#include "OdeSystemInformation.hpp"
#include "VectorHelperFunctions.hpp"

#ifdef CHASTE_CVODE
bool ParameterisedCvode::fakeSecondParameter = false;
bool ParameterisedCvode::noParameterDefaults = false;

ParameterisedCvode::ParameterisedCvode() : AbstractCvodeSystem(1) // 1 here is the number of variables
{
    this->mpSystemInfo = OdeSystemInformation<ParameterisedCvode>::Instance();
    Init();
    if (!noParameterDefaults)
    {
        CreateVectorIfEmpty(mParameters,1);
        SetVectorComponent(mParameters, 0, 0.0);
        if (fakeSecondParameter)
        {
            SetVectorComponent(mParameters, 1, -1.0);
        }
    }
}

void ParameterisedCvode::EvaluateYDerivatives(double time, const N_Vector rY, N_Vector rDY)
{
    NV_Ith_S(rDY, 0) = GetVectorComponent(mParameters,0);
}

N_Vector ParameterisedCvode::ComputeDerivedQuantities(double time,
                                                      const N_Vector& rState)
{
    N_Vector derived_quantities = NULL;
    CreateVectorIfEmpty(derived_quantities,1);
    SetVectorComponent(derived_quantities, 0,
                       2*GetVectorComponent(mParameters,0) + GetVectorComponent(rState,0));
    return derived_quantities;
}

template<>
void OdeSystemInformation<ParameterisedCvode>::Initialise()
{
    this->mSystemName = "ParameterisedCvode";

    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mParameterNames.push_back("a");
    this->mParameterUnits.push_back("dimensionless");

    this->mDerivedQuantityNames.push_back("2a_plus_y");
    this->mDerivedQuantityUnits.push_back("dimensionless");

    this->mAttributes["attr"] = 1.1;

    this->mInitialised = true;
}
#endif // CHASTE_CVODE
