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

#include "ParameterisedOde.hpp"
#include "OdeSystemInformation.hpp"

bool ParameterisedOde::fakeSecondParameter = false;
bool ParameterisedOde::noParameterDefaults = false;

ParameterisedOde::ParameterisedOde() : AbstractOdeSystem(1) // 1 here is the number of variables
{
    mpSystemInfo = OdeSystemInformation<ParameterisedOde>::Instance();
    ResetToInitialConditions();
    if (!noParameterDefaults)
    {
        mParameters.push_back(0);
        if (fakeSecondParameter)
        {
            mParameters.push_back(-1.0);
        }
    }
}

void ParameterisedOde::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    rDY[0] = mParameters[0];
}

std::vector<double> ParameterisedOde::ComputeDerivedQuantities(double time,
                                                               const std::vector<double>& rState)
{
    std::vector<double> dqs;
    dqs.push_back(2*mParameters[0] + rState[0]);
    return dqs;
}


#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ParameterisedOde);


template<>
void OdeSystemInformation<ParameterisedOde>::Initialise()
{
    this->mSystemName = "ParameterisedOde";

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
