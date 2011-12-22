/*

Copyright (C) University of Oxford, 2005-2011

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



#ifndef _ODETHIRDORDER_HPP
#define _ODETHIRDORDER_HPP

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

/**
 * Concrete OdeThirdOrder class
 */
class OdeThirdOrder : public AbstractOdeSystem
{
public:
    OdeThirdOrder()
            : AbstractOdeSystem(3) // 3 here is the number of unknowns
    {
        mpSystemInfo = OdeSystemInformation<OdeThirdOrder>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0]-rY[1]+rY[2];
        rDY[1] = rY[1]-rY[2];
        rDY[2] = 2*rY[1]-rY[2];
    }
};

template<>
void OdeSystemInformation<OdeThirdOrder>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Variable_2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Variable_3");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

#endif //_ODETHIRDORDER_HPP
