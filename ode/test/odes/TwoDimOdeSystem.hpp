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


#ifndef TWODIMODESYSTEM_HPP_
#define TWODIMODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

class TwoDimOdeSystem : public AbstractOdeSystem
{
public:
    TwoDimOdeSystem() : AbstractOdeSystem(2)
    {
        mpSystemInfo = OdeSystemInformation<TwoDimOdeSystem>::Instance();
        mStateVariables.push_back(3);
        mStateVariables.push_back(4);
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY.assign(rY.begin(), rY.end());
    }
};

template<>
void OdeSystemInformation<TwoDimOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1);

    this->mVariableNames.push_back("Variable_2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(2);

    this->mInitialised = true;
}


#endif /*TWODIMODESYSTEM_HPP_*/
