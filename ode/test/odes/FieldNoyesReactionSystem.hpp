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


/**
 * Concrete FieldNoyesReactionSystem class
 */
#ifndef FIELDNOYESREACTIONSYSTEM_HPP_
#define FIELDNOYESREACTIONSYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

class FieldNoyesReactionSystem : public AbstractOdeSystem
{
public:
    FieldNoyesReactionSystem() : AbstractOdeSystem(3)
    {
        mpSystemInfo = OdeSystemInformation<FieldNoyesReactionSystem>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        const double epsilon = 0.05;
        const double p = 6.7;
        const double q = 1e-4;
        const double f = 0.5;

        rDY[0] = (rY[1] - rY[0]*rY[1] + rY[0]*(1 - q*rY[0]))/epsilon;
        rDY[1] = -rY[1] - rY[0]*rY[1] + 2*f*rY[2];
        rDY[2] = (rY[0] - rY[2])/p;
    }
};

template<>
void OdeSystemInformation<FieldNoyesReactionSystem>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Variable_2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Variable_3");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}

#endif /*FIELDNOYESREACTIONSYSTEM_HPP_*/
