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


#ifndef _ODE5JACOBIAN_HPP_
#define _ODE5JACONIAN_HPP_

#include "AbstractOdeSystemWithAnalyticJacobian.hpp"
#include "OdeSystemInformation.hpp"

/**
 * dy/dt = 100 y (1-y), y(0) = 0.2.
 */
class Ode5Jacobian : public AbstractOdeSystemWithAnalyticJacobian
{
private:

    double mAlpha;

public:

    Ode5Jacobian() : AbstractOdeSystemWithAnalyticJacobian(1)  // 1 here is the number of unknowns
    {
        mpSystemInfo = OdeSystemInformation<Ode5Jacobian>::Instance();
        mAlpha = 100;
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = mAlpha*rY[0]*(1-rY[0]);
    }

    void AnalyticJacobian(const std::vector<double>& rSolutionGuess, double** jacobian, double time, double timeStep)
    {
        jacobian[0][0] = 1.0 - timeStep*mAlpha + 2.0*timeStep*mAlpha*rSolutionGuess[0];
    }
};

template<>
void OdeSystemInformation<Ode5Jacobian>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.2);

    this->mInitialised = true;
}

#endif /*_ODE5JACOBIAN_HPP_*/
