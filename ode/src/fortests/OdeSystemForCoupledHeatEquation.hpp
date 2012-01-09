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

#ifndef _ODESYSTEMFORCOUPLEDHEATEQUATION_HPP_
#define _ODESYSTEMFORCOUPLEDHEATEQUATION_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "OdeSystemInformation.hpp"

/**
 * An ODE system given by
 *
 * dv/dt = a*u, v(0) = 1,
 *
 * that couples with the heat equation
 *
 * u_t = div (grad u),
 *
 * which is defined in the separate class HeatEquationForCoupledOdeSystem.
 */
class OdeSystemForCoupledHeatEquation : public AbstractOdeSystemForCoupledPdeSystem
{
private:

    /** A parameter for use in the ODE system. */
    double mA;

public:

    /**
     * Constructor.
     *
     * @param a the value of the parameter mA
     */
    OdeSystemForCoupledHeatEquation(double a)
        : AbstractOdeSystemForCoupledPdeSystem(1,1),
          mA(a)
    {
        mpSystemInfo = OdeSystemInformation<OdeSystemForCoupledHeatEquation>::Instance();
        ResetToInitialConditions();
    }

    /**
     * Method to evaluate the derivatives of the system.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     * @param rDY  storage for the derivatives of the system; will be filled in on return
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        assert(mPdeSolutionSize == 1);
        double u = mPdeSolution[0];
        rDY[0] = mA * u;
    }
};

template<>
void OdeSystemInformation<OdeSystemForCoupledHeatEquation>::Initialise()
{
    this->mVariableNames.push_back("v");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);
    this->mInitialised = true;
}

#endif //_ODESYSTEMFORCOUPLEDHEATEQUATION_HPP_
