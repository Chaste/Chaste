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

#ifndef _ODESYSTEMFORCOUPLEDHEATEQUATIONWITHSOURCE_HPP_
#define _ODESYSTEMFORCOUPLEDHEATEQUATIONWITHSOURCE_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "OdeSystemInformation.hpp"

/**
 * An ODE system given by
 *
 * dv/dt = a*v, v(0) = 1,
 *
 * that 'couples' with (although in this case, is actually independent of) the heat equation with source term
 *
 * u_t = del^2 u + v,
 *
 * which is defined in the separate class HeatEquationWithSourceForCoupledOdeSystem.
 */
class OdeSystemForCoupledHeatEquationWithSource : public AbstractOdeSystemForCoupledPdeSystem
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
    OdeSystemForCoupledHeatEquationWithSource(double a)
        : AbstractOdeSystemForCoupledPdeSystem(1,1),
          mA(a)
    {
        mpSystemInfo = OdeSystemInformation<OdeSystemForCoupledHeatEquationWithSource>::Instance();
        SetStateVariables(GetInitialConditions());
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
        rDY[0] = -mA * rY[0];
    }

    /**
     * @return mA
     */
    double GetA()
    {
        return mA;
    }
};

template<>
void OdeSystemInformation<OdeSystemForCoupledHeatEquationWithSource>::Initialise()
{
    this->mVariableNames.push_back("v");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);
    this->mInitialised = true;
}

#endif //_ODESYSTEMFORCOUPLEDHEATEQUATIONWITHSOURCE_HPP_
