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
