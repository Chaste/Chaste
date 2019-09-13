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


#ifndef _ODE2_HPP_
#define _ODE2_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

/**
 * dy/dt = y*t, y(0) = 4.
 */
class Ode2 : public AbstractOdeSystem
{
public:

    Ode2() : AbstractOdeSystem(1) // 1 here is the number of unknowns
    {
        mpSystemInfo = OdeSystemInformation<Ode2>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0]*time;
    }
};

template<>
void OdeSystemInformation<Ode2>::Initialise()
{
    // These two lines are commented out, to show that variable names & units can
    // be left unspecified, provided they are not actually needed.  This also
    // ensures coverage of OdeSolution::WriteToFile, which contains a case for when
    // variable names are not given.
//    this->mVariableNames.push_back("Variable_1");
//    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(4.0);

    this->mInitialised = true;
}

#endif //_ODE2_HPP_
