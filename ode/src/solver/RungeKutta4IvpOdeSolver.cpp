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

#include "RungeKutta4IvpOdeSolver.hpp"

void RungeKutta4IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& rCurrentYValues,
                                                  std::vector<double>& rNextYValues)
{
    /*
     * Apply Runge-Kutta 4th order method for each timestep in AbstractOneStepIvpSolver.
     * Calculates a vector containing the next Y value from the current one for each
     * equation in the system.
     */

    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    if (num_equations != k1.size())
    {
        k1.resize(num_equations);
        k2.resize(num_equations);
        k3.resize(num_equations);
        k4.resize(num_equations);
        yki.resize(num_equations);
    }

    std::vector<double>& dy = rNextYValues; // re-use memory

    pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k1[i] = timeStep*dy[i];
        yki[i] = rCurrentYValues[i] + 0.5*k1[i];
    }

    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, yki, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k2[i] = timeStep*dy[i];
        yki[i] = rCurrentYValues[i] + 0.5*k2[i];
    }

    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, yki, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k3[i] = timeStep*dy[i];
        yki[i] = rCurrentYValues[i] + k3[i];
    }

    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, yki, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        k4[i] = timeStep*dy[i];
        rNextYValues[i] = rCurrentYValues[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RungeKutta4IvpOdeSolver)
