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

//Heun's Method
//Megan Lewis, University of Saskatchewan
//January 2009
#include "HeunIvpOdeSolver.hpp"

void HeunIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& rCurrentYValues,
                                                  std::vector<double>& rNextYValues)
{
    /*
     * Apply Heun 2nd order method for each time step in AbstractOneStepIvpSolver.
     * Calculates a vector containing the next Y value from the current one for each
     * equation in the system.
     */

    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    std::vector<double> k1(num_equations);
    std::vector<double> k2(num_equations);
    std::vector<double>& dy = rNextYValues; // re-use memory

    // Work out k1
    pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, k1);

    // Add current y values to k1 values
    for (unsigned i=0; i<num_equations; i++)
    {
        dy[i] = timeStep*k1[i] + rCurrentYValues[i];
    }
    //Work out k2
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, dy, k2);

    // New solution
    for (unsigned i=0; i<num_equations; i++)
    {
        rNextYValues[i] = rCurrentYValues[i] + timeStep*0.5*(k1[i] + k2[i]);
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HeunIvpOdeSolver)
