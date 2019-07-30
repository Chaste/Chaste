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

/*
Megan E. Marsh, Raymond J. Spiteri
Numerical Simulation Laboratory
University of Saskatchewan
December 2011
Partial support provided by research grants from the National
Science and Engineering Research Council (NSERC) of Canada
and the MITACS/Mprime Canadian Network of Centres of Excellence.
*/
#include <cmath>
#include "GRL1IvpOdeSolver.hpp"

void GRL1IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                           double timeStep,
                                           double time,
                                           std::vector<double>& rCurrentYValues,
                                           std::vector<double>& rNextYValues)
{
    /*
     * Apply GRL1 first-order method for each time step in AbstractOneStepIvpSolver.
     * Calculates a vector containing the next Y value from the current one for each
     * equation in the system.
     */
    const double delta = 1.0e-8; // The step for numerical Jacobian calculation

    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    if (mEvalF.size() != num_equations)
    {
        mEvalF.resize(num_equations);
        mPartialF.resize(num_equations);
        mTemp.resize(num_equations);
    }

    double tempY;

    pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, mEvalF);
    for (unsigned i=0; i<num_equations; i++)
    {
        tempY= rCurrentYValues[i];
        rCurrentYValues[i]=tempY+delta;
        pAbstractOdeSystem->EvaluateYDerivatives(time, rCurrentYValues, mTemp);
        mPartialF[i]=(mTemp[i]-mEvalF[i])/delta;
        rCurrentYValues[i]=tempY;
    }
    // New solution
    for (unsigned i=0; i<num_equations; i++)
    {
        // std::cout<<"Partial F"<< mPartialF[i]<<std::endl;
        if (fabs(mPartialF[i]) < delta)
        {
            rNextYValues[i] = rCurrentYValues[i]+mEvalF[i]*timeStep;
        }
        else
        {
            rNextYValues[i] = rCurrentYValues[i]+(mEvalF[i]/mPartialF[i])*(exp(mPartialF[i]*timeStep)-1);
        }
    }
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(GRL1IvpOdeSolver)
