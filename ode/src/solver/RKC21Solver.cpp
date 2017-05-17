/*

Copyright (c) 2005-2016, University of Oxford.
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


/**
 * A concrete one step ODE solver class that employs the Runge Kutta Chebyshev 1st order
 * 2 stages solver.
 * Wenxian Guo and Raymond Spiteri
 * University of Saskatchewan
 * May 2017
 */

#include "RKC21Solver.hpp"

void RKC21Solver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& rCurrentYValues,
                                                  std::vector<double>& rNextYValues)
{
    /*
     * Apply Runge-Kutta-Chebyshev 1st order 2 stages method for each timestep in AbstractOneStepIvpSolver.
     * Calculates a vector containing the next Y value from the current one for each
     * equation in the system.
     */


     /*
     RKC scheme:
     w1 = w0+mu1_tilde*dt*F(tn,w0)
     w2 = (1-muj-nuj)*w0 + muj*w1+nuj*w0 + muj_tilde*dt*F(tn,w1) + gammaj_tilde*dt*F(tn,w0)
     where F(nj) = F(tn+cj *dt, wj) cj = j^2/s^2
     A decent reference on this method refer to Hundsdorfer, Willem, and Jan G. Verwer. Numerical solution of time-dependent advection-diffusion-reaction equations. Vol. 33. Springer Science & Business Media, 2013. P429
     */

     //RKC coefficients
    const double mu1_tilde = 0.256134735558604;
    const double mu2 = 1.952097590002976;
    //const double nu2 = -b2/b0; canceled out in equation
    const double mu2_tilde =0.500000000000000;
    //const double gamma2_tilde = 0; never used in first order RKC

    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();

    std::vector<double>& w0 = rCurrentYValues; // alias
    std::vector<double> w1(num_equations);
    std::vector<double>& w2 = rNextYValues; // alias
    std::vector<double> F0(num_equations);
    std::vector<double>& F1 = rNextYValues;

    // Work out w1
    pAbstractOdeSystem->EvaluateYDerivatives(time, w0, F0);

    for (unsigned i=0; i<num_equations; i++)
    {
        w1[i] = w0[i] + mu1_tilde*timeStep*F0[i];
    } 

    // Work next step
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep/4.0, w1, F1);
    for (unsigned i=0; i<num_equations; i++)
    {
        w2[i] = (1-mu2)*w0[i] + mu2*w1[i] + mu2_tilde*timeStep*F1[i];// + gamma2_tilde*timeStep*F0[i];
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RKC1Order2StagesSolver)
