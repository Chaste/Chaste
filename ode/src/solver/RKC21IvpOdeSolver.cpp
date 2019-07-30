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


/**
 * A concrete one step ODE solver class that employs the Runge Kutta Chebyshev 1st order
 * 2 stages solver.
 * Wenxian Guo and Raymond Spiteri
 * University of Saskatchewan
 * May 2017
 */

#include "RKC21IvpOdeSolver.hpp"

void RKC21IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
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
     * The General RKC scheme:

    * w_{n0} = w_n;
    * w_{n1} = w_n + \tilde{\mu_1} * dt * F_{n0};
    * w_{nj} = (1 - mu_j - nu_j) w_n + mu_j*w_{n,j-1} + nu_j*w_{n,j-2} + \tilde{\mu_j}*dt*F_{n,j-1} + \tilde{\gamma_j}*dt*F_{n,0},
    * w_{n+1} = w_{ns}.
    * where:
    *       j = 2,...,s;
    *       F_{nk} = F(t_n + c_k*dt, w_{nk});
    *       c_k = T_s(w_0)*Tp_k(w_0) / (Tp_s(w_0)*T_k(w_0)) with T_k/Tp_k are k-th Chebyshev polynomials/ derivatives of them
    *       w_0 is specific parameters
    * In RKC21, c_2 = \tilde{\mu_1}


     A decent reference on this method refer to Hundsdorfer, Willem, and Jan G. Verwer. Numerical solution of time-dependent advection-diffusion-reaction equations. Vol. 33. Springer Science & Business Media, 2013. P429
     */

     //RKC coefficients, hard-coded //TODO how to compute/where to store arbitrary RKC coefficients?

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
        w1[i] = w0[i] + mu1_tilde * timeStep * F0[i];
    }

    // Work next step
    pAbstractOdeSystem->EvaluateYDerivatives(time + mu1_tilde * timeStep, w1, F1);
    for (unsigned i=0; i<num_equations; i++)
    {
        w2[i] = (1-mu2) * w0[i] + mu2 * w1[i] + mu2_tilde * timeStep * F1[i];// + gamma2_tilde*timeStep*F0[i]; -> never used
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RKC21IvpOdeSolver)
