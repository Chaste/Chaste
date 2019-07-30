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

#ifndef NHSMODELWITHBACKWARDSOLVER_HPP_
#define NHSMODELWITHBACKWARDSOLVER_HPP_

#include "NhsContractionModel.hpp"
#include "TimeStepper.hpp"
#include <cmath>


/**
 *
 *  The full backward Euler method on the NHS system compute the solution w^{n+1} in R^5
 *  by:
 *   w^{n+1}  - dt g(w^{n+1})  =  w^n
 *  where the ODE system is dw/dt = g(w).
 *
 *  This would involve solving a 5d nonlinear system at each timestep. However, only g1 and g2
 *  (ie dCatrop/dt and dz/dt) are nonlinear, g3, g4 and g5 (corresponding to dQi/dt) are linear
 *  and uncoupled. Therefore the backward euler solutions of Q1,Q2,Q3 can be computed immediately,
 *  leaving a 2D nonlinear system to be solved using newton's method
 */
class NhsModelWithBackwardSolver : public NhsContractionModel
{
private:
    /** Tolerance for solving nonlinear system which require newton's method */
    static const double mTolerance;

    /** Timestep for the ODEs solving */
    double mDt;

    /**
     *  Solve for Q1,Q2,Q3 (and therefore Q) implicitly using backward euler.
     *  These can be done directly as the rhs is linear in Qi
     *  @return Q=Q1+Q2+Q3
     */
    double ImplicitSolveForQ();

    /**
     *  The same as EvaluateYDerivatives in NhsContractionModel, but doesn't use std::vectors (for
     *  efficiency, and because this class is hardcoded and hand-optimised for backward euler),
     *  and just returns the derivatives of the first two components
     *  @param calciumTroponin Calcium troponin
     *  @param z z
     *  @param Q Q=Q1+Q2+Q3
     *  @param dCaTrop the returned value of dCaTrop/dt
     *  @param dz the returned value of dz/dt
     */
    void CalculateCaTropAndZDerivatives(double calciumTroponin, double z, double Q, double& dCaTrop, double& dz);

    /**
     *  Compute the residual function for the 2D nonlinear system when Backward Euler is used
     *  The backward Euler discretisation is: w^{n+1}  - dt g(w^{n+1}) =  w^n
     *  where: w = (ca_trop, z)
     *  and the ODE system is: dw/dt = g(w)
     *
     *  so the residual is: f = w^{n+1} - dt g(w^{n+1}) -  w^n
     *  @param calciumTroponin Current guess for Ca_trop value
     *  @param z current guess for z
     *  @param Q = Q1+Q2+Q3, where Qi already computed at next timestep
     *  @param residualComponent1 Returned value - first component of residual
     *  @param residualComponent2 Returned value - second component of residual
     */
    void CalculateBackwardEulerResidual(double calciumTroponin, double z, double Q,
                                        double& residualComponent1, double& residualComponent2);

public :
    /**
     *  Constructor
     */
    NhsModelWithBackwardSolver();


    /**
     *  Solves for the new state variables at the given end time using the implicit
     *  method. Note that the internal state variables are not altered, the solution
     *  is saved instead. Call UpdateStateVariables() to update, and
     *  GetNextActiveTension() to get the solved active tension
     *
     *  The state variables are not updated because this solve will be called as part
     *  of the newton iteration (ie guess stretch, see what the new active tension is)
     *  in a fully implicit method.
     *
     *  Note: overloaded from the method in AbstractOdeBasedContractionModel, which
     *  just does a simple Euler solve
     *
     *  @param startTime
     *  @param endTime
     *  @param timestep
     */
    void RunDoNotUpdate(double startTime, double endTime, double timestep);

    /**
     *  @return the active tension corresponding to the stored state variables computed
     *  from the last RunDoNotUpdate(), ie the active tension at the next time.
     *  Note that calling GetActiveTension() on the base class will use the internal
     *  state variables and return the active tension at the last time, if
     *  RunDoNotUpdate() has been called but UpdateStateVariables() has not
     */
    double GetNextActiveTension();

    /**
     *  Overload the RunAndUpdate() method too, as that would use the base class's default (euler) solver
     *
     *  @param startTime
     *  @param endTime
     *  @param timestep
     */
    void RunAndUpdate(double startTime, double endTime, double timestep);
};

#endif /*NHSMODELWITHBACKWARDSOLVER_HPP_*/
