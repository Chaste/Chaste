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
     *  Get the active tension corresponding to the stored state variables computed
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
