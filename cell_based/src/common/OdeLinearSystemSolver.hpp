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

#ifndef ODELINEARSYSTEMSOLVER_HPP_
#define ODELINEARSYSTEMSOLVER_HPP_

#include "LinearSystem.hpp"

/**
 * Solve large systems of ODEs of the form
 *
 * M dr/dt = f(t,r),
 *
 * where M is a matrix, typically large and sparse, and r and f are vectors.
 *
 * This differs from the OdeSolver classes in the ode component, which are for
 * systems of the form dr/dt = f(t,r) and for small numbers of unknowns.
 *
 * The solver uses forward euler to dicretise as M r^{n+1} = M r^{n} + dt f
 * and solves this linear system.
 *
 * The calling code is responsible with setting up M and f each timestep.
 */
class OdeLinearSystemSolver
{
private:

    /** Timestep for solver. */
    double mTimeStep;

    /** The LHS matrix and the force vector. Solve() is called on this. */
    LinearSystem mLinearSystem;

    /** Force vector (f in M dr/dt = f). */
    Vec mForceVector;

    /** Solution at current timestep. */
    Vec mCurrentSolution;

public:

    /**
     * Constructor.
     *
     * @param systemSize size of the ODE system
     * @param timeStep the time step used to integrate the ODE system
     */
    OdeLinearSystemSolver(unsigned systemSize, double timeStep);

    /** Destructor. */
    ~OdeLinearSystemSolver();

    /** @return the timestep for the solver. */
    double GetTimeStep();

    /** @return the matrix. */
    Mat& rGetLhsMatrix();

    /** @return the force vector. */
    Vec& rGetForceVector();

    /**
     * Set the initial conditions.
     *
     * @param initialConditionsVector the initial condition
     */
    void SetInitialConditionVector(Vec initialConditionsVector);

    /** Solve the ODE system over one time step.
     * @return solution vector
     */
    Vec SolveOneTimeStep();
};

#endif /*ODELINEARSYSTEMSOLVER_HPP_*/
