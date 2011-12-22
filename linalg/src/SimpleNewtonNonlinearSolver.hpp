/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef SIMPLENEWTONNONLINEARSOLVER_HPP_
#define SIMPLENEWTONNONLINEARSOLVER_HPP_

#include "LinearSystem.hpp"
#include "AbstractNonlinearSolver.hpp"
#include <vector>

/**
 * A solver class that uses Newton's method with damping.
 */
class SimpleNewtonNonlinearSolver : public AbstractNonlinearSolver
{
private:

    double mLinearSolverRelativeTolerance;   /**< The linear solver relative tolerance. */
    double mTolerance;                       /**< The tolerance (set to 1e-5 in constructor). */
    bool mWriteStats;                        /**< Whether the solver writes details as it solves (set to false in constructor). */

    std::vector<double> mTestDampingValues;  /**< Vector of possible damping factors (set in the constructor). */

public:

    /**
     * Constructor.
     *
     * @param linearSolverRelativeTolerance defaults to 1e-6
     */
    SimpleNewtonNonlinearSolver(double linearSolverRelativeTolerance = 1e-6);

    /**
     * Destructor.
     */
    virtual ~SimpleNewtonNonlinearSolver();

    /**
     * Solve().
     *
     * Solve a nonlinear system using Newton's method with damping. Newton's algorithm
     * is
     *
     * x_new = x_old - J^{-1} f
     *
     * where J is the Jacobian matrix evaluated at x_old and f the residual evaluated at
     * x_old. The Newton method with damping is
     *
     * x_new = x_old - s J^{-1} f
     *
     * where s is some damping factor. Here s is chosen by just looked at a fixed set of
     * possible damping factors and choosing the one which gives the best x_new (the one
     * for which the residual evaluated at x_new has the lowest norm).
     *
     * The solver quits once the ||f||/numVariables
     *
     * @param pComputeResidual points to the function which
     * computes the residual, it must take arguments SNES (a PETSc nonlinear solver
     * object), Vec (current guess - a vector of the correct size), Vec (a Vec of the
     * correct size in which the residual is returned), void* (a pointer to
     * anything you may need to refer to when calculating the residual)
     *
     * @param pComputeJacobian points to
     * the function which computes the Jacobian, it must take arguments SNES (a PETSc
     * nonlinear solver * object), Mat* (a pointer to the Jacobian matrix) ,Mat* (a pointer
     * to a preconditioner matrix), MatStructure* (points to the PETSc matrix type e.g. AIJ), void* (a pointer to
     * anything you may need to refer to when calculating the residual).
     *
     * @param initialGuess A PETSc Vec of the correct size, containing initial guesses
     * for the nonlinear solver.
     *
     * @param pContext [optional] A pointer to a class that may have to be used in the
     *  ComputeResidual and ComputeJacobian functions
     *
     * @param fill the expected maximum number of nonzeros in a row of the Jacobian matrix
     *
     * @return Returns a PETSc Vec of the solution.
     *
     * To be used in the form:
     * Vec answer=solver->Solve(&ComputeResidual, &ComputeJacobian, initialGuess, NULL);
     *
     * In the same file, but outside this class the functions ComputeResidual and
     * ComputeJacobian must sit, using the input arguments specified above.
     */
    virtual Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                      Vec initialGuess,
                      unsigned fill,
                      void* pContext);

    /**
     * Set a tolerance other than the default.
     *
     * @param tolerance
     */
    void SetTolerance(double tolerance);

    /**
     * Call to set the solver to write details as it solves.
     *
     * @param writeStats defaults to true
     */
    void SetWriteStats(bool writeStats = true)
    {
        mWriteStats = writeStats;
    }
};

#endif /*SIMPLENEWTONNONLINEARSOLVER_HPP_*/
