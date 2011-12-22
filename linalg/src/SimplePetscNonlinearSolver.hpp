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

#ifndef _SIMPLENONLINEARSOLVER_HPP_
#define _SIMPLENONLINEARSOLVER_HPP_

/**
 * Concrete Simple Nonlinear PDE system solver.
 */

#include <petscsnes.h>
#include "AbstractNonlinearSolver.hpp"

class SimplePetscNonlinearSolver : public AbstractNonlinearSolver
{
public:

/**
 * Simple Nonlinear PDE system solver, uses the PETSc SNES Solver, which uses Newton's method.
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
 *  for the nonlinear solver.
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
 *
 */
    Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
              PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
              Vec initialGuess,
              unsigned fill,
              void* pContext);
};

#endif // _SIMPLENONLINEARSOLVER_HPP_
