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

/**
 * Concrete Simple Nonlinear PDE system solver.
 */

#include "SimplePetscNonlinearSolver.hpp"
#include "Exception.hpp"
#include "petscsnes.h"
#include "PetscTools.hpp"
#include <sstream>

Vec SimplePetscNonlinearSolver::Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                                      Vec initialGuess,
                                      unsigned fill,
                                      void* pContext)
{
    SNES snes;

    // Create the residual vector by copying the structure of the initial guess
    Vec residual;
    VecDuplicate(initialGuess, &residual);

    Mat jacobian; // Jacobian Matrix

    PetscInt N; // number of elements

    // Get the size of the Jacobian from the residual
    VecGetSize(initialGuess,&N);

    PetscTools::SetupMat(jacobian, N, N, fill);

    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFunction(snes, residual, pComputeResidual, pContext);
    SNESSetJacobian(snes, jacobian, jacobian, pComputeJacobian, pContext);
    SNESSetType(snes,SNESLS);
    SNESSetTolerances(snes,1.0e-5,1.0e-5,1.0e-5,PETSC_DEFAULT,PETSC_DEFAULT);

    // x is the iteration vector SNES uses when solving, set equal to initialGuess to start with
    Vec x;
    VecDuplicate(initialGuess, &x);
    VecCopy(initialGuess, x);

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    SNESSolve(snes, x);
#else
    SNESSolve(snes, PETSC_NULL, x);
#endif

    VecDestroy(residual);
    MatDestroy(jacobian); // Free Jacobian

    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);
#define COVERAGE_IGNORE
    if (reason < 0)
    {
        std::stringstream reason_stream;
        reason_stream << reason;
        VecDestroy(x); // Since caller can't free the memory in this case
        SNESDestroy(snes);
        EXCEPTION("Nonlinear Solver did not converge. PETSc reason code:"
                  +reason_stream.str()+" .");
    }
#undef COVERAGE_IGNORE
    SNESDestroy(snes);

    return x;
}
