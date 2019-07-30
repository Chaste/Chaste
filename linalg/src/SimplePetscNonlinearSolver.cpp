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
 * Concrete Simple Nonlinear PDE system solver.
 */

#include <sstream>
#include "petscsnes.h"
#include "SimplePetscNonlinearSolver.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

Vec SimplePetscNonlinearSolver::Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
                                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat,Mat,void*),
#else
                                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
#endif
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
    VecGetSize(initialGuess, &N);

    // Note that the Jacobian matrix might involve new non-zero elements in the course of a SNES solve
    PetscTools::SetupMat(jacobian, N, N, fill, PETSC_DECIDE, PETSC_DECIDE, true, false /*malloc flag*/);

    SNESCreate(PETSC_COMM_WORLD, &snes);

    SNESSetFunction(snes, residual, pComputeResidual, pContext);
    SNESSetJacobian(snes, jacobian, jacobian, pComputeJacobian, pContext);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 4) //PETSc 3.4 or later
    SNESSetType(snes, SNESNEWTONLS);
#else
    SNESSetType(snes, SNESLS);
#endif

    SNESSetTolerances(snes,1.0e-5,1.0e-5,1.0e-5,PETSC_DEFAULT,PETSC_DEFAULT);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 3) //PETSc 3.3
    SNESLineSearch linesearch;
    SNESGetSNESLineSearch(snes, &linesearch);
    SNESLineSearchSetType(linesearch, "bt"); //Use backtracking search as default
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 4) //PETSc 3.4 or later
    SNESLineSearch linesearch;
    SNESGetLineSearch(snes, &linesearch);
    SNESLineSearchSetType(linesearch, "bt"); //Use backtracking search as default
#endif

    // x is the iteration vector SNES uses when solving, set equal to initialGuess to start with
    Vec x;

    VecDuplicate(initialGuess, &x);
    VecCopy(initialGuess, x);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 5)
    // Seems to want the preconditioner to be explicitly set to none now
    // Copied this from the similar PETSc example at:
    // http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex1.c
    // Which got it to work...
    KSP ksp;
    SNESGetKSP(snes,&ksp);
    PC pc;
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCNONE);
#endif

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    SNESSolve(snes, x);
#else
    SNESSolve(snes, PETSC_NULL, x);
#endif
    PetscTools::Destroy(residual);
    PetscTools::Destroy(jacobian); // Free Jacobian

    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);
// LCOV_EXCL_START
    if (reason < 0)
    {
        std::stringstream reason_stream;
        reason_stream << reason;
        PetscTools::Destroy(x); // Since caller can't free the memory in this case
        SNESDestroy(PETSC_DESTROY_PARAM(snes));
        EXCEPTION("Nonlinear Solver did not converge. PETSc reason code:"
                  +reason_stream.str()+" .");
    }
// LCOV_EXCL_STOP
    SNESDestroy(PETSC_DESTROY_PARAM(snes));

    return x;
}
