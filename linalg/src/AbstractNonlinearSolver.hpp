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

#ifndef _ABSTRACTNONLINEARSOLVER_HPP_
#define _ABSTRACTNONLINEARSOLVER_HPP_

/**
 * Abstract Nonlinear equation system solver, dictates that each solver must have a
 * Solve function.
 */

#include <petscsnes.h>

class AbstractNonlinearSolver
{
public:

    /**
     * Pure virtual Solve() method.
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
     */
    virtual Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat,Mat,void*),
#else
                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
#endif
                      Vec initialGuess,
                      unsigned fill,
                      void* pContext)=0;

    /**
     * Destructor.
     */
    virtual ~AbstractNonlinearSolver()
    {}
};

#endif // _ABSTRACTNONLINEARSOLVER_HPP_
