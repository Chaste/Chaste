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

#ifndef _TESTNONLINEARSOLVERS_HPP_
#define _TESTNONLINEARSOLVERS_HPP_

#include <cxxtest/TestSuite.h>
#include "SimpleNewtonNonlinearSolver.hpp"
#include "SimplePetscNonlinearSolver.hpp"
#include <iostream>
#include <cmath>
#include "ReplicatableVector.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "limits.h"

PetscErrorCode ComputeTestResidual(SNES snes, Vec solution_guess, Vec residual, void* pContext);
PetscErrorCode ComputeTestResidual3d(SNES snes, Vec solution_guess, Vec residual, void* pContext);
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
PetscErrorCode ComputeTestJacobian(SNES snes, Vec input, Mat jacobian, Mat preconditioner, void* pContext);
PetscErrorCode ComputeTestJacobian3d(SNES snes, Vec input, Mat jacobian, Mat preconditioner, void* pContext);
#else
PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat* pJacobian ,Mat* pPreconditioner,MatStructure* pMatStructure ,void* pContext);
PetscErrorCode ComputeTestJacobian3d(SNES snes,Vec input,Mat* pJacobian ,Mat* pPreconditioner,MatStructure* pMatStructure ,void* pContext);
#endif

class TestNonlinearSolvers : public CxxTest::TestSuite
{
public:
    void TestNonlinearProblemException()
    {
        SimpleNewtonNonlinearSolver solver_newton;

        // Set up solution guess for residuals
        int length = 2;

        // Set up initial Guess
        Vec initial_guess = PetscTools::CreateVec(length);
        PetscVecTools::SetElement(initial_guess, 0, -1e16);
        PetscVecTools::SetElement(initial_guess, 1, -1e8);
        PetscVecTools::Finalise(initial_guess);
        TS_ASSERT_THROWS_THIS(solver_newton.Solve(&ComputeTestResidual, &(ComputeTestJacobian), initial_guess, length, NULL),
                "Iteration 27, unable to find damping factor such that residual decreases in update direction");
        PetscTools::Destroy(initial_guess);
    }

    void TestOn2dNonlinearProblem()
    {
        SimplePetscNonlinearSolver solver_petsc;
        SimpleNewtonNonlinearSolver solver_newton;

        // Set up solution guess for residuals
        int length = 2;

        // Set up initial Guess
        Vec initial_guess=PetscTools::CreateVec(length);
        PetscVecTools::SetElement(initial_guess, 0, 1.0);
        PetscVecTools::SetElement(initial_guess, 1, 1.0);
        PetscVecTools::Finalise(initial_guess);

        // Solve using petsc solver
        Vec answer_petsc = solver_petsc.Solve(&ComputeTestResidual, &ComputeTestJacobian,
                                              initial_guess, length, NULL);

        // Solve using newton method
        Vec answer_newton = solver_newton.Solve(&ComputeTestResidual, &ComputeTestJacobian,
                                                initial_guess, length, NULL);

        // Replicate the answers so we can access them without worrying about parallel stuff
        ReplicatableVector answer_petsc_repl(answer_petsc);
        ReplicatableVector answer_newton_repl(answer_newton);

        double tol = 1e-4;

        for (int i=0; i<2; i++)
        {
            // the solution is x = 1/sqrt(2.0), y = 1/sqrt(2.0)
            TS_ASSERT_DELTA(answer_petsc_repl[i] ,1/sqrt(2.0),tol);
            TS_ASSERT_DELTA(answer_newton_repl[i],1/sqrt(2.0),tol);
        }

        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer_petsc);
        PetscTools::Destroy(answer_newton);
    }

    void TestOn3dNonlinearProblem()
    {
        SimplePetscNonlinearSolver solver_petsc;
        SimpleNewtonNonlinearSolver solver_newton;

        // Set up solution guess for residuals
        int length = 3;

        // Set up initial Guess
        Vec initial_guess=PetscTools::CreateVec(length);
        PetscVecTools::SetElement(initial_guess, 0, 1);
        PetscVecTools::SetElement(initial_guess, 1, 1);
        PetscVecTools::SetElement(initial_guess, 2, 1);
        PetscVecTools::Finalise(initial_guess);

        // Solve using petsc solver
        Vec answer_petsc = solver_petsc.Solve(&ComputeTestResidual3d, &ComputeTestJacobian3d,
                                              initial_guess, length, NULL);

        // Solve using newton method
        solver_newton.SetTolerance(1e-10);                      // to cover this method
        solver_newton.SetWriteStats();                          // to cover this method
        Vec answer_newton = solver_newton.Solve(&ComputeTestResidual3d, &ComputeTestJacobian3d,
                                                initial_guess, length, NULL);

        // Replicate the answers so we can access them without worrying about parallel stuff
        ReplicatableVector answer_petsc_repl(answer_petsc);
        ReplicatableVector answer_newton_repl(answer_newton);

        double tol = 1e-6;

        for (int i=0; i<3; i++)
        {
            // the solution is x = 1/sqrt(3.0), y = 1/sqrt(3.0),  z = 1/sqrt(3.0)
            TS_ASSERT_DELTA(answer_petsc_repl[i] ,1/sqrt(3.0),tol);
            TS_ASSERT_DELTA(answer_newton_repl[i],1/sqrt(3.0),tol);
        }

        // Check the residual really did have scaled norm within the tolerance
        Vec residual;
        VecDuplicate(answer_newton, &residual);
        ComputeTestResidual3d(NULL, answer_newton, residual, NULL);
        double norm;
        VecNorm(residual, NORM_2, &norm);
        TS_ASSERT_LESS_THAN(norm/length, 1e-10);

        PetscTools::Destroy(residual);
        PetscTools::Destroy(initial_guess);
        PetscTools::Destroy(answer_petsc);
        PetscTools::Destroy(answer_newton);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////
// global functions called by nonlinear solvers
///////////////////////////////////////////////////////////////////////////////////////////////

PetscErrorCode ComputeTestResidual(SNES snes, Vec solution_guess, Vec residual, void* pContext)
{
    double x, y;

    ReplicatableVector solution_guess_replicated;
    solution_guess_replicated.ReplicatePetscVector(solution_guess);
    x = solution_guess_replicated[0];
    y = solution_guess_replicated[1];

    PetscVecTools::SetElement(residual,0,x*x+y*y-1);
    PetscVecTools::SetElement(residual,1,x-y);
    PetscVecTools::Finalise(residual);
    return 0;
}

PetscErrorCode ComputeTestResidual3d(SNES snes, Vec solution_guess, Vec residual, void* pContext)
{
    double x, y, z;

    ReplicatableVector solution_guess_replicated;
    solution_guess_replicated.ReplicatePetscVector(solution_guess);

    x = solution_guess_replicated[0];
    y = solution_guess_replicated[1];
    z = solution_guess_replicated[2];

    PetscVecTools::SetElement(residual,0,x*x+y*y+z*z-1);
    PetscVecTools::SetElement(residual,1,x-y);
    PetscVecTools::SetElement(residual,2,y-z);
    PetscVecTools::Finalise(residual);

    return 0;
}

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
PetscErrorCode ComputeTestJacobian(SNES snes, Vec input, Mat jacobian, Mat preconditioner, void* pContext)
{
#else
PetscErrorCode ComputeTestJacobian(SNES snes, Vec input, Mat* pJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext)
{
    Mat jacobian = *pJacobian;
#endif
    double x, y;

    ReplicatableVector input_replicated;
    input_replicated.ReplicatePetscVector(input);
    x = input_replicated[0];
    y = input_replicated[1];

    PetscMatTools::SetElement(jacobian, 0 , 0 , 2.0*x );
    PetscMatTools::SetElement(jacobian, 0 , 1 , 2.0*y);
    PetscMatTools::SetElement(jacobian, 1 , 0 , 1.0);
    PetscMatTools::SetElement(jacobian, 1 , 1 , -1.0);
    PetscMatTools::Finalise(jacobian);

    return 0;
}


#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
PetscErrorCode ComputeTestJacobian3d(SNES snes, Vec input, Mat jacobian, Mat preconditioner, void* pContext)
{
#else
PetscErrorCode ComputeTestJacobian3d(SNES snes, Vec input, Mat* pJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext)
{
    Mat jacobian = *pJacobian;
#endif
    double x, y, z;

    ReplicatableVector input_replicated;
    input_replicated.ReplicatePetscVector(input);

    x = input_replicated[0];
    y = input_replicated[1];
    z = input_replicated[2];

    PetscMatTools::SetElement(jacobian, 0 , 0 , 2.0*x );
    PetscMatTools::SetElement(jacobian, 0 , 1 , 2.0*y);
    PetscMatTools::SetElement(jacobian, 0 , 2 , 2.0*z);
    PetscMatTools::SetElement(jacobian, 1 , 0 , 1.0);
    PetscMatTools::SetElement(jacobian, 1 , 1 , -1.0);
    PetscMatTools::SetElement(jacobian, 1 , 2 , 0.0);
    PetscMatTools::SetElement(jacobian, 2 , 0 , 0.0);
    PetscMatTools::SetElement(jacobian, 2 , 1 , 1.0);
    PetscMatTools::SetElement(jacobian, 2 , 2 , -1.0);
    PetscMatTools::Finalise(jacobian);

    return 0;
}

#endif //_TESTNONLINEARSOLVERS_HPP_
