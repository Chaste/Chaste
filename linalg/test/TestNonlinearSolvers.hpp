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

PetscErrorCode ComputeTestResidual(SNES snes,Vec solution_guess,Vec residual,void* pContext);
PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat* pJacobian ,Mat* pPreconditioner,MatStructure* pMatStructure ,void* pContext);
PetscErrorCode ComputeTestResidual3d(SNES snes,Vec solution_guess,Vec residual,void* pContext);
PetscErrorCode ComputeTestJacobian3d(SNES snes,Vec input,Mat* pJacobian ,Mat* pPreconditioner,MatStructure* pMatStructure ,void* pContext);

class TestNonlinearSolvers : public CxxTest::TestSuite
{
public:
    void TestNonlinearProblemException() throw (Exception)
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
        VecDestroy(initial_guess);
    }

    void TestOn2dNonlinearProblem() throw (Exception)
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
            // the solution is x = 1/sqrt(2), y = 1/sqrt(2)
            TS_ASSERT_DELTA(answer_petsc_repl[i] ,1/sqrt(2.0),tol);
            TS_ASSERT_DELTA(answer_newton_repl[i],1/sqrt(2.0),tol);
        }

        VecDestroy(initial_guess);
        VecDestroy(answer_petsc);
        VecDestroy(answer_newton);
    }

    void TestOn3dNonlinearProblem() throw (Exception)
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
            // the solution is x = 1/sqrt(3), y = 1/sqrt(3),  z = 1/sqrt(3)
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

        VecDestroy(residual);
        VecDestroy(initial_guess);
        VecDestroy(answer_petsc);
        VecDestroy(answer_newton);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////
// global functions called by nonlinear solvers
///////////////////////////////////////////////////////////////////////////////////////////////

PetscErrorCode ComputeTestResidual(SNES snes,Vec solution_guess,Vec residual,void* pContext)
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

PetscErrorCode ComputeTestJacobian(SNES snes,Vec input,Mat* pJacobian ,Mat* pPreconditioner,MatStructure* pMatStructure ,void* pContext)
{
    double x, y;

    ReplicatableVector input_replicated;
    input_replicated.ReplicatePetscVector(input);
    x = input_replicated[0];
    y = input_replicated[1];

    PetscMatTools::SetElement(*pJacobian, 0 , 0 , 2.0*x );
    PetscMatTools::SetElement(*pJacobian, 0 , 1 , 2.0*y);
    PetscMatTools::SetElement(*pJacobian, 1 , 0 , 1.0);
    PetscMatTools::SetElement(*pJacobian, 1 , 1 , -1.0);
    MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);

    return 0;
}

PetscErrorCode ComputeTestResidual3d(SNES snes,Vec solution_guess,Vec residual,void* pContext)
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

PetscErrorCode ComputeTestJacobian3d(SNES snes,Vec input,Mat* pJacobian ,Mat* pPreconditioner,MatStructure* pMatStructure ,void* pContext)
{
    double x, y, z;

    ReplicatableVector input_replicated;
    input_replicated.ReplicatePetscVector(input);

    x = input_replicated[0];
    y = input_replicated[1];
    z = input_replicated[2];

    PetscMatTools::SetElement(*pJacobian, 0 , 0 , 2.0*x );
    PetscMatTools::SetElement(*pJacobian, 0 , 1 , 2.0*y);
    PetscMatTools::SetElement(*pJacobian, 0 , 2 , 2.0*z);
    PetscMatTools::SetElement(*pJacobian, 1 , 0 , 1.0);
    PetscMatTools::SetElement(*pJacobian, 1 , 1 , -1.0);
    PetscMatTools::SetElement(*pJacobian, 1 , 2 , 0.0);
    PetscMatTools::SetElement(*pJacobian, 2 , 0 , 0.0);
    PetscMatTools::SetElement(*pJacobian, 2 , 1 , 1.0);
    PetscMatTools::SetElement(*pJacobian, 2 , 2 , -1.0);
    MatAssemblyBegin(*pJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*pJacobian,MAT_FINAL_ASSEMBLY);

    return 0;
}

#endif //_TESTNONLINEARSOLVERS_HPP_
