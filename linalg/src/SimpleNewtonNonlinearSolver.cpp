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

#include "SimpleNewtonNonlinearSolver.hpp"
#include <iostream>
#include <cassert>
#include "Exception.hpp"

SimpleNewtonNonlinearSolver::SimpleNewtonNonlinearSolver()
        : mTolerance(1e-5),
        mWriteStats(false)
{
    mTestDampingValues.push_back(-0.1);
    mTestDampingValues.push_back(0.05);
    for (unsigned i=1; i<=12; i++)
    {
        double val = double(i)/10;
        mTestDampingValues.push_back(val);
    }
}

SimpleNewtonNonlinearSolver::~SimpleNewtonNonlinearSolver()
{}

Vec SimpleNewtonNonlinearSolver::Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
                                       PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat,Mat,void*),
#else
                                       PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
#endif
                                       Vec initialGuess,
                                       unsigned fill,
                                       void* pContext)
{
    PetscInt size;
    VecGetSize(initialGuess, &size);

    Vec current_solution;
    VecDuplicate(initialGuess, &current_solution);
    VecCopy(initialGuess, current_solution);

    // The "false" says that we are allowed to do new mallocs without PETSc 3.3 causing an error
    LinearSystem linear_system(current_solution, fill, false);

    (*pComputeResidual)(nullptr, current_solution, linear_system.rGetRhsVector(), pContext);


    double residual_norm;
    VecNorm(linear_system.rGetRhsVector(), NORM_2, &residual_norm);
    double scaled_residual_norm = residual_norm/size;

    if (mWriteStats)
    {
        std::cout << "Newton's method:\n  Initial ||residual||/N = " << scaled_residual_norm
                  << "\n  Attempting to solve to tolerance " << mTolerance << "..\n";
    }

    double old_scaled_residual_norm;
    unsigned counter = 0;
    while (scaled_residual_norm > mTolerance)
    {
        counter++;

        // Store the old norm to check with the new later
        old_scaled_residual_norm = scaled_residual_norm;

        // Compute Jacobian and solve J dx = f for the (negative) update dx, (J the jacobian, f the residual)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
        (*pComputeJacobian)(nullptr, current_solution, (linear_system.rGetLhsMatrix()), nullptr, pContext);
#else
        (*pComputeJacobian)(nullptr, current_solution, &(linear_system.rGetLhsMatrix()), nullptr, nullptr, pContext);
#endif

        Vec negative_update = linear_system.Solve();


        Vec test_vec;
        VecDuplicate(initialGuess, &test_vec);

        double best_damping_factor = 1.0;
        double best_scaled_residual = 1e20; // large

        // Loop over all the possible damping value and determine which gives smallest residual
        for (unsigned i=0; i<mTestDampingValues.size(); i++)
        {
            // Note: WAXPY calls VecWAXPY(w,a,x,y) which computes w = ax+y
            PetscVecTools::WAXPY(test_vec,-mTestDampingValues[i],negative_update,current_solution);

            // Compute new residual
            linear_system.ZeroLinearSystem();
            (*pComputeResidual)(nullptr, test_vec, linear_system.rGetRhsVector(), pContext);
            VecNorm(linear_system.rGetRhsVector(), NORM_2, &residual_norm);
            scaled_residual_norm = residual_norm/size;

            if (scaled_residual_norm < best_scaled_residual)
            {
                best_scaled_residual = scaled_residual_norm;
                best_damping_factor = mTestDampingValues[i];
            }
        }
        PetscTools::Destroy(test_vec);

        // Check the smallest residual was actually smaller than the previous; if not, quit
        if (best_scaled_residual > old_scaled_residual_norm)
        {
            // Free memory
            PetscTools::Destroy(current_solution);
            PetscTools::Destroy(negative_update);

            // Raise error
            EXCEPTION("Iteration " << counter << ", unable to find damping factor such that residual decreases in update direction");
        }

        if (mWriteStats)
        {
            std::cout << "    Best damping factor = " << best_damping_factor << "\n";
        }

        // Update solution: current_guess = current_solution - best_damping_factor*negative_update
        PetscVecTools::AddScaledVector(current_solution, negative_update, -best_damping_factor);
        scaled_residual_norm = best_scaled_residual;
        PetscTools::Destroy(negative_update);

        // Compute best residual vector again and store in linear_system for next Solve()
        linear_system.ZeroLinearSystem();
        (*pComputeResidual)(nullptr, current_solution, linear_system.rGetRhsVector(), pContext);

        if (mWriteStats)
        {
            std::cout << "    Iteration " << counter << ": ||residual||/N = " << scaled_residual_norm << "\n";
        }
    }

    if (mWriteStats)
    {
        std::cout << "  ..solved!\n\n";
    }

    return current_solution;
}

void SimpleNewtonNonlinearSolver::SetTolerance(double tolerance)
{
    assert(tolerance > 0);
    mTolerance = tolerance;
}
