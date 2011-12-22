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

#ifndef TESTODELINEARSYSTEMSOLVER_HPP_
#define TESTODELINEARSYSTEMSOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"
#include "ReplicatableVector.hpp"
#include "OdeLinearSystemSolver.hpp"

// The header file below must be included in any test that uses Petsc
#include "PetscSetupAndFinalize.hpp"

class TestOdeLinearSystemSolver : public CxxTest::TestSuite
{
public:

    // Solve a trivial ODE linear system M dr/dt = f
    // where M = [1 0 ; 0 2] and f = [1 3];
    void TestWithTrivial2dProblem()
    {
        // Declare solver and give the size of the system and timestep
        unsigned system_size = 2;
        double dt = 0.01;
        OdeLinearSystemSolver solver(system_size, dt);

        // Set up the matrix
        Mat& r_matrix = solver.rGetLhsMatrix();
        PetscMatTools::SetElement(r_matrix, 0, 0, 1.0);
        PetscMatTools::SetElement(r_matrix, 1, 0, 0.0);
        PetscMatTools::SetElement(r_matrix, 0, 1, 0.0);
        PetscMatTools::SetElement(r_matrix, 1, 1, 2.0);
        PetscMatTools::Finalise(r_matrix);

        // Initial condition
        Vec initial_condition = PetscTools::CreateAndSetVec(2, 0.0);
        PetscVecTools::SetElement(initial_condition, 0, 10.0);
        PetscVecTools::SetElement(initial_condition, 1, 11.0);

        solver.SetInitialConditionVector(initial_condition);

        // Then an rGetVector for RHS
        Vec& r_force_vector = solver.rGetForceVector();
        PetscVecTools::SetElement(r_force_vector, 0, 1.0);
        PetscVecTools::SetElement(r_force_vector, 1, 3.0);

        // Solve to get solution at next timestep
        Vec soln_next_timestep = solver.SolveOneTimeStep();

        ReplicatableVector soln_next_timestep_repl(soln_next_timestep);

        TS_ASSERT_DELTA(soln_next_timestep_repl[0], 10.0 + dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl[1], 11.0 + 1.5*dt, 1e-6);

        // Solve again, with the same force
        soln_next_timestep = solver.SolveOneTimeStep();

        ReplicatableVector soln_next_timestep_repl2(soln_next_timestep);

        TS_ASSERT_DELTA(soln_next_timestep_repl2[0], 10.0 + 2*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl2[1], 11.0 + 3*dt, 1e-6);

        // Tidy up
        VecDestroy(initial_condition);
    }

    // Solve a simple ODE linear system M dr/dt = f
    // where M = [0 1; 1 0] and f = [1 2];
    void TestWithSimple2dProblem()
    {
        // Declare solver and give the size of the system and timestep
        unsigned system_size = 2;
        double dt = 0.01;
        OdeLinearSystemSolver solver(system_size, dt);

        // Set up the matrix
        Mat& r_matrix = solver.rGetLhsMatrix();
        PetscMatTools::SetElement(r_matrix, 0, 0, 0.0);
        PetscMatTools::SetElement(r_matrix, 1, 0, 1.0);
        PetscMatTools::SetElement(r_matrix, 0, 1, 1.0);
        PetscMatTools::SetElement(r_matrix, 1, 1, 0.0);
        PetscMatTools::Finalise(r_matrix);

        // Initial condition
        Vec initial_condition = PetscTools::CreateAndSetVec(2, 0.0);
        PetscVecTools::SetElement(initial_condition, 0, 10.0);
        PetscVecTools::SetElement(initial_condition, 1, 11.0);

        solver.SetInitialConditionVector(initial_condition);

        // Then an rGetForceVector for RHS
        Vec& r_vector = solver.rGetForceVector();
        PetscVecTools::SetElement(r_vector, 0, 1.0);
        PetscVecTools::SetElement(r_vector, 1, 2.0);

        // Solve to get solution at next timestep
        Vec soln_next_timestep = solver.SolveOneTimeStep();

        ReplicatableVector soln_next_timestep_repl(soln_next_timestep);

        TS_ASSERT_DELTA(soln_next_timestep_repl[0], 10.0 + 2*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl[1], 11.0 +   dt, 1e-6);

        // Solve again, with the same force
        soln_next_timestep = solver.SolveOneTimeStep();

        ReplicatableVector soln_next_timestep_repl2(soln_next_timestep);

        TS_ASSERT_DELTA(soln_next_timestep_repl2[0], 10.0 + 4*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl2[1], 11.0 + 2*dt, 1e-6);

        // Tidy up
        VecDestroy(initial_condition);
    }

    // Solve a simple ODE linear system M dr/dt = f
    // where M = [1 2 3 4; 4 1 2 3; 3 2 2 1; 3 2 1 1] and f = [1 2 3 4];
    void TestWithSimple4dProblem()
    {
        // Declare solver and give the size of the system and timestep
        unsigned system_size = 4;
        double dt = 0.01;
        OdeLinearSystemSolver solver(system_size, dt);

        // Set up the matrix
        Mat& r_matrix = solver.rGetLhsMatrix();
        PetscMatTools::SetElement(r_matrix, 0, 0, 1.0);
        PetscMatTools::SetElement(r_matrix, 0, 1, 2.0);
        PetscMatTools::SetElement(r_matrix, 0, 2, 3.0);
        PetscMatTools::SetElement(r_matrix, 0, 3, 4.0);
        PetscMatTools::SetElement(r_matrix, 1, 0, 4.0);
        PetscMatTools::SetElement(r_matrix, 1, 1, 1.0);
        PetscMatTools::SetElement(r_matrix, 1, 2, 2.0);
        PetscMatTools::SetElement(r_matrix, 1, 3, 3.0);
        PetscMatTools::SetElement(r_matrix, 2, 0, 3.0);
        PetscMatTools::SetElement(r_matrix, 2, 1, 2.0);
        PetscMatTools::SetElement(r_matrix, 2, 2, 2.0);
        PetscMatTools::SetElement(r_matrix, 2, 3, 1.0);
        PetscMatTools::SetElement(r_matrix, 3, 0, 3.0);
        PetscMatTools::SetElement(r_matrix, 3, 1, 2.0);
        PetscMatTools::SetElement(r_matrix, 3, 2, 1.0);
        PetscMatTools::SetElement(r_matrix, 3, 3, 1.0);
        PetscMatTools::Finalise(r_matrix);

        // Initial condition
        Vec initial_condition = PetscTools::CreateAndSetVec(4, 0.0);
        PetscVecTools::SetElement(initial_condition, 0, 10.0);
        PetscVecTools::SetElement(initial_condition, 1, 11.0);
        PetscVecTools::SetElement(initial_condition, 2, 12.0);
        PetscVecTools::SetElement(initial_condition, 3, 13.0);

        solver.SetInitialConditionVector(initial_condition);

        // Then an rGetForceVector for RHS
        Vec& r_vector = solver.rGetForceVector();
        PetscVecTools::SetElement(r_vector, 0, 1.0);
        PetscVecTools::SetElement(r_vector, 1, 2.0);
        PetscVecTools::SetElement(r_vector, 2, 3.0);
        PetscVecTools::SetElement(r_vector, 3, 4.0);

        // Solve to get solution at next timestep
        Vec soln_next_timestep = solver.SolveOneTimeStep();

        ReplicatableVector soln_next_timestep_repl(soln_next_timestep);

        TS_ASSERT_DELTA(soln_next_timestep_repl[0], 10.0 + 0.56*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl[1], 11.0 + 1.64*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl[2], 12.0 - 1.00*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl[3], 13.0 + 0.04*dt, 1e-6);

        // Solve again, with the same force
        soln_next_timestep = solver.SolveOneTimeStep();

        ReplicatableVector soln_next_timestep_repl2(soln_next_timestep);

        TS_ASSERT_DELTA(soln_next_timestep_repl2[0], 10.0 + 2*0.56*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl2[1], 11.0 + 2*1.64*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl2[2], 12.0 - 2*1.00*dt, 1e-6);
        TS_ASSERT_DELTA(soln_next_timestep_repl2[3], 13.0 + 2*0.04*dt, 1e-6);

        // Tidy up
        VecDestroy(initial_condition);
    }
};

#endif /*TESTODELINEARSYSTEMSOLVER_HPP_*/
