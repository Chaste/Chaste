/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTPCBLOCKDIAGONAL_HPP_
#define TESTPCBLOCKDIAGONAL_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "Timer.hpp"
#include "DistributedVectorFactory.hpp"
#include <cstring>
#include "Warnings.hpp"

/*
 * Warning: these tests do not inform PETSc about the nullspace of the matrix. Therefore, convergence might be
 * different compared to a real cardiac simulation. Do not take conclusions about preconditioner/solver performance
 * based on these tests only.
 */
class TestChebyshevIteration : public CxxTest::TestSuite
{
public:

    void TestChebyshevVsCG() throw (Exception)
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        unsigned cg_its;
        unsigned chebyshev_its;

        Timer::Reset();
        {
            Mat system_matrix;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

            Vec system_rhs;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetMatrixIsSymmetric();
            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");
            ls.SetPcType("bjacobi");

            Vec solution = ls.Solve();

            cg_its = ls.GetNumIterations();

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::PrintAndReset("CG");

        {
            Mat system_matrix;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

            Vec system_rhs;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetMatrixIsSymmetric();
            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("chebychev");
            ls.SetPcType("bjacobi");

            Vec solution = ls.Solve();

            chebyshev_its = ls.GetNumIterations();

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::Print("Chebyshev");

        TS_ASSERT_LESS_THAN(cg_its, 15u); // Takes 14 iterations with 16 cores
        TS_ASSERT_LESS_THAN(chebyshev_its, 17u); // Takes 16 iterations with 16 cores

        VecDestroy(parallel_layout);
    }

    void TestChebyshevAdaptiveVsNoAdaptive() throw (Exception)
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        // Solving with zero guess for coverage
        Vec zero_guess = factory.CreateVec(2);
        double zero = 0.0;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2)
        VecSet(&zero, zero_guess);
#else
        VecSet(zero_guess, zero);
#endif

        Mat system_matrix;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        Vec system_rhs;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

        // Make sure we are not inheriting a non-default number of iterations from previous test
        std::stringstream num_it_str;
        num_it_str << 1000;
        PetscOptionsSetValue("-ksp_max_it", num_it_str.str().c_str());

        try
        {
            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetMatrixIsSymmetric();
            // Solve to relative convergence for coverage
            ls.SetRelativeTolerance(1e-6);
            ls.SetPcType("jacobi");
            ls.SetKspType("chebychev");
            ls.SetUseFixedNumberIterations(true, 64);

            // Solving with zero guess for coverage.
            Vec solution = ls.Solve(zero_guess);
            unsigned chebyshev_adaptive_its = ls.GetNumIterations();

            TS_ASSERT_EQUALS(chebyshev_adaptive_its, 40u);
            TS_ASSERT_DELTA(ls.mEigMin, 0.0124, 1e-4);
            TS_ASSERT_DELTA(ls.mEigMax, 1.8810, 1e-4);

            VecDestroy(solution);
        }
        catch (Exception& e)
        {
            if (e.GetShortMessage() == "Chebyshev with fixed number of iterations is known to be broken in PETSc <= 2.3.2")
            {
                WARNING(e.GetShortMessage());
            }
            else
            {
                TS_FAIL(e.GetShortMessage());
            }
        }


        // Make sure we are not inheriting a non-default number of iterations from previous test
        PetscOptionsSetValue("-ksp_max_it", num_it_str.str().c_str());
        {
            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetMatrixIsSymmetric();
            ls.SetRelativeTolerance(1e-6);
            ls.SetPcType("jacobi");
            ls.SetKspType("chebychev");

            Vec solution = ls.Solve(zero_guess);
            unsigned chebyshev_no_adaptive_its = ls.GetNumIterations();

            TS_ASSERT_LESS_THAN(chebyshev_no_adaptive_its, 100u); // Normally 88, but 99 on maverick & natty
            TS_ASSERT_DELTA(ls.mEigMin, 0.0124, 1e-4);
            TS_ASSERT_DELTA(ls.mEigMax, 1.8841, 1e-4);

            VecDestroy(solution);
        }

        // Make sure we are not inheriting a non-default number of iterations from previous test
        PetscOptionsSetValue("-ksp_max_it", num_it_str.str().c_str());
        {
            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetMatrixIsSymmetric();
            ls.SetRelativeTolerance(1e-6);
            ls.SetPcType("jacobi");
            ls.SetKspType("cg");
            Vec solution = ls.Solve(zero_guess);
            unsigned cg_its = ls.GetNumIterations();

            TS_ASSERT_EQUALS(cg_its, 40u);
            TS_ASSERT_EQUALS(ls.mEigMin, DBL_MAX);
            TS_ASSERT_EQUALS(ls.mEigMax, DBL_MIN);

            VecDestroy(solution);
        }

        MatDestroy(system_matrix);
        VecDestroy(system_rhs);

        VecDestroy(parallel_layout);
        VecDestroy(zero_guess);
    }
};

#endif /*TESTPCBLOCKDIAGONAL_HPP_*/
