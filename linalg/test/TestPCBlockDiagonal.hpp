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

#ifndef TESTPCBLOCKDIAGONAL_HPP_
#define TESTPCBLOCKDIAGONAL_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "Timer.hpp"
#include "DistributedVectorFactory.hpp"
#include <cstring>


/*
 * Warning: these tests do not inform PETSc about the nullspace of the matrix. Therefore, convergence might be
 * different compared to a real cardiac simulation. Do not take conclusions about preconditioner performance
 * based on these tests only.
 */
class TestPCBlockDiagonal : public CxxTest::TestSuite
{
public:

    void TestBasicFunctionality() throw (Exception)
    {
        /*
         * We need to make sure here that the matrix is loaded with the appropriate parallel layout. Petsc's
         * default puts 1331 rows in each processor. This wouldn't be possible in a real bidomain simulation
         * because implies that equations V_665 an Phi_e_665 are solved in different processors.
         */
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        Mat system_matrix;
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        VecDestroy(parallel_layout);

        // Set rhs = A * [1 0 1 0 ... 1 0]'
        Vec one_zeros = factory.CreateVec(2);
        Vec rhs = factory.CreateVec(2);

        for (unsigned node_index=0; node_index<2*num_nodes; node_index+=2)
        {
            PetscVecTools::SetElement(one_zeros, node_index, 1.0);
            PetscVecTools::SetElement(one_zeros, node_index+1, 0.0);
        }
        PetscVecTools::Finalise(one_zeros);

        MatMult(system_matrix, one_zeros, rhs);
        VecDestroy(one_zeros);

        LinearSystem ls = LinearSystem(rhs, system_matrix);

        ls.SetAbsoluteTolerance(1e-9);
        ls.SetKspType("cg");
        ls.SetPcType("none");

        ls.AssembleFinalLinearSystem();

        Vec solution = ls.Solve();

        DistributedVector distributed_solution = factory.CreateDistributedVector(solution);
        DistributedVector::Stripe vm(distributed_solution, 0);
        DistributedVector::Stripe phi_e(distributed_solution, 1);

        for (DistributedVector::Iterator index = distributed_solution.Begin();
             index!= distributed_solution.End();
             ++index)
        {
            /*
             * Although we're trying to enforce the solution to be [1 0 ... 1 0], the system is singular and
             * therefore it has infinite solutions. I've (migb) found that the use of different preconditioners
             * lead to different solutions ([0.8 -0.2 ... 0.8 -0.2], [0.5 -0.5 ... 0.5 -0.5], ...)
             *
             * If we were using PETSc null space, it would find the solution that satisfies x'*v=0,
             * being v the null space of the system (v=[1 1 ... 1])
             */
            TS_ASSERT_DELTA(vm[index] - phi_e[index], 1.0, 1e-6);
        }

        // Coverage (setting PC type after first solve)
        ls.SetPcType("blockdiagonal");

        MatDestroy(system_matrix);
        VecDestroy(rhs);
        VecDestroy(solution);

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        const PCType pc;
#else
        PCType pc;
#endif
        PC prec;
        KSPGetPC(ls.mKspSolver, &prec);
        PCGetType(prec, &pc);
        // Although we call it "blockdiagonal", PETSc considers this PC a generic SHELL preconditioner
        TS_ASSERT( strcmp(pc,"shell")==0 );

    }

    void TestBetterThanNoPreconditioning()
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        unsigned point_jacobi_its;
        unsigned block_diag_its;

        Timer::Reset();
        {
            Mat system_matrix;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

            Vec system_rhs;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");
            ls.SetPcType("none");

            Vec solution = ls.Solve();

            point_jacobi_its = ls.GetNumIterations();

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::PrintAndReset("No preconditioning");

        {
            Mat system_matrix;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

            Vec system_rhs;
            // Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");
            ls.SetPcType("blockdiagonal");

            Vec solution = ls.Solve();

            block_diag_its = ls.GetNumIterations();

            // Coverage (setting PC type after using blockdiagonal solve)
            ls.SetPcType("blockdiagonal");

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::Print("Block diagonal preconditioner");

        std::cout << block_diag_its << " " << point_jacobi_its << std::endl;
        TS_ASSERT_LESS_THAN_EQUALS(block_diag_its, point_jacobi_its);

        VecDestroy(parallel_layout);
    }
};

#endif /*TESTPCBLOCKDIAGONAL_HPP_*/
