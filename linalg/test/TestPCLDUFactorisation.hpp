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
class TestPCLDUFactorisation : public CxxTest::TestSuite
{
public:

    void TestBasicFunctionality()
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

        PetscTools::Destroy(parallel_layout);

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
        PetscTools::Destroy(one_zeros);

        LinearSystem ls = LinearSystem(rhs, system_matrix);

        ls.SetAbsoluteTolerance(1e-9);
        ls.SetKspType("cg");
        ls.SetPcType("ldufactorisation");

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
        ls.SetPcType("ldufactorisation");

        PetscTools::Destroy(system_matrix);
        PetscTools::Destroy(rhs);
        PetscTools::Destroy(solution);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 3) //PETSc 3.0 to PETSc 3.3
        //The PETSc developers changed this one, but later changed it back again!
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

            PetscTools::Destroy(system_matrix);
            PetscTools::Destroy(system_rhs);
            PetscTools::Destroy(solution);
        }
        Timer::PrintAndReset("No preconditioning");

        {
            Mat system_matrix;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

            Vec system_rhs;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");
            ls.SetPcType("ldufactorisation");

            Vec solution = ls.Solve();

            block_diag_its = ls.GetNumIterations();

            PetscTools::Destroy(system_matrix);
            PetscTools::Destroy(system_rhs);
            PetscTools::Destroy(solution);
        }
        Timer::Print("LDU factorisation preconditioner");

        std::cout << block_diag_its << " " << point_jacobi_its << std::endl;
        TS_ASSERT_LESS_THAN_EQUALS(block_diag_its, point_jacobi_its);

        PetscTools::Destroy(parallel_layout);
    }
};

#endif /*TESTPCBLOCKDIAGONAL_HPP_*/
