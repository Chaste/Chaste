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

#ifndef TESTPETSCTOOLS_HPP_
#define TESTPETSCTOOLS_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>
#include <cstring>
#include "DistributedVectorFactory.hpp"
#include "ReplicatableVector.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "OutputFileHandler.hpp"
#include "DistributedVector.hpp"

class TestPetscTools : public CxxTest::TestSuite
{
public:

    void TestMostOfPetscTools()
    {
        TS_ASSERT(PetscTools::IsInitialised());
        PetscInt my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        TS_ASSERT_EQUALS(PetscTools::GetMyRank(), (unsigned)my_rank);
        bool am_master = (my_rank == 0);
        TS_ASSERT_EQUALS( PetscTools::AmMaster(), am_master);

        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        TS_ASSERT_EQUALS( PetscTools::GetNumProcs(), (unsigned)num_procs);
        bool is_sequential = (num_procs==1);
        TS_ASSERT_EQUALS( PetscTools::IsSequential(), is_sequential);
        TS_ASSERT_EQUALS( PetscTools::IsParallel(), !is_sequential);
        bool am_right = (my_rank == num_procs - 1 );
        TS_ASSERT_EQUALS( PetscTools::AmTopMost(), am_right);

        std::cout << "These should be ordered:" << std::endl;
        PetscTools::BeginRoundRobin();
        std::cout << "  Process " << PetscTools::GetMyRank() << std::endl;
        PetscTools::EndRoundRobin();

        // Test CreateVec which returns a vec of constants
        Vec vec1 = PetscTools::CreateAndSetVec(10, 3.41);
        ReplicatableVector vec1_repl(vec1);

        TS_ASSERT_EQUALS(vec1_repl.GetSize(), 10u);
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(vec1_repl[i], 3.41, 1e-12);
        }

        // Test CreateVec which uses a std::vector of data
        std::vector<double> data(10);
        for (unsigned i=0; i<10; i++)
        {
            data[i] = i+0.45;
        }

        Vec vec2 = PetscTools::CreateVec(data);

        ReplicatableVector vec2_repl(vec2);

        TS_ASSERT_EQUALS(vec2_repl.GetSize(), 10u);
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(vec2_repl[i], i+0.45, 1e-12);
        }

        // Test SetupMatrix
        Mat mat;
        PetscTools::SetupMat(mat, 10, 11, 11);
        int m,n;
        MatGetSize(mat, &m, &n);
        TS_ASSERT_EQUALS(m, 10);
        TS_ASSERT_EQUALS(n, 11);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 3) //PETSc 3.0 to PETSc 3.3
        //The PETSc developers changed this one, but later changed it back again!
        const MatType type;
#else
        MatType type;
#endif
        MatGetType(mat,&type);
        if (PetscTools::IsSequential())
        {
            TS_ASSERT(strcmp(type, MATSEQAIJ)==0);
        }
        else
        {
            TS_ASSERT(strcmp(type, MATMPIAIJ)==0);
        }

        PetscTools::Destroy(vec1);
        PetscTools::Destroy(vec2);
        PetscTools::Destroy(mat);

        // Test SetupMatrix with non-default preallocation
        Mat mat2;
        PetscTools::SetupMat(mat2, 12, 10, 4, PETSC_DECIDE, PETSC_DECIDE, false, false);
        MatGetSize(mat2, &m, &n);
        TS_ASSERT_EQUALS(m, 12);
        TS_ASSERT_EQUALS(n, 10);

        MatGetType(mat2,&type);
        if (PetscTools::IsSequential())
        {
            TS_ASSERT(strcmp(type, MATSEQAIJ)==0);
        }
        else
        {
            TS_ASSERT(strcmp(type, MATMPIAIJ)==0);
        }

        MatInfo info;
        unsigned nonzeros_allocated;

        MatGetInfo(mat2,MAT_LOCAL,&info);
        nonzeros_allocated = (unsigned) info.nz_allocated;

        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS( nonzeros_allocated, 4*12u );
        }
        else
        {
            // Total number of nozeros that should be allocated is 36 (4*12 (12 = number of rows) in diagonal part,
            // plus 4*12 in the off-diagonal part. These are then split between the number of processors. So, a
            // processor that owns n rows should have 8*n nonzeros allocated.
            PetscInt lo, hi;
            MatGetOwnershipRange(mat2, &lo, &hi);
            TS_ASSERT_EQUALS( nonzeros_allocated, (unsigned)(2*4*(hi-lo)) );
        }

        PetscTools::Destroy(mat2);

        Mat mat_over_allocate;
        PetscTools::SetupMat(mat_over_allocate, 12, 12, 13);
        PetscTools::Destroy(mat_over_allocate);

        // coverage
        Mat mat3;
        PetscTools::SetupMat(mat3, 1, 1, 0);
        PetscTools::Destroy(mat3);
    }

    void TestBarrier()
    {
        /*
         * Testing the barrier method is kind of tricky, since we really want to check
         * if it also works when PETSc isn't set up.  So see TestPetscTools2.hpp!
         */
        PetscTools::Barrier("TestBarrier");
    }

    void TestReplicateBool()
    {
        bool my_flag = false;
        if (PetscTools::AmMaster())
        {
            my_flag = true;
        }

        TS_ASSERT(PetscTools::ReplicateBool(my_flag));
    }

    void TestReplicateException()
    {
        DistributedVectorFactory factory(1);
        if (factory.IsGlobalIndexLocal(0))
        {
            TS_ASSERT_THROWS_NOTHING(PetscTools::ReplicateException(true));
        }
        else
        {
            TS_ASSERT_THROWS_THIS(PetscTools::ReplicateException(false), "Another process threw an exception; bailing out.");
        }
    }

    void TestProcessIsolation()
    {
        TS_ASSERT_EQUALS(PetscTools::GetWorld(), PETSC_COMM_WORLD); // No isolation at first
        TS_ASSERT(!PetscTools::IsIsolated());
        bool really_parallel = PetscTools::IsParallel();
        unsigned my_rank = PetscTools::GetMyRank();
        unsigned num_procs = PetscTools::GetNumProcs();
        PetscTools::IsolateProcesses();
        TS_ASSERT_EQUALS(PetscTools::GetWorld(), PETSC_COMM_SELF); // Each process is its own world
        TS_ASSERT(PetscTools::IsIsolated());
        TS_ASSERT(PetscTools::AmMaster());  // All processes are masters
        TS_ASSERT(PetscTools::AmTopMost()); // All processes are top
        TS_ASSERT_EQUALS(PetscTools::GetMyRank(), my_rank); // Rank and process count are still accurate though
        TS_ASSERT_EQUALS(PetscTools::GetNumProcs(), num_procs);
        // Note: this will deadlock in parallel if IsolateProcesses doesn't work
        if (PetscTools::GetMyRank() == 0u)
        {
            // Only the rank 0 process does this
            PetscTools::Barrier("TestProcessIsolation");
        }
        bool am_top = (PetscTools::GetMyRank() == PetscTools::GetNumProcs() - 1);
        bool any_is_top = PetscTools::ReplicateBool(am_top); // Replication is a no-op
        TS_ASSERT_EQUALS(am_top, any_is_top);
        if (PetscTools::GetMyRank() == 0u)
        {
            TS_ASSERT_THROWS_NOTHING(PetscTools::ReplicateException(true));
        }
        else
        {
            TS_ASSERT_THROWS_NOTHING(PetscTools::ReplicateException(false));
        }
        TS_ASSERT(PetscTools::IsSequential()); // IsParallel() will tell the truth, however
        TS_ASSERT_EQUALS(PetscTools::IsParallel(), really_parallel);
        PetscTools::IsolateProcesses(false);
        TS_ASSERT(!PetscTools::IsIsolated());
    }

    void TestDumpPetscObjects()
    {
        Mat matrix;
        Vec vector;

        PetscTools::SetupMat(matrix, 10, 10, 10);

        vector = PetscTools::CreateVec(10);

        PetscInt lo, hi;
        VecGetOwnershipRange(vector, &lo, &hi);

        for (int row=0; row<10; row++)
        {
            if (row >= lo && row < hi)
            {
                for (int col=0; col<10; col++)
                {
                    MatSetValue(matrix, row, col, (double) 10*row+col+1, INSERT_VALUES);
                }

                double value = row;
                VecSetValues(vector, 1, &row, &value, INSERT_VALUES);
            }
        }

        MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(vector);
        VecAssemblyEnd(vector);

        OutputFileHandler handler("DumpPetscObjects");
        std::string output_dir = handler.GetOutputDirectoryFullPath();

        PetscTools::DumpPetscObject(matrix, output_dir + "ten_times_ten.mat");
        PetscTools::DumpPetscObject(vector, output_dir + "ten_times_ten.vec");

        PetscTools::Destroy(matrix);
        PetscTools::Destroy(vector);

        Mat matrix_read;
        Vec vector_read;

        PetscTools::ReadPetscObject(matrix_read, output_dir + "ten_times_ten.mat");
        PetscTools::ReadPetscObject(vector_read, output_dir + "ten_times_ten.vec");

        double* p_vector_read;
        VecGetArray(vector_read, &p_vector_read);

        for (PetscInt row=0; row<10; row++)
        {
            if (lo<=row && row<hi)
            {
                for (PetscInt col=0; col<10; col++)
                {
                    double value;
                    MatGetValues(matrix_read, 1, &row, 1, &col, &value);
                    TS_ASSERT_EQUALS(value, (double) 10*row+col+1);
                }

            unsigned local_index = row-lo;
            TS_ASSERT_EQUALS(p_vector_read[local_index], (double)row);
            }
        }

        VecRestoreArray(vector_read, &p_vector_read);

        PetscTools::Destroy(matrix_read);
        PetscTools::Destroy(vector_read);
    }

    /*
     * This test reuses the 10x10 matrix written to disc in the previous test. It reads it
     * back in with a different parallel layout. For p=2 it is partitioned in 6 and 4 rows,
     * for p=3 4, 4, and 2.
     */
    void TestReadWithNonDefaultParallelLayout()
    {
        DistributedVectorFactory factory(5);


        //This is like factory.CreateVec(2) but doesn't change the block size
        Vec parallel_layout;
        unsigned small_locals =  factory.GetHigh() - factory.GetLow();
        VecCreateMPI(PETSC_COMM_WORLD, 2*small_locals, 10, &parallel_layout);

        PetscInt lo, hi;
        VecGetOwnershipRange(parallel_layout, &lo, &hi);

        Mat matrix_read;
        Vec vector_read;

        OutputFileHandler handler("DumpPetscObjects", false);
        std::string output_dir = handler.GetOutputDirectoryFullPath();

        PetscTools::ReadPetscObject(matrix_read, output_dir + "ten_times_ten.mat", parallel_layout);
        PetscTools::ReadPetscObject(vector_read, output_dir + "ten_times_ten.vec", parallel_layout);

        double* p_vector_read;
        VecGetArray(vector_read, &p_vector_read);

        for (PetscInt row=0; row<10; row++)
        {
            if (lo<=row && row<hi)
            {
                for (PetscInt col=0; col<10; col++)
                {
                    double value;
                    MatGetValues(matrix_read, 1, &row, 1, &col, &value);
                    TS_ASSERT_EQUALS(value, (double) 10*row+col+1);
                }

            unsigned local_index = row-lo;
            TS_ASSERT_EQUALS(p_vector_read[local_index], (double)row);
            }
        }

        VecRestoreArray(vector_read, &p_vector_read);

        PetscTools::Destroy(matrix_read);
        PetscTools::Destroy(vector_read);
        PetscTools::Destroy(parallel_layout);
    }

    void TestUnevenCreation()
    {
        /*
         * Uneven test (as in TestDistributedVectorFactory).
         * Calculate total number of elements in the vector.
         */
        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned total_elements = (num_procs+1)*num_procs/2;
        unsigned my_rank=PetscTools::GetMyRank();

        Vec petsc_vec_uneven = PetscTools::CreateVec(total_elements, my_rank+1);

        int petsc_lo, petsc_hi;
        VecGetOwnershipRange(petsc_vec_uneven, &petsc_lo, &petsc_hi);

        unsigned expected_lo = (my_rank+1)*my_rank/2;
        unsigned expected_hi = (my_rank+2)*(my_rank+1)/2;

        TS_ASSERT_EQUALS((unsigned)petsc_lo, expected_lo);
        TS_ASSERT_EQUALS((unsigned)petsc_hi, expected_hi);

        PetscTools::Destroy(petsc_vec_uneven);
    }

    void TestHasParMetis()
    {
        //This just covers the method, as there is no other way to test if ParMetis is available.
        std::cout << "Testing to see if Petsc is configured with ParMetis support. " << std::endl;
        PetscTools::HasParMetis();
    }

    void TestExceptionMacros()
    {
        // Should not throw
        TS_ASSERT_THROWS_NOTHING(TRY_IF_MASTER(std::cout << "No exception\n"));

        // Should not throw - since only master should run the contents of the macro
        TS_ASSERT_THROWS_NOTHING(TRY_IF_MASTER(if (!PetscTools::AmMaster()) EXCEPTION("Should not occur")));

        // Should get thrown by all processes as exception is replicated.
        TS_ASSERT_THROWS_CONTAINS(TRY_IF_MASTER(EXCEPTION("master; bailing out")),
                "; bailing out"); // both the replicated and original should contain this phrase
    }

    void TestSetOptionWithLogging()
    {
        // See #2933: we need to cover either PetscLogBegin() or PetscLogDefaultBegin()
        PetscTools::SetOption("-log_summary", "");
    }
};

#endif /*TESTPETSCTOOLS_HPP_*/
