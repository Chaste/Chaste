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

#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
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

        VecDestroy(vec1);
        VecDestroy(vec2);
        MatDestroy(mat);

        // Test SetupMatrix with non-default preallocation
        Mat mat2;
        PetscTools::SetupMat(mat2, 12, 10, 4);
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
            // plus 0.5*4*12 in the off-diagonal part. These are then split between the number of processors. So, a
            // processor that owns n rows should have 6*n nonzeros allocated.
            PetscInt lo, hi;
            MatGetOwnershipRange(mat2, &lo, &hi);
            TS_ASSERT_EQUALS( nonzeros_allocated, (unsigned)(6*(hi-lo)) );
        }

        MatDestroy(mat2);

        Mat mat_over_allocate;
        PetscTools::SetupMat(mat_over_allocate, 12, 12, 13);
        MatDestroy(mat_over_allocate);

        // coverage
        Mat mat3;
        PetscTools::SetupMat(mat3, 1, 1, 0);
        MatDestroy(mat3);
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
        PetscTools::IsolateProcesses();
        TS_ASSERT(PetscTools::AmMaster());
        // Note: this will deadlock in parallel if IsolateProcesses doesn't work
        if (PetscTools::AmTopMost())
        {
            PetscTools::Barrier("TestProcessIsolation");
        }
        bool am_top = PetscTools::AmTopMost();
        bool any_is_top = PetscTools::ReplicateBool(am_top);
        TS_ASSERT_EQUALS(am_top, any_is_top);
        if (PetscTools::AmTopMost())
        {
            TS_ASSERT_THROWS_NOTHING(PetscTools::ReplicateException(true));
        }
        else
        {
            TS_ASSERT_THROWS_NOTHING(PetscTools::ReplicateException(false));
        }
        PetscTools::IsolateProcesses(false);
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

        MatDestroy(matrix);
        VecDestroy(vector);

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

        MatDestroy(matrix_read);
        VecDestroy(vector_read);
    }

    /*
     * This test reuses the 10x10 matrix written to disc in the previous test. It reads it
     * back in with a different parallel layout. For p=2 it is partitioned in 6 and 4 rows,
     * for p=3 4, 4, and 2.
     */
    void TestReadWithNonDefaultParallelLayout()
    {
        DistributedVectorFactory factory(5);
        Vec parallel_layout = factory.CreateVec(2);

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

        MatDestroy(matrix_read);
        VecDestroy(vector_read);
        VecDestroy(parallel_layout);
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

        VecDestroy(petsc_vec_uneven);
    }
};

#endif /*TESTPETSCTOOLS_HPP_*/
