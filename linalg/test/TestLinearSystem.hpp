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

#ifndef _TESTLINEARSYSTEM_HPP_
#define _TESTLINEARSYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <cmath>
#include <iostream>
#include <cstring>

#include "LinearSystem.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"
#include "OutputFileHandler.hpp"
#include "DistributedVectorFactory.hpp"
#include "ReplicatableVector.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Timer.hpp"

/**
 * Tests the LinearSystem class, and some methods in the PETSc helper classes PetscVecTools and PetscMatTools.
 */
class TestLinearSystem : public CxxTest::TestSuite
{
public:

   void TestLinearSystem1()
    {
        TS_ASSERT_THROWS_THIS(LinearSystem too_big_to_be_dense(20), "You must provide a rowPreallocation argument for a large sparse system");

        const unsigned size_u = 3u;
        const int size = (int) size_u;
        LinearSystem ls(size_u);
        TS_ASSERT_EQUALS(ls.GetSize(), size_u);
        TS_ASSERT_EQUALS(PetscVecTools::GetSize(ls.GetRhsVector()), size_u);
        TS_ASSERT_EQUALS(PetscMatTools::GetSize(ls.GetLhsMatrix()), size_u);

        ls.SetMatrixIsConstant(true);

        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double) row*3+col+1);
            }
        }

        ls.SetRhsVectorElement(0, 1400000.0);
        ls.SetRhsVectorElement(1, 3200000.0);
        ls.SetRhsVectorElement(2, 5000000.0);

        ls.AssembleFinalLinearSystem();

        // For coverage
        ls.DisplayMatrix();
        ls.DisplayRhs();

        Vec solution_vector;
        solution_vector = ls.Solve();

        int lo, hi;
        VecGetOwnershipRange(solution_vector, &lo, &hi);
        PetscScalar* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);
        for (int global_index=0; global_index<size; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 100000.0*(global_index+1), 1e-8);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        PetscTools::Destroy(solution_vector);

        // SetRelativeTolerance
        ls.SetRelativeTolerance(1e-2);
        TS_ASSERT_THROWS_NOTHING(solution_vector = ls.Solve());

        KSPConvergedReason reason;
        KSPGetConvergedReason(ls.mKspSolver, &reason);
        TS_ASSERT_EQUALS(reason, KSP_CONVERGED_RTOL);

        VecGetOwnershipRange(solution_vector,&lo,&hi);
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<size; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 100000.0*(global_index+1), 2e4);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        PetscTools::Destroy(solution_vector);

        /*
         * Reset KSP stuff. This doesn't need to be done, but we're making sure
         * we cover this method (and check that it works).
         */
        ls.ResetKspSolver();

        // SetAbsoluteTolerance
        ls.SetAbsoluteTolerance(1e-8);
        solution_vector = ls.Solve();
        KSPGetConvergedReason(ls.mKspSolver, &reason);
        TS_ASSERT_EQUALS(reason, KSP_CONVERGED_ATOL);

        // Check that it converged for the right reason
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<size; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 100000.0*(global_index+1), 1e-8);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);
        PetscTools::Destroy(solution_vector);
    }

    void TestZeroingLinearSystem()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
            ls.SetRhsVectorElement(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for (int row=lo; row<hi; row++)
        {
            TS_ASSERT_EQUALS(ls.GetRhsVectorElement(row), row);
            for (int col=0; col<5; col++)
            {
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), row);
            }
        }

        ls.ZeroLinearSystem();

        for (int row=lo; row<hi; row++)
        {
            TS_ASSERT_EQUALS(ls.GetRhsVectorElement(row), 0.0);
            for (int col=0; col<5; col++)
            {
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 0.0);
            }
        }


        ls.SetMatrixRow(2, 125.0);

        // Note: this method is collective. All processes MUST do it together.
        ls.AssembleFinalLinearSystem();

        if (lo <=2 && 2<hi)
        {
            TS_ASSERT_EQUALS(ls.GetMatrixElement(2, 1), 125.0);
        }
    }

    void TestZeroMatrixRowsWithValueOnDiagonal()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
            ls.SetRhsVectorElement(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        std::vector<unsigned> rows;
        rows.push_back(2);
        rows.push_back(3);
        rows.push_back(4);

        ls.ZeroMatrixRowsWithValueOnDiagonal(rows, 3.14);

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for (int row=2; row<5; row++)
        {
            if (lo<=row && row<hi)
            {
                for (int i=0; i<5; (i+1==row? i+=2 : i++)) // for i=0,1..,row-1,row+1,..,5
                {
                    TS_ASSERT_EQUALS(ls.GetMatrixElement(row,i), 0.0);
                }
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row,row), 3.14);
            }
        }
    }

    void TestAddingNonzeroesLater()
    {
        // Make a linear system from which the vector becomes a template to use in a later constructor
        LinearSystem ls_template(5);
        for (int i=0; i<5; i++)
        {
            ls_template.SetMatrixElement(i, i, 3.0);
        }
        ls_template.SetMatrixElement(1, 3, 3.0);
        ls_template.AssembleFinalLinearSystem();

        // The "false" says that we are allowed to do new mallocs without PETSc 3.3 causing an error
        LinearSystem ls(ls_template.rGetRhsVector(), 1, false);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixElement(i, i, 3.0);
        }
        ls.AssembleFinalLinearSystem();
        //Adding a new non-zero element gives an error in PETSc 3.3
        ls.SetMatrixElement(1, 3, 3.0);
        ls.AssembleFinalLinearSystem();

        TS_ASSERT(PetscMatTools::CheckEquality(ls_template.rGetLhsMatrix(), ls.rGetLhsMatrix()));
    }

    void TestZeroingLinearSystemByColumn()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixElement(i, i, 3.0);
        }
        ls.SetMatrixElement(0, 1, 4.0);

        ls.AssembleFinalLinearSystem();

        for (unsigned col=0; col<5; col++)
        {
            ls.ZeroMatrixColumn(col);
        }
        ls.AssembleFinalLinearSystem();

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for (int row=lo; row<hi; row++)
        {
            TS_ASSERT_EQUALS(ls.GetRhsVectorElement(row), 0);
            for (int col=0; col<5; col++)
            {
                TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 0);
            }
        }

        MatInfo info;
        double num_nonzeros;

        MatGetInfo(ls.rGetLhsMatrix(),MAT_GLOBAL_SUM,&info);

        num_nonzeros = info.nz_used;

        TS_ASSERT_EQUALS(int(num_nonzeros),6);
    }

    void TestZeroMatrixRowsAndColumnsWithValueOnDiagonal()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        std::vector<unsigned> rows(3);
        rows[0] = 2;
        rows[1] = 3;
        rows[2] = 4;

        ls.ZeroMatrixRowsAndColumnsWithValueOnDiagonal(rows, 3.1);

        int lo, hi;
        ls.GetOwnershipRange(lo, hi);
        for (int row=lo; row<hi; row++)
        {
            for (int col=0; col<5; col++)
            {
                if ((col>=2) || (row>=2))
                {
                    // The altered values
                    if (row != col)
                    {
                        TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 0);
                    }
                    else
                    {
                        TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), 3.1);
                    }
                }
                else
                {
                    // Unaltered values
                    TS_ASSERT_EQUALS(ls.GetMatrixElement(row, col), row);
                }
            }
        }
    }

    void TestGetMatrixRowDistributed()
    {
        LinearSystem ls(5);
        for (int i=0; i<5; i++)
        {
            ls.SetMatrixRow(i, (double)i);
        }
        ls.AssembleFinalLinearSystem();

        Vec third_row = ls.GetMatrixRowDistributed(3);

        DistributedVectorFactory factory(third_row);

        DistributedVector distributed_third_row = factory.CreateDistributedVector(third_row);
        for (DistributedVector::Iterator index = distributed_third_row.Begin();
             index!= distributed_third_row.End();
             ++index)
        {
            TS_ASSERT_EQUALS(distributed_third_row[index], 3);
        }

        PetscTools::Destroy(third_row);
    }


    void TestCreateFromVector()
    {
        const int SIZE = 5;
        Vec test_vec=PetscTools::CreateVec(SIZE);

        LinearSystem ls(test_vec, 5);

        // Check ownership ranges match
        int lo1, hi1, lo2, hi2;
        VecGetOwnershipRange(test_vec, &lo1, &hi1);
        ls.GetOwnershipRange(lo2, hi2);
        TS_ASSERT_EQUALS(lo1, lo2);
        TS_ASSERT_EQUALS(hi1, hi2);

        PetscTools::Destroy(test_vec);
    }

    void TestLinearSystem2()
    {
        LinearSystem ls(2);
        ls.SetMatrixRow(0, 1.0);
        ls.SetMatrixRow(1, 3.0);
        ls.AssembleIntermediateLinearSystem();

        ls.AddToMatrixElement(0, 1, 1.0);
        ls.AddToMatrixElement(1, 1, 1.0);

        ls.AddToRhsVectorElement(0, 3.0);
        ls.AddToRhsVectorElement(1, 7.0);
        ls.AssembleFinalLinearSystem();

        Vec solution_vector;
        solution_vector = ls.Solve();


        int lo,hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<2; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], 1.0, 0.000001);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);

        PetscTools::Destroy(solution_vector);
    }

    void TestLinearSystem3()
    {
        LinearSystem ls(7);
        ls.SetMatrixIsConstant(true);

        for (unsigned i = 0; i < 7; ++i)
        {
            ls.SetMatrixElement(i, i, static_cast<double>(i + 1));
        }

        ls.AssembleIntermediateLinearSystem();

        ls.AddToMatrixElement(3, 1, sqrt(5.0));

        /*
         * Matrix looks like:
         *  1 0 0 0 0 0 0
         *  0 2 0 0 0 0 0
         *  0 0 3 0 0 0 0
         *  0 x 0 4 0 0 0
         *  0 0 0 0 5 0 0
         *  0 0 0 0 0 6 0
         *  0 0 0 0 0 0 7  where x = sqrt(5)
         */

        ls.SetRhsVectorElement(0, 1.23);
        ls.SetRhsVectorElement(1, 2.34);
        ls.SetRhsVectorElement(2, 3.45);
        ls.SetRhsVectorElement(3, 4.56);
        ls.SetRhsVectorElement(4, 5.67);
        ls.SetRhsVectorElement(5, 6.78);
        ls.SetRhsVectorElement(6, 7.89);
        ls.AssembleFinalLinearSystem();

        Vec solution_vector;
        solution_vector = ls.Solve();

        std::vector<double> hand_cooked_solution(7, 0.0);
        hand_cooked_solution[0] = 1.23;
        hand_cooked_solution[1] = 1.17;
        hand_cooked_solution[2] = 1.15;
        hand_cooked_solution[3] = 0.485950116581311;
        hand_cooked_solution[4] = 1.134;
        hand_cooked_solution[5] = 1.13;
        hand_cooked_solution[6] = 1.127142857142857;

        int lo;
        int hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<7; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], hand_cooked_solution[global_index], 1e-8);
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);

        PetscTools::Destroy(solution_vector);
    }

    /**
     * This is a stub test for coverage purposes.
     */
    void TestNullBasis()
    {
 #ifndef NDEBUG //Only do this test in debug mode, since functionality is skipped in optimized code
        unsigned size = 10;

        // Test it throws if one of the vectors in the base is not normal
        {
            std::vector<double> data(size,1.0);
            Vec non_orthonormal = PetscTools::CreateVec(data);

            LinearSystem ls((PetscInt) size);
            TS_ASSERT_THROWS_THIS(ls.SetNullBasis(&non_orthonormal, 1),
                    "One of the vectors in the null space is not normalised");
            PetscTools::Destroy(non_orthonormal);
        }

        // Test it throws if the vectors in the base are not orthogonal
        {
            std::vector<double> data(size,0.0);
            data[0] = 1.0;
            Vec one_zeros = PetscTools::CreateVec(data);

            std::vector<double> data2(size,0.0);
            data2[1] = 1.0;
            Vec zero_one_zeros = PetscTools::CreateVec(data2);

            Vec null_basis[] = {one_zeros, zero_one_zeros, one_zeros};

            LinearSystem ls((PetscInt) size);
            TS_ASSERT_THROWS_THIS(ls.SetNullBasis(null_basis, 3),"The null space is not orthogonal.");

            PetscTools::Destroy(one_zeros);
            PetscTools::Destroy(zero_one_zeros);
        }


        // Test it doesn't throws if a non-orthonormal basis is passed
        {
            std::vector<double> data(size,0.0);
            data[0] = 1.0;
            Vec one_zeros = PetscTools::CreateVec(data);

            std::vector<double> data2(size,0.0);
            data2[1] = 1.0;
            Vec zero_one_zeros = PetscTools::CreateVec(data2);

            Vec null_basis[] = {one_zeros, zero_one_zeros};

            LinearSystem ls((PetscInt) size);
            TS_ASSERT_THROWS_NOTHING(ls.SetNullBasis(null_basis, 2));

            PetscTools::Destroy(one_zeros);
            PetscTools::Destroy(zero_one_zeros);
        }
#endif
    }

   void TestRemoveNullSpace()
    {
        LinearSystem ls(3);
        ls.SetMatrixIsConstant(true);

        TS_ASSERT_EQUALS(ls.GetSize(), 3U);

        for (int row=0; row<3; row++)
        {
            ls.SetMatrixElement(row, row, (double) row+1);
        }

        ls.SetRhsVectorElement(0, 1.0);
        ls.SetRhsVectorElement(1, 2.0);
        ls.SetRhsVectorElement(2, 3.0);
        ls.AssembleFinalLinearSystem();

        std::vector<double> data(3,0.0);
        data[0] = 1.0;
        Vec one_zeros = PetscTools::CreateVec(data);

        Vec null_basis[] = {one_zeros};

        ls.SetNullBasis(null_basis, 1);

        Vec wrong_solution = ls.Solve();
        ReplicatableVector replicated_wrong_solution(wrong_solution);

        // Wrong solution since wrong null space was provided.
        TS_ASSERT_DELTA(replicated_wrong_solution[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(replicated_wrong_solution[1], 1.0, 1e-8);
        TS_ASSERT_DELTA(replicated_wrong_solution[2], 1.0, 1e-8);

        // Now remove the null space and we will hopefully get the right solution
        ls.RemoveNullSpace();
        Vec solution = ls.Solve();
        ReplicatableVector replicated_solution(solution);

        TS_ASSERT_DELTA(replicated_solution[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(replicated_solution[1], 1.0, 1e-8);
        TS_ASSERT_DELTA(replicated_solution[2], 1.0, 1e-8);

        PetscTools::Destroy(one_zeros);
        PetscTools::Destroy(wrong_solution);
        PetscTools::Destroy(solution);

    }

    // Test the 3rd constructor
    void TestCreateWithPetscObjects()
    {
        // First try giving it just a vector
        unsigned size = 5u;
        DistributedVectorFactory factory(size);
        Vec test_vec = factory.CreateVec();

        DistributedVector dist_vec = factory.CreateDistributedVector(test_vec);
        double test_val = -1.0;
        if (dist_vec.Begin() != dist_vec.End())
        {
            dist_vec[dist_vec.Begin()] = test_val;
        }
        dist_vec.Restore();

        LinearSystem lsv(test_vec, (Mat) NULL);
        TS_ASSERT_EQUALS(lsv.GetSize(), size);

        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsv.GetRhsVectorElement(dist_vec.Begin().Global),
                             test_val);
        }

        // Change the Vec and see if the linear system reflects the change
        double test_val2 = 2.0;
        if (dist_vec.Begin() != dist_vec.End())
        {
            dist_vec[dist_vec.Begin()] = test_val2;
        }
        dist_vec.Restore();
        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsv.GetRhsVectorElement(dist_vec.Begin().Global),
                             test_val2);
        }

        // Now try with just a matrix
        Mat m;
        PetscTools::SetupMat(m, size, size, size);

        if (dist_vec.Begin() != dist_vec.End())
        {
            PetscMatTools::SetElement(m, dist_vec.Begin().Global, 0, test_val);
        }

        LinearSystem lsm(NULL, m);
        TS_ASSERT_EQUALS(lsm.GetSize(), size);
        lsm.FinaliseLhsMatrix();

        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsm.GetMatrixElement(dist_vec.Begin().Global, 0),
                             test_val);
        }

        // Change the Mat and see if the linear system reflects the change
        if (dist_vec.Begin() != dist_vec.End())
        {
            PetscMatTools::SetElement(m, dist_vec.Begin().Global, 0, test_val2);
        }
        MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
        if (dist_vec.Begin() != dist_vec.End())
        {
            TS_ASSERT_EQUALS(lsm.GetMatrixElement(dist_vec.Begin().Global, 0),
                             test_val2);
        }

        PetscTools::Destroy(test_vec);
        PetscTools::Destroy(m);
    }

    void TestLinearSystem1WithIntialGuess()
    {
        LinearSystem ls(3);

        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double) row*3+col+1);
            }
        }

        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);
        ls.AssembleFinalLinearSystem();
//      /* This will mess up the following because the *first call* to Solve has to contain an initial guess. */
//        Vec sol1 = ls.Solve();
//        PetscTools::Destroy(sol1);
//        unsigned num_iters_from_default = ls.GetNumIterations();

        Vec zero_guess=PetscTools::CreateVec(3);
        PetscVecTools::SetElement(zero_guess, 0, 0.0);
        PetscVecTools::SetElement(zero_guess, 1, 0.0);
        PetscVecTools::SetElement(zero_guess, 2, 0.0);
        PetscVecTools::Finalise(zero_guess);
        Vec sol2 = ls.Solve(zero_guess);
        PetscTools::Destroy(zero_guess);
        PetscTools::Destroy(sol2);

        unsigned num_iters_from_zero = ls.GetNumIterations();

        // Set the correct answer for the intial guess
        Vec good_guess=PetscTools::CreateVec(3);
        PetscVecTools::SetElement(good_guess, 0, 1.0);
        PetscVecTools::SetElement(good_guess, 1, 2.0);
        PetscVecTools::SetElement(good_guess, 2, 3.0);
        PetscVecTools::Finalise(good_guess);


        Vec solution_vector;
        solution_vector = ls.Solve(good_guess);
        unsigned num_iters_from_perfect = ls.GetNumIterations();
        TS_ASSERT_EQUALS(num_iters_from_perfect, 0u);
        TS_ASSERT_LESS_THAN(num_iters_from_perfect, num_iters_from_zero);

        int lo, hi;
        VecGetOwnershipRange(solution_vector,&lo,&hi);
        PetscScalar* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);

        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_EQUALS(p_solution_elements_array[local_index], global_index+1.0);
                // Zero tolerance
            }
        }
        VecRestoreArray(solution_vector, &p_solution_elements_array);

        // Set the a bad intial guess
        Vec bad_guess;
        VecDuplicate(good_guess, &bad_guess);
        PetscScalar too_big = 1e5;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        VecSet(&too_big, bad_guess);
#else
        VecSet(bad_guess, too_big);
#endif
        TS_ASSERT_THROWS_CONTAINS(solution_vector = ls.Solve(bad_guess),
                "DIVERGED_DTOL in function");
        PetscTools::Destroy(solution_vector);
        PetscTools::Destroy(good_guess);
        PetscTools::Destroy(bad_guess);
    }

    void TestAddMultipleValues()
    {
        LinearSystem syst = LinearSystem(3);

        c_matrix<double, 2, 2> small_matrix;
        c_vector<double, 2> small_vector;

        small_matrix(0,0) = 1;
        small_matrix(0,1) = 2;
        small_matrix(1,0) = 3;
        small_matrix(1,1) = 4;

        small_vector(0) = -1;
        small_vector(1) = -2;

        unsigned large_matrix_indices[2]={0,2};

        syst.AddLhsMultipleValues(large_matrix_indices, small_matrix);
        syst.AddRhsMultipleValues(large_matrix_indices, small_vector);

        syst.AssembleFinalLinearSystem();

        PetscInt lo, hi;
        syst.GetOwnershipRange(lo, hi);

        if (lo <=0 && 0<hi)
        {
            TS_ASSERT_EQUALS(syst.GetMatrixElement(0,0), 1);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(0,1), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(0,2), 2);

            TS_ASSERT_EQUALS(syst.GetRhsVectorElement(0), -1);
        }
        if (lo <=1 && 1<hi)
        {
            TS_ASSERT_EQUALS(syst.GetMatrixElement(1,0), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(1,1), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(1,2), 0);

            TS_ASSERT_EQUALS(syst.GetRhsVectorElement(1), 0);
        }
        if (lo <=2 && 2<hi)
        {
            TS_ASSERT_EQUALS(syst.GetMatrixElement(2,0), 3);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(2,1), 0);
            TS_ASSERT_EQUALS(syst.GetMatrixElement(2,2), 4);

            TS_ASSERT_EQUALS(syst.GetRhsVectorElement(2), -2);
        }
    }

    void TestSymmetricMatrix()
    {
        LinearSystem ls = LinearSystem(3);

        TS_ASSERT(!ls.IsMatrixSymmetric());
        ls.SetMatrixIsSymmetric();
        TS_ASSERT(ls.IsMatrixSymmetric());

        // Enter symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double)(abs(row-col)));
            }
        }

        // arbitrary
        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);
        ls.AssembleFinalLinearSystem();

        Vec solution_vector;
        solution_vector = ls.Solve();

        double expected_solution[3]={25.0,0.0,7.0};
        PetscInt lo, hi;
        ls.GetOwnershipRange(lo, hi);
        double* p_solution_elements_array;
        VecGetArray(solution_vector, &p_solution_elements_array);
        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], expected_solution[global_index], 1e-5);
            }
        }

        VecRestoreArray(solution_vector, &p_solution_elements_array);

        PetscTools::Destroy(solution_vector);

        // coverage
        ls.SetMatrixIsSymmetric(false);
        TS_ASSERT(!ls.IsMatrixSymmetric());
    }

    void TestNonSymmetricMatrix()
    {
        LinearSystem ls = LinearSystem(3);

        // Enter non-symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double)(10+row-col));
            }
        }

        // Arbitrary
        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);
        ls.AssembleFinalLinearSystem();

        // Solving should be fine
        Vec solution_vector;
        solution_vector = ls.Solve();

        LinearSystem ls2 = LinearSystem(3);
        ls2.SetMatrixIsSymmetric();

        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls2.SetMatrixElement(row, col, (double)(10+row-col));
            }
        }
        ls2.AssembleFinalLinearSystem();

        // What happens when we solve?
        Vec solution_vector2;
        solution_vector2 = ls2.Solve();

        // Check answers
        double expected_solution[3]={-68.0,6.0,80.0};
        PetscInt lo, hi;
        ls.GetOwnershipRange(lo, hi);
        double* p_solution_elements_array,* p_solution_elements_array2;
        VecGetArray(solution_vector, &p_solution_elements_array);
        VecGetArray(solution_vector2, &p_solution_elements_array2);
        for (int global_index=0; global_index<3; global_index++)
        {
            int local_index = global_index-lo;
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_DELTA(p_solution_elements_array[local_index], expected_solution[global_index], 1e-5);
                TS_ASSERT_LESS_THAN(2, fabs(p_solution_elements_array2[local_index] - expected_solution[global_index]));
                //Diverges from expected by more than 2
            }
        }

        VecRestoreArray(solution_vector, &p_solution_elements_array);
        VecRestoreArray(solution_vector2, &p_solution_elements_array2);
        PetscTools::Destroy(solution_vector2);
        PetscTools::Destroy(solution_vector);
    }

    void TestGetSetKSP()
    {
        {
            // Test that small systems don't have preconditioning
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 3) //PETSc 3.0 to PETSc 3.3
        //The PETSc developers changed this one, but later changed it back again!
            const PCType pc_small;
#else
            PCType pc_small;
#endif
            LinearSystem ls_small = LinearSystem(6);
            ls_small.SetPcType("jacobi"); //Will be over-ridden because the system is small
            ls_small.AssembleFinalLinearSystem();
            Vec solution_vector_small;
            solution_vector_small = ls_small.Solve();
            PetscTools::Destroy(solution_vector_small);
            PC prec_small;
            KSPGetPC(ls_small.mKspSolver, &prec_small);
            PCGetType(prec_small, &pc_small);
            TS_ASSERT( strcmp(pc_small,"none")==0 );
        }
        // Set relative tolerance before first solve
        LinearSystem ls = LinearSystem(7);
        ls.SetRelativeTolerance(1e-3);
        ls.SetKspType("cg");
        ls.SetPcType("jacobi"); //Will not be over-ridden because the system is larger...
        ls.AssembleFinalLinearSystem();
        Vec solution_vector;
        solution_vector = ls.Solve();
        PetscTools::Destroy(solution_vector);
        PetscReal rtol, atol, dtol;
        int maxits;
        KSPGetTolerances(ls.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, 1e-3);

        // Others should be their PETSc defaults (unless we've done different)
        TS_ASSERT_EQUALS(atol, 1e-50);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 1000); /// \todo #1695 Test against member variable in LinearSystem. At the moment this value depends on whether any previous test called ResetKspSolver() (maxits=1000) or not (maxits=10000).

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 3) //PETSc 3.0 to PETSc 3.3
        //The PETSc developers changed this one, but later changed it back again!
        const KSPType solver;
        const PCType pc;
#else
        KSPType solver;
        PCType pc;
#endif

        PC prec;
        KSPGetType(ls.mKspSolver, &solver);
        KSPGetPC(ls.mKspSolver, &prec);
        PCGetType(prec, &pc);
        TS_ASSERT( strcmp(solver,"cg")==0 );
        TS_ASSERT( strcmp(pc,"jacobi")==0 );
        ls.SetKspType("gmres");
        ls.SetPcType("ilu");

        // Test that we can change the solver type after its first use
        KSPGetType(ls.mKspSolver, &solver);
        PCGetType(prec, &pc);
        TS_ASSERT( strcmp(solver,"gmres")==0 );
        TS_ASSERT( strcmp(pc,"ilu")==0 );

        /////////////////////////////////
        // Set relative tolerance after first solve
        /////////////////////////////////
        ls.SetRelativeTolerance(1e-4);
        KSPGetTolerances(ls.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, 1e-4);
        TS_ASSERT_EQUALS(atol, 1e-50);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 1000); /// \todo #1695 Test against member variable in LinearSystem. At the moment this value depends on whether any previous test called ResetKspSolver() (maxits=1000) or not (maxits=10000).

        /////////////////////////////////
        // Set abs tolerance before first solve
        //////////////////////////////////
        LinearSystem ls2 = LinearSystem(5);
        ls2.SetAbsoluteTolerance(1e-3);
        ls2.AssembleFinalLinearSystem();
        Vec solution_vector2;
        solution_vector2 = ls2.Solve();
        PetscTools::Destroy(solution_vector2);
        KSPGetTolerances(ls2.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, DBL_EPSILON);
        TS_ASSERT_EQUALS(atol, 1e-3);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 1000); /// \todo #1695 Test against member variable in LinearSystem. At the moment this value depends on whether any previous test called ResetKspSolver() (maxits=1000) or not (maxits=10000).

        ///////////////////////////////////
        // Set abs tolerance after first solve
        ////////////////////////////////////
        ls2.SetAbsoluteTolerance(1e-2);
        Vec solution_vector3;
        solution_vector3 = ls2.Solve();
        PetscTools::Destroy(solution_vector3);
        KSPGetTolerances(ls2.mKspSolver, &rtol, &atol, &dtol, &maxits);
        TS_ASSERT_EQUALS(rtol, DBL_EPSILON);
        TS_ASSERT_EQUALS(atol, 1e-2);
        TS_ASSERT_EQUALS(dtol, 10000.0);
        TS_ASSERT_EQUALS(maxits, 1000); /// \todo #1695 Test against member variable in LinearSystem

    }

    void TestPetscSaveAndLoad()
    {
        // Archive
        OutputFileHandler handler("archive", false);
        std::string archive_filename_lhs, archive_filename_rhs;
        archive_filename_lhs = handler.GetOutputDirectoryFullPath() + "direct_lhs.mat";
        archive_filename_rhs = handler.GetOutputDirectoryFullPath() + "direct_rhs.vec";

        // Make a linear system
        LinearSystem ls = LinearSystem(3);
        ls.SetMatrixIsSymmetric();

        // Enter symmetric data
        for (int row=0; row<3; row++)
        {
            for (int col=0; col<3; col++)
            {
                ls.SetMatrixElement(row, col, (double)(abs(row-col)));
            }
        }

        // Arbitrary
        ls.SetRhsVectorElement(0, 14.0);
        ls.SetRhsVectorElement(1, 32.0);
        ls.SetRhsVectorElement(2, 50.0);

        ls.AssembleFinalLinearSystem();

        // SAVE
        {
            PetscViewer vec_viewer;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
            PetscViewerFileType type = PETSC_FILE_CREATE;
#else
            PetscFileMode type = FILE_MODE_WRITE;
#endif
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_rhs.c_str(), type, &vec_viewer);

            VecView(ls.GetRhsVector(), vec_viewer);
            PetscViewerDestroy(PETSC_DESTROY_PARAM(vec_viewer));

            PetscViewer mat_viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_lhs.c_str(), type, &mat_viewer);

            MatView(ls.GetLhsMatrix(), mat_viewer);
            PetscViewerDestroy(PETSC_DESTROY_PARAM(mat_viewer));
        }
        // LOAD
        {

            PetscViewer vec_viewer;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
            PetscViewerFileType type = PETSC_FILE_RDONLY;
#else
            PetscFileMode type = FILE_MODE_READ;
#endif
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_rhs.c_str(), type, &vec_viewer);
            Vec new_vec;

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
            VecCreate(PETSC_COMM_WORLD, &new_vec);
            VecSetType(new_vec, PETSC_NULL);
            VecLoad(new_vec, vec_viewer);
#else
            VecLoad(vec_viewer, PETSC_NULL, &new_vec);
#endif
            PetscViewerDestroy(PETSC_DESTROY_PARAM(vec_viewer));

            int lo, hi;
            VecGetOwnershipRange(new_vec, &lo, &hi);
            std::vector<double> answer;
            answer.push_back(14.0);
            answer.push_back(32.0);
            answer.push_back(50.0);

            double* p_vec_values;
            VecGetArray(new_vec, &p_vec_values);

            for ( int i = lo; i < hi; i++ )
            {
                TS_ASSERT_DELTA(p_vec_values[i-lo], answer[i], 1e-9);
            }

            PetscTools::Destroy(new_vec);

            PetscViewer mat_viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_lhs.c_str(), type, &mat_viewer);
            Mat new_mat;

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
            MatCreate(PETSC_COMM_WORLD, &new_mat);
            MatSetType(new_mat, PETSC_NULL);
            MatLoad(new_mat, mat_viewer);
#else
            MatLoad(mat_viewer, PETSC_NULL, &new_mat);
#endif

            PetscViewerDestroy(PETSC_DESTROY_PARAM(mat_viewer));

            for (int row=lo; row<hi; row++)
            {
                // Get a whole row out of the matrix and check it
                PetscInt row_as_array[1];
                row_as_array[0] = row;
                PetscInt col_as_array[3];
                for (int col=0; col<3; col++)
                {
                   col_as_array[col] = col;
                }
                double ret_array[3];
                MatGetValues(new_mat, 1, row_as_array, 3, col_as_array, ret_array);

                for (int col=0; col<3; col++)
                {
                    TS_ASSERT_DELTA(ret_array[col], (double)(abs(row-col)), 1e-9);
                }
            }

            PetscTools::Destroy(new_mat);
        }
    }

    void TestSaveAndLoadLinearSystem()
    {
        // Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "linear_system.arch";

        int lo, hi;
        unsigned size = 5;
        std::vector<double> rhs_values;
        rhs_values.push_back(14.0);
        rhs_values.push_back(32.0);
        rhs_values.push_back(50.0);
        rhs_values.push_back(50.0);
        rhs_values.push_back(50.0);
        TS_ASSERT_EQUALS(rhs_values.size(), size);
        // SAVE
        {
            LinearSystem ls = LinearSystem(size);
            Mat temp_mat=ls.GetLhsMatrix();
            PetscBool symm_set, is_symmetric;
            is_symmetric = PETSC_FALSE;
            MatIsSymmetricKnown(temp_mat, &symm_set, &is_symmetric);
            TS_ASSERT_EQUALS(symm_set, PETSC_FALSE);
            TS_ASSERT_EQUALS(is_symmetric, PETSC_FALSE);

            ls.SetMatrixIsSymmetric();

            MatIsSymmetricKnown(temp_mat, &symm_set, &is_symmetric);
            TS_ASSERT_EQUALS(symm_set, PETSC_TRUE);
            TS_ASSERT_EQUALS(is_symmetric, PETSC_TRUE);

            // Enter symmetric data
            for (unsigned row=0; row<size; row++)
            {
                for (unsigned col=0; col<size; col++)
                {
                    ls.SetMatrixElement(row, col, (row == col)?(row+1.0):0.0);
                }
            }
            ls.AssembleFinalLinearSystem();

            for (unsigned i=0; i<size; i++)
            {
                ls.SetRhsVectorElement(i, rhs_values[i]);
            }
            ls.SetKspType("cg");
            ls.SetPcType("none");

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            LinearSystem* const p_linear_system = &ls;
            output_arch << p_linear_system;

            TS_ASSERT_EQUALS(p_linear_system->GetSize(), size);
            VecGetOwnershipRange(p_linear_system->GetRhsVector(), &lo, &hi);

            for ( int i = lo; i < hi; i++ )
            {
                TS_ASSERT_DELTA(p_linear_system->GetRhsVectorElement(i), rhs_values[i], 1e-9);
            }

            for (unsigned row=(unsigned)lo; row<(unsigned)hi; row++)
            {
                for (unsigned col=0; col<size; col++)
                {
                    TS_ASSERT_DELTA(p_linear_system->GetMatrixElement(row, col), (row == col)?(row+1.0):0.0, 1e-9);
                }
            }

        }
        // LOAD
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // LinearSystem linear_system(3);
            LinearSystem* p_linear_system;//=&linear_system;
            input_arch >> p_linear_system;

            // Check that structural symmetry is preserved
            PetscBool symm_set, is_symmetric;
            is_symmetric=PETSC_FALSE;
            MatIsSymmetricKnown(p_linear_system->GetLhsMatrix(), &symm_set, &is_symmetric);
            TS_ASSERT_EQUALS(symm_set, PETSC_TRUE);
            TS_ASSERT_EQUALS(is_symmetric, PETSC_TRUE);


            TS_ASSERT_EQUALS(p_linear_system->GetSize(), size);

            int saved_lo, saved_hi;
            VecGetOwnershipRange(p_linear_system->GetRhsVector(), &saved_lo, &saved_hi);

            TS_ASSERT_EQUALS(hi, saved_hi);
            TS_ASSERT_EQUALS(lo, saved_lo);

            for ( int i = lo; i < hi; i++ )
            {
                TS_ASSERT_DELTA(p_linear_system->GetRhsVectorElement(i), rhs_values[i], 1e-9);
            }

            for (unsigned row=(unsigned)lo; row<(unsigned)hi; row++)
            {
                for (unsigned col=0; col<size; col++)
                {
                    TS_ASSERT_DELTA(p_linear_system->GetMatrixElement(row, col), (row == col)?(row+1.0):0.0, 1e-9);
                }
            }

            // Check archiving of KSP/PC types
            Vec solution_vector3;
            solution_vector3 = p_linear_system->Solve();
            PetscTools::Destroy(solution_vector3);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 3) //PETSc 3.0 to PETSc 3.3
            //The PETSc developers changed this one, but later changed it back again!
            const KSPType solver;
            const PCType pc;
#else
            KSPType solver;
            PCType pc;
#endif
            PC prec;
            KSPGetType(p_linear_system->mKspSolver, &solver);
            KSPGetPC(p_linear_system->mKspSolver, &prec);
            PCGetType(prec, &pc);

            TS_ASSERT(strcmp(solver, "cg") == 0);
            TS_ASSERT(strcmp(pc, "none") == 0);
            delete p_linear_system;
        }
    }
    // This test causes test-suite to exit abnormally with PETSc 3.2
    // and with MPICH-1.
    // Note: it's on the LDU solve.  Investigate Hypre.
    void TestConsecutiveSolvesDifferentPreconditioner()
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        Mat system_matrix;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        Vec system_rhs;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

        PetscTools::Destroy(parallel_layout);

        LinearSystem ls = LinearSystem(system_rhs, system_matrix);

        ls.SetAbsoluteTolerance(1e-9);
        ls.SetKspType("cg");

        ls.SetPcType("bjacobi");
        Vec solution = ls.Solve(/*no guess provided*/);
        unsigned block_jacobi_its = ls.GetNumIterations();
        PetscTools::Destroy(solution);

        ls.SetPcType("ldufactorisation");
        solution = ls.Solve(/*no guess provided*/);
        unsigned ldu_its = ls.GetNumIterations();
        PetscTools::Destroy(solution);

        ls.SetPcType("bjacobi");
        solution = ls.Solve(/*no guess provided*/);
        unsigned second_block_jacobi_its = ls.GetNumIterations();
        PetscTools::Destroy(solution);

        TS_ASSERT_DIFFERS(block_jacobi_its, ldu_its)
        TS_ASSERT_DIFFERS(ldu_its, second_block_jacobi_its)

        PetscTools::Destroy(system_matrix);
        PetscTools::Destroy(system_rhs);
    }

//    void TestSingularSolves()
//    {
//        LinearSystem ls(2);
//
//        ls.SetMatrixElement(0, 0, 2);
//        ls.SetMatrixElement(0, 1, 2);
//        ls.SetMatrixElement(1, 0, 2);
//        ls.SetMatrixElement(1, 1, 2);
//
//        ls.SetRhsVectorElement(0, 100.0);
//        ls.SetRhsVectorElement(1, 100.0);
//
//        ls.AssembleFinalLinearSystem();
//
//        Vec x = ls.Solve();
//        ReplicatableVector xx(x);
//
//        std::cout << xx[0] << " " << xx[1] << "\n"; //solves fine without null space?
//    }

    void TestDifferentMatrixForPreconditioning()
    {

        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        Mat system_matrix;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        Vec system_rhs;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

        PetscTools::Destroy(parallel_layout);

        unsigned num_it_same_mat=0, num_it_diff_mat=1;

        {
            LinearSystem ls(system_rhs, system_matrix);
            ls.SetKspType("cg");

            Vec solution = ls.Solve();
            num_it_same_mat = ls.GetNumIterations();

            PetscTools::Destroy(solution);
        }

        /*
         * Basic test, we pretend matrix for preconditioning assembly is different
         * from LHS but then we set LHS as preconditioning matrix. Number of iterations
         * should agree.
         */
        {
            LinearSystem ls_diff_precond(system_rhs, system_matrix);
            ls_diff_precond.SetKspType("cg");

            // For coverage
            TS_ASSERT_THROWS_THIS(ls_diff_precond.rGetPrecondMatrix(), "LHS matrix used for preconditioner construction");

            ls_diff_precond.SetPrecondMatrixIsDifferentFromLhs();
            MatCopy(ls_diff_precond.GetLhsMatrix(), ls_diff_precond.rGetPrecondMatrix(), DIFFERENT_NONZERO_PATTERN);

            Vec solution = ls_diff_precond.Solve();
            num_it_diff_mat = ls_diff_precond.GetNumIterations();

            PetscTools::Destroy(solution);
        }

        TS_ASSERT_EQUALS(num_it_diff_mat, num_it_same_mat);

        /*
         * Setting the identity matrix as a preconditioner is equivalent to no preconditioning
         */
        {
            LinearSystem ls_diff_precond(system_rhs, system_matrix);
            ls_diff_precond.SetKspType("cg");

            ls_diff_precond.SetPrecondMatrixIsDifferentFromLhs();
            Mat& r_identity_matrix=ls_diff_precond.rGetPrecondMatrix();

            for (unsigned row_col=0; row_col<num_nodes; row_col++)
            {
                PetscMatTools::SetElement(r_identity_matrix, row_col, row_col, 1.0);
            }
            ls_diff_precond.FinalisePrecondMatrix();

            Vec solution = ls_diff_precond.Solve();
            num_it_diff_mat = ls_diff_precond.GetNumIterations();

            PetscTools::Destroy(solution);
        }

        TS_ASSERT_EQUALS(num_it_diff_mat, 80u); // It takes 80 iterations if you run with ls.SetPcType("none");
        TS_ASSERT_LESS_THAN(num_it_same_mat, num_it_diff_mat);

        PetscTools::Destroy(system_rhs);
        PetscTools::Destroy(system_matrix);
    }

    void TestFixedNumberOfIterations()
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        Mat system_matrix;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        Vec system_rhs;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

        LinearSystem ls = LinearSystem(system_rhs, system_matrix);

        ls.SetMatrixIsSymmetric();
        ls.SetKspType("cg");
        ls.SetPcType("jacobi");
        ls.SetAbsoluteTolerance(1e-4);

        Vec guess;
        VecDuplicate(parallel_layout, &guess);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2)
        PetscScalar zero = 0.0;
        VecSet(&zero, guess);
#else
        VecSet(guess, 0.0);
#endif

        /*
         * Use fixed number of iterations, updating the number of iterations to perform every other solve.
         */
        TS_ASSERT_THROWS_NOTHING(ls.SetUseFixedNumberIterations(true, 2));

        Vec solution = ls.Solve(guess);

        unsigned chebyshev_its = ls.GetNumIterations();
        TS_ASSERT_EQUALS(chebyshev_its, 52u);

        Vec new_solution;
        Vec difference;
        VecDuplicate(parallel_layout, &difference);
        PetscReal l_inf_norm;

        /*
         * Solve using previous solution as new guess. If we were checking convergence
         * normally it would take 0 iterations to solve. Since we set fixed number of
         * iterations based on first solve, it will take the same number as above.
         */
        new_solution = ls.Solve(solution);
        chebyshev_its = ls.GetNumIterations();

        TS_ASSERT_EQUALS(chebyshev_its, 52u);

        PetscVecTools::WAXPY(difference, -1.0, new_solution, solution);
        VecNorm(difference, NORM_INFINITY, &l_inf_norm);
        TS_ASSERT_DELTA(l_inf_norm, 0.0, 4e-4);
        PetscTools::Destroy(new_solution);

        /*
         * Solve using previous solution as new guess takes 0 iterations as
         * it is performed with tolerance-based stop criteria.
         */
        new_solution = ls.Solve(solution);
        chebyshev_its = ls.GetNumIterations();

        TS_ASSERT_EQUALS(chebyshev_its, 0u);

        PetscVecTools::WAXPY(difference, -1.0, new_solution, solution);
        VecNorm(difference, NORM_INFINITY, &l_inf_norm);
        TS_ASSERT_DELTA(l_inf_norm, 0.0, 4e-4);
        PetscTools::Destroy(new_solution);

        /*
         * Solve with initial guess should take 52 iterations but the solver
         * uses the same number of iterations as the previous call.
         */
        new_solution = ls.Solve(guess);
        chebyshev_its = ls.GetNumIterations();

#if ((PETSC_VERSION_MAJOR==3) || (PETSC_VERSION_MAJOR==2 && PETSC_VERSION_MINOR==3 && PETSC_VERSION_SUBMINOR==3))
        TS_ASSERT_EQUALS(chebyshev_its, 0u);
#else
        TS_ASSERT_EQUALS(chebyshev_its, 1u);
#endif

        PetscVecTools::WAXPY(difference, -1.0, new_solution, solution);
        VecNorm(difference, NORM_INFINITY, &l_inf_norm);

#if ((PETSC_VERSION_MAJOR==3) || (PETSC_VERSION_MAJOR==2 && PETSC_VERSION_MINOR==3 && PETSC_VERSION_SUBMINOR==3))
        TS_ASSERT_DELTA(l_inf_norm, 60.47, 2.0);
#else
        TS_ASSERT_DELTA(l_inf_norm, 22.43, 2.0);
#endif

        PetscTools::Destroy(new_solution);

        /*
         * Solve with initial guess and tolerance-based stop criteria takes 52 iterations
         */
        new_solution = ls.Solve(guess);
        chebyshev_its = ls.GetNumIterations();

        TS_ASSERT_EQUALS(chebyshev_its, 52u);

        PetscVecTools::WAXPY(difference, -1.0, new_solution, solution);
        VecNorm(difference, NORM_INFINITY, &l_inf_norm);
        TS_ASSERT_DELTA(l_inf_norm, 0.0, 4e-4);

        PetscTools::Destroy(solution);
        PetscTools::Destroy(new_solution);
        PetscTools::Destroy(difference);

        PetscTools::Destroy(system_matrix);
        PetscTools::Destroy(system_rhs);
        PetscTools::Destroy(parallel_layout);
        PetscTools::Destroy(guess);
    }

    void TestFixedNumberOfIterationsRelativeToleranceCoverage()
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        Mat system_matrix;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        Vec system_rhs;
        // Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

        LinearSystem ls = LinearSystem(system_rhs, system_matrix);

        ls.SetMatrixIsSymmetric();
        ls.SetKspType("cg");
        ls.SetPcType("jacobi");
        ls.SetRelativeTolerance(1e-6); // Covering the use of relative tolerance for solution

        Vec guess;
        VecDuplicate(parallel_layout, &guess);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2)
        PetscScalar zero = 0.0;
        VecSet(&zero, guess);
#else
        VecSet(guess, 0.0);
#endif

        /*
         *  Use fixed number of iterations, updating the number of iterations to perform every other solve.
         */
        TS_ASSERT_THROWS_NOTHING(ls.SetUseFixedNumberIterations(true, 2));

        Vec solution = ls.Solve(guess);

        unsigned chebyshev_its = ls.GetNumIterations();
        TS_ASSERT_EQUALS(chebyshev_its, 40u);

        Vec new_solution;
        Vec difference;
        VecDuplicate(parallel_layout, &difference);
        PetscReal l_inf_norm;

        /*
         * Solve using previous solution as new guess. If we were checking convergence
         * normally it would take 0 iterations to solve. Since we set fixed number of
         * iterations based on first solve, it will take the same number as above.
         */
        new_solution = ls.Solve(solution);
        chebyshev_its = ls.GetNumIterations();

        TS_ASSERT_EQUALS(chebyshev_its, 40u);

        PetscVecTools::WAXPY(difference, -1.0, new_solution, solution);
        VecNorm(difference, NORM_INFINITY, &l_inf_norm);
        TS_ASSERT_DELTA(l_inf_norm, 0.0, 1e-3);
        PetscTools::Destroy(new_solution);

        PetscTools::Destroy(solution);
        PetscTools::Destroy(difference);

        PetscTools::Destroy(system_matrix);
        PetscTools::Destroy(system_rhs);
        PetscTools::Destroy(parallel_layout);
        PetscTools::Destroy(guess);
    }

    void TestSolveZerosInitialGuessForSmallRhs()
    {
        LinearSystem ls(2);

        ls.SetMatrixElement(0, 0, 1.0);
        ls.SetMatrixElement(0, 1, 0.0);
        ls.SetMatrixElement(1, 0, 0.0);
        ls.SetMatrixElement(1, 1, 1.0);

        ls.SetRhsVectorElement(0, 1e-17);
        ls.SetRhsVectorElement(1, 0.0);

        Vec init_cond = PetscTools::CreateAndSetVec(2, 0.0);
        PetscVecTools::SetElement(init_cond, 0, 0.01);

        ls.AssembleFinalLinearSystem();

        Vec solution_vector;
        solution_vector = ls.Solve(init_cond);

        ReplicatableVector solution_vector_repl(solution_vector);
        TS_ASSERT_DELTA(solution_vector_repl[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(solution_vector_repl[1], 0.0, 1e-6);

        PetscTools::Destroy(init_cond);
        PetscTools::Destroy(solution_vector);
    }

    /** See #1834 */
    void xxxxTestSolveShowingPetscRubbishInEvenMoreWays()
    {
        LinearSystem ls(2);

        ls.SetMatrixElement(0, 0, 1.0);
        ls.SetMatrixElement(0, 1, 0.0);
        ls.SetMatrixElement(1, 0, 0.0);
        ls.SetMatrixElement(1, 1, 1.0);

        ls.SetRhsVectorElement(0, log(0.0)); // -inf
        ls.SetRhsVectorElement(1, 0.0);

        ls.AssembleFinalLinearSystem();

        ls.DisplayRhs();

        Vec solution_vector;
        solution_vector = ls.Solve(); // SHOULD fail with an error, but doesn't
        PetscTools::Destroy(solution_vector);
    }

    // This test should be the last in the suite
    void TestSetFromOptions()
    {
        LinearSystem ls = LinearSystem(5);
        PetscTools::SetOption("-ksp_type", "gmres");
        PetscTools::SetOption("-pc_type", "jacobi");
        ls.AssembleFinalLinearSystem();

        ls.SetKspType("cg"); //Not really -- see above
        ls.SetPcType("ilu"); //Not really -- see above
        Vec solution_vector3;
        solution_vector3 = ls.Solve();
        PetscTools::Destroy(solution_vector3);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 3) //PETSc 3.0 to PETSc 3.3
        //The PETSc developers changed this one, but later changed it back again!
        const KSPType solver;
        const PCType pc;
#else
        KSPType solver;
        PCType pc;
#endif
        PC prec;
        KSPGetType(ls.mKspSolver, &solver);
        KSPGetPC(ls.mKspSolver, &prec);
        PCGetType(prec, &pc);

        TS_ASSERT( strcmp(solver,"gmres")==0 );
        TS_ASSERT( strcmp(pc,"jacobi")==0 );
    }
    // The above test should be last in the suite

};
#endif //_TESTLINEARSYSTEM_HPP_
