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

#ifndef TESTPETSCVECTOOLS_HPP_
#define TESTPETSCVECTOOLS_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscVecTools.hpp" // Includes Ublas so must come before PETSc

#include <petscvec.h>

#include "DistributedVectorFactory.hpp"
#include "ReplicatableVector.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedVector.hpp"

/**
 * Tests methods in the PETSc helper class PetscVecTools that are not directly related to LinearSystem
 */
class TestPetscVecTools : public CxxTest::TestSuite
{
public:

    void TestInterleavedVecScatter()
    {
        // Vectors will be twice PROBLEM_SIZE, since this is to be used in bidomain code.
        const unsigned PROBLEM_SIZE = 10;
        DistributedVectorFactory factory(PROBLEM_SIZE);

        // Source vector = [-1 1 -2 2 ... -PROBLEM_SIZE/2+1 PROBLEM_SIZE/2+1]
        Vec interleaved_vec=factory.CreateVec(2);
        DistributedVector interleaved_dist_vec = factory.CreateDistributedVector(interleaved_vec);
        DistributedVector::Stripe first_variable(interleaved_dist_vec, 0);
        DistributedVector::Stripe second_variable(interleaved_dist_vec, 1);
        for (DistributedVector::Iterator index = interleaved_dist_vec.Begin();
             index!= interleaved_dist_vec.End();
             ++index)
        {
            first_variable[index]  = -1.0 * (index.Global+1);
            second_variable[index] =        (index.Global+1);
        }

        // Destination vectors. It is important they have a compatible parallel layout with the source vector, hence the 2nd param of CreateVec()
        PetscInt local_num_rows;
        VecGetLocalSize(interleaved_vec, &local_num_rows);
        Vec first_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);
        Vec second_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);

        // Setup and perform scatter operation
        VecScatter first_variable_context;
        VecScatter second_variable_context;
        PetscVecTools::SetupInterleavedVectorScatterGather(interleaved_vec, first_variable_context, second_variable_context);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 6) //PETSc 3.6 or later
        // When used in the preconditioner code within a KSP solve, the main (bidomain) vector will be locked.
        // Lock the vector so that modern PETSc (3.6) won't want to change it
        EXCEPT_IF_NOT(VecLockPush(interleaved_vec) == 0);
#endif

        PetscVecTools::DoInterleavedVecScatter(interleaved_vec, first_variable_context, first_variable_vec, second_variable_context, second_variable_vec);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 6) //PETSc 3.6 or later
        // Unlock the vector (for symmetry)
        EXCEPT_IF_NOT(VecLockPop(interleaved_vec) == 0);
#endif

        // Check destination vectors are [-1 -2 -3 ...] and [1 2 3 ...] respectively.
        DistributedVector dist_1st_var_vec = factory.CreateDistributedVector(first_variable_vec);
        DistributedVector dist_2nd_var_vec = factory.CreateDistributedVector(second_variable_vec);
        for (DistributedVector::Iterator index = dist_1st_var_vec.Begin();
             index!= dist_1st_var_vec.End();
             ++index)
        {
            TS_ASSERT_EQUALS(dist_1st_var_vec[index], -1.0 * (index.Global+1));
            TS_ASSERT_EQUALS(dist_2nd_var_vec[index],         index.Global+1);
        }

        PetscTools::Destroy(interleaved_vec);
        PetscTools::Destroy(first_variable_vec);
        PetscTools::Destroy(second_variable_vec);
        VecScatterDestroy(PETSC_DESTROY_PARAM(first_variable_context));
        VecScatterDestroy(PETSC_DESTROY_PARAM(second_variable_context));
    }

    void TestInterleavedVecGather()
    {
        // Vectors will be twice PROBLEM_SIZE, since this is to be used in bidomain code.
        const unsigned PROBLEM_SIZE = 10;
        DistributedVectorFactory factory(PROBLEM_SIZE);

        // Destination vector
        Vec interleaved_vec=factory.CreateVec(2);

        // Source vectors. It is important they have a compatible parallel layout with the destination vector, hence the 2nd param of CreateVec()
        PetscInt local_num_rows;
        VecGetLocalSize(interleaved_vec, &local_num_rows);
        Vec first_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);
        Vec second_variable_vec = PetscTools::CreateVec(PROBLEM_SIZE, local_num_rows/2);

        // Fill in source vectors.
        DistributedVector dist_first_variable_vec = factory.CreateDistributedVector(first_variable_vec);
        DistributedVector dist_second_variable_vec = factory.CreateDistributedVector(second_variable_vec);
        for (DistributedVector::Iterator index = dist_first_variable_vec.Begin();
             index!= dist_first_variable_vec.End();
             ++index)
        {
            dist_first_variable_vec[index]  = -1.0 * (index.Global+1);
            dist_second_variable_vec[index] =        (index.Global+1);
        }

        // Setup and perform gather operation
        VecScatter first_variable_context;
        VecScatter second_variable_context;
        PetscVecTools::SetupInterleavedVectorScatterGather(interleaved_vec, first_variable_context, second_variable_context);
        PetscVecTools::DoInterleavedVecGather(interleaved_vec, first_variable_context, first_variable_vec, second_variable_context, second_variable_vec);

        // Check that the destination vector looks like [ -1 1 -2 2 ... -PROBLEM_SIZE/2+1 PROBLEM_SIZE/2+1 ]
        DistributedVector dist_interleaved_vec = factory.CreateDistributedVector(interleaved_vec);
        DistributedVector::Stripe dist_inter_vec_1st_var(dist_interleaved_vec, 0);
        DistributedVector::Stripe dist_inter_vec_2nd_var(dist_interleaved_vec, 1);

        for (DistributedVector::Iterator index = dist_interleaved_vec.Begin();
             index!= dist_interleaved_vec.End();
             ++index)
        {
            TS_ASSERT_EQUALS(dist_inter_vec_1st_var[index], -1.0 * (index.Global+1));
            TS_ASSERT_EQUALS(dist_inter_vec_2nd_var[index],         index.Global+1 );
        }

        PetscTools::Destroy(interleaved_vec);
        PetscTools::Destroy(first_variable_vec);
        PetscTools::Destroy(second_variable_vec);
        VecScatterDestroy(PETSC_DESTROY_PARAM(first_variable_context));
        VecScatterDestroy(PETSC_DESTROY_PARAM(second_variable_context));
    }
};

#endif /* TESTPETSCVECTOOLS_HPP_ */
