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
        PetscVecTools::DoInterleavedVecScatter(interleaved_vec, first_variable_context, first_variable_vec, second_variable_context, second_variable_vec);

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
