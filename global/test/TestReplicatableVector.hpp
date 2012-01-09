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

#ifndef TESTREPLICATABLEVECTOR_HPP_
#define TESTREPLICATABLEVECTOR_HPP_
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "PetscTools.hpp"

static const int VEC_SIZE=10;

class TestReplicatableVector : public CxxTest::TestSuite
{
public:

    void TestBasics()
    {
        ReplicatableVector rep_vector(VEC_SIZE);
        rep_vector[0] = 15;
        rep_vector[1] = 20;

        TS_ASSERT_EQUALS(rep_vector[0], 15);
        TS_ASSERT_EQUALS(rep_vector[1], 20);
        TS_ASSERT_EQUALS(rep_vector.GetSize(), (unsigned) VEC_SIZE);

        rep_vector.Resize(5);
        TS_ASSERT_EQUALS(rep_vector.GetSize(), 5u);
    }

    void TestReplication()
    {
        for (int vec_size=0; vec_size<10; vec_size++)
        {
            int lo, hi;
            Vec temp_vec = PetscTools::CreateVec(vec_size);
            VecGetOwnershipRange(temp_vec,&lo,&hi);
            VecDestroy(temp_vec); // vector no longer needed

            ReplicatableVector rep_vector(vec_size);
            for (int global_index=0; global_index<vec_size; global_index++)
            {
                rep_vector[global_index]=lo;
            }

            rep_vector.Replicate(lo, hi);

            for (int global_index=0; global_index<vec_size; global_index++)
            {
                if (lo<=global_index && global_index<hi)
                {
                    TS_ASSERT_EQUALS(rep_vector[global_index], lo);
                }
                else
                {
                    TS_ASSERT_DIFFERS(rep_vector[global_index], lo);
                }
            }
        }
    }

    void TestPetscReplication()
    {
        int lo, hi;
        Vec petsc_vec = PetscTools::CreateVec(VEC_SIZE);
        VecGetOwnershipRange(petsc_vec,&lo,&hi);

        double* p_petsc_vec;

        VecGetArray(petsc_vec, &p_petsc_vec);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_petsc_vec[local_index]=lo;
        }
        VecRestoreArray(petsc_vec, &p_petsc_vec);
        VecAssemblyBegin(petsc_vec);
        VecAssemblyEnd(petsc_vec);

        ReplicatableVector rep_vec;
        rep_vec.ReplicatePetscVector(petsc_vec);

        for (int global_index=0; global_index<VEC_SIZE; global_index++)
        {
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_EQUALS(rep_vec[global_index], lo);
            }
            else
            {
                TS_ASSERT_DIFFERS(rep_vec[global_index], lo);
            }
        }

        VecDestroy(petsc_vec);
    }

    void TestPetscReplicationUsingAlternativeConstructor()
    {
        int lo, hi;
        Vec petsc_vec = PetscTools::CreateVec(VEC_SIZE);

        VecGetOwnershipRange(petsc_vec,&lo,&hi);

        double* p_petsc_vec;

        VecGetArray(petsc_vec, &p_petsc_vec);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_petsc_vec[local_index]=lo;
        }
        VecRestoreArray(petsc_vec, &p_petsc_vec);
        VecAssemblyBegin(petsc_vec);
        VecAssemblyEnd(petsc_vec);

        ReplicatableVector rep_vec(petsc_vec);

        for (int global_index=0; global_index<VEC_SIZE; global_index++)
        {
            if (lo<=global_index && global_index<hi)
            {
                TS_ASSERT_EQUALS(rep_vec[global_index], lo);
            }
            else
            {
                TS_ASSERT_DIFFERS(rep_vec[global_index], lo);
            }
        }

        VecDestroy(petsc_vec);
    }
};

#endif /*TESTREPLICATABLEVECTOR_HPP_*/
