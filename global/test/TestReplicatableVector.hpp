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
            PetscTools::Destroy(temp_vec); // vector no longer needed

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

        PetscTools::Destroy(petsc_vec);
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

        PetscTools::Destroy(petsc_vec);
    }
};

#endif /*TESTREPLICATABLEVECTOR_HPP_*/
