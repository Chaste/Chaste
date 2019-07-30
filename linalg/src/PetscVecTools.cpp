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

#include "PetscVecTools.hpp"
#include "PetscTools.hpp"
#include <petscviewer.h>
#include <cassert>
#include "DistributedVectorFactory.hpp"
#include "DistributedVector.hpp"
#include "PetscException.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

void PetscVecTools::Finalise(Vec vector)
{
    VecAssemblyBegin(vector);
    VecAssemblyEnd(vector);
}

void PetscVecTools::SetElement(Vec vector, PetscInt row, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(vector, lo, hi);

    if (row >= lo && row < hi)
    {
        VecSetValues(vector, 1, &row, &value, INSERT_VALUES);
    }
}

void PetscVecTools::AddToElement(Vec vector, PetscInt row, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(vector, lo, hi);

    if (row >= lo && row < hi)
    {
        VecSetValues(vector, 1, &row, &value, ADD_VALUES);
    }
}

void PetscVecTools::Display(Vec vector)
{
    //Give full precision, scientific notation
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7) // PETSc 3.7+
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
#else
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
#endif
    VecView(vector, PETSC_VIEWER_STDOUT_WORLD);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7) // PETSc 3.7+
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
#endif
}

void PetscVecTools::Zero(Vec vector)
{
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PetscScalar zero = 0.0;
    VecSet(&zero, vector);
#else
    VecZeroEntries(vector);
#endif
}

unsigned PetscVecTools::GetSize(Vec vector)
{
    PetscInt size;
    VecGetSize(vector, &size);
    return (unsigned) size;
}

void PetscVecTools::GetOwnershipRange(Vec vector, PetscInt& lo, PetscInt& hi)
{
    VecGetOwnershipRange(vector, &lo, &hi);
}

double PetscVecTools::GetElement(Vec vector, PetscInt row)
{
    PetscInt lo, hi;
    GetOwnershipRange(vector, lo, hi);
    assert(lo <= row && row < hi);

    double* p_vector;
    PetscInt local_index = row-lo;
    VecGetArray(vector, &p_vector);
    double answer = p_vector[local_index];
    VecRestoreArray(vector, &p_vector);

    return answer;
}

void PetscVecTools::AddScaledVector(Vec y, Vec x, double scaleFactor)
{
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    VecAXPY(&scaleFactor, x, y);
#else
    VecAXPY(y, scaleFactor, x);
#endif
}

void PetscVecTools::Scale(Vec vector, double scaleFactor)
{
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PETSCEXCEPT( VecScale(&scaleFactor, vector) );
#else
    PETSCEXCEPT( VecScale(vector, scaleFactor) );
#endif
}

void PetscVecTools::WAXPY(Vec w, double a, Vec x, Vec y)
{
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PETSCEXCEPT( VecWAXPY(&a, x, y, w) );
#else
    PETSCEXCEPT( VecWAXPY(w, a, x, y) );
#endif
}

void PetscVecTools::SetupInterleavedVectorScatterGather(Vec interleavedVec, VecScatter& rFirstVariableScatterContext, VecScatter& rSecondVariableScatterContext)
{
    PetscInt num_rows, num_local_rows;

    VecGetSize(interleavedVec, &num_rows);
    VecGetLocalSize(interleavedVec, &num_local_rows);

    IS A11_rows, A22_rows;
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 2, &A11_rows);
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 1, 2, &A22_rows);

    IS all_vector;
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 1, &all_vector);

    unsigned subvector_num_rows = num_rows/2;
    unsigned subvector_local_rows = num_local_rows/2;
    Vec x1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    Vec x2_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);

    VecScatterCreate(interleavedVec, A11_rows, x1_subvector, all_vector, &rFirstVariableScatterContext);
    VecScatterCreate(interleavedVec, A22_rows, x2_subvector, all_vector, &rSecondVariableScatterContext);

    PetscTools::Destroy(x1_subvector);
    PetscTools::Destroy(x2_subvector);

    ISDestroy(PETSC_DESTROY_PARAM(A11_rows));
    ISDestroy(PETSC_DESTROY_PARAM(A22_rows));
    ISDestroy(PETSC_DESTROY_PARAM(all_vector));
}

void PetscVecTools::DoInterleavedVecScatter(Vec interleavedVec, VecScatter firstVariableScatterContext, Vec firstVariableVec, VecScatter secondVariableScatterContefirstVariableScatterContextt, Vec secondVariableVec)
{
    PetscScalar *p_interleaved_vec;
    PetscScalar *p_1st_variable_vec;
    PetscScalar *p_2nd_variable_vec;

    //  // Note that the first vector is normally locked read-only when this is used
    //    PetscInt lock;
    //    VecLockGet(interleavedVec, &lock);
    //    assert(lock == 1);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        // Request read-only access properly
        VecGetArrayRead(interleavedVec, (const PetscScalar**)&p_interleaved_vec);
#else
        VecGetArray(interleavedVec, &p_interleaved_vec);
#endif
    VecGetArray(firstVariableVec, &p_1st_variable_vec);
    VecGetArray(secondVariableVec, &p_2nd_variable_vec);

    PetscInt vec_local_size;
    VecGetLocalSize(interleavedVec, &vec_local_size);
    assert(vec_local_size%2 == 0);

    for (PetscInt local_index=0; local_index<vec_local_size/2; local_index++)
    {
        p_1st_variable_vec[local_index] = p_interleaved_vec[2*local_index];
        p_2nd_variable_vec[local_index] = p_interleaved_vec[2*local_index+1];
    }

    VecRestoreArray(interleavedVec, &p_interleaved_vec);
    VecRestoreArray(firstVariableVec, &p_1st_variable_vec);
    VecRestoreArray(secondVariableVec, &p_2nd_variable_vec);

//    DistributedVectorFactory factory(vec_size/2);
//
//    DistributedVector dist_inter_vec = factory.CreateDistributedVector(interleavedVec);
//    DistributedVector::Stripe dist_inter_vec_1st_var(dist_inter_vec, 0);
//    DistributedVector::Stripe dist_inter_vec_2nd_var(dist_inter_vec, 1);
//
//    DistributedVector dist_1st_var_vec = factory.CreateDistributedVector(firstVariableVec);
//    DistributedVector dist_2nd_var_vec = factory.CreateDistributedVector(secondVariableVec);
//
//    for (DistributedVector::Iterator index = dist_1st_var_vec.Begin();
//         index!= dist_1st_var_vec.End();
//         ++index)
//    {
//        dist_1st_var_vec[index] = dist_inter_vec_1st_var[index];
//        dist_2nd_var_vec[index] = dist_inter_vec_2nd_var[index];
//    }


//    //PETSc-3.x.x or PETSc-2.3.3
//#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
//    VecScatterBegin(firstVariableScatterContext, interleavedVec, firstVariableVec, INSERT_VALUES, SCATTER_FORWARD);
//    VecScatterEnd(firstVariableScatterContext, interleavedVec, firstVariableVec, INSERT_VALUES, SCATTER_FORWARD);
//#else
//    VecScatterBegin(interleavedVec, firstVariableVec, INSERT_VALUES, SCATTER_FORWARD, firstVariableScatterContext);
//    VecScatterEnd(interleavedVec, firstVariableVec, INSERT_VALUES, SCATTER_FORWARD, firstVariableScatterContext);
//#endif
//
//    //PETSc-3.x.x or PETSc-2.3.3
//#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
//    VecScatterBegin(secondVariableScatterContefirstVariableScatterContextt, interleavedVec, secondVariableVec, INSERT_VALUES, SCATTER_FORWARD);
//    VecScatterEnd(secondVariableScatterContefirstVariableScatterContextt, interleavedVec, secondVariableVec, INSERT_VALUES, SCATTER_FORWARD);
//#else
//    VecScatterBegin(interleavedVec, secondVariableVec, INSERT_VALUES, SCATTER_FORWARD, secondVariableScatterContefirstVariableScatterContextt);
//    VecScatterEnd(interleavedVec, secondVariableVec, INSERT_VALUES, SCATTER_FORWARD, secondVariableScatterContefirstVariableScatterContextt);
//#endif
}

void PetscVecTools::DoInterleavedVecGather(Vec interleavedVec, VecScatter firstVariableScatterContext, Vec firstVariableVec, VecScatter secondVariableScatterContext, Vec secondVariableVec)
{
    PetscScalar *p_interleaved_vec;
    PetscScalar *p_1st_variable_vec;
    PetscScalar *p_2nd_variable_vec;

    VecGetArray(interleavedVec, &p_interleaved_vec);
    VecGetArray(firstVariableVec, &p_1st_variable_vec);
    VecGetArray(secondVariableVec, &p_2nd_variable_vec);

    PetscInt vec_local_size;
    VecGetLocalSize(interleavedVec, &vec_local_size);
    assert(vec_local_size%2 == 0);

    for (PetscInt local_index=0; local_index<vec_local_size/2; local_index++)
    {
        p_interleaved_vec[2*local_index] = p_1st_variable_vec[local_index];
        p_interleaved_vec[2*local_index+1] = p_2nd_variable_vec[local_index];
    }

    VecRestoreArray(interleavedVec, &p_interleaved_vec);
    VecRestoreArray(firstVariableVec, &p_1st_variable_vec);
    VecRestoreArray(secondVariableVec, &p_2nd_variable_vec);

//#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
//    VecScatterBegin(firstVariableScatterContext, firstVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE);
//    VecScatterEnd(firstVariableScatterContext, firstVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE);
//#else
//    VecScatterBegin(firstVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE, firstVariableScatterContext);
//    VecScatterEnd(firstVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE, firstVariableScatterContext);
//#endif
//
////PETSc-3.x.x or PETSc-2.3.3
//#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
//    VecScatterBegin(secondVariableScatterContext, secondVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE);
//    VecScatterEnd(secondVariableScatterContext, secondVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE);
//#else
//    VecScatterBegin(secondVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE, secondVariableScatterContext);
//    VecScatterEnd(secondVariableVec, interleavedVec, INSERT_VALUES, SCATTER_REVERSE, secondVariableScatterContext);
//#endif
}
