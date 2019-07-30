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

#include "ReplicatableVector.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

#include <cassert>
#include <iostream>

// Private methods

void ReplicatableVector::RemovePetscContext()
{
    if (mToAll != nullptr)
    {
        VecScatterDestroy(PETSC_DESTROY_PARAM(mToAll));
        mToAll = nullptr;
    }

    if (mReplicated != nullptr)
    {
        PetscTools::Destroy(mReplicated);
        mReplicated = nullptr;
    }

    if (mpData != nullptr)
    {
        delete[] mpData;
        mpData = nullptr;
    }
}

// Constructors & destructors

ReplicatableVector::ReplicatableVector()
    : mpData(nullptr),
      mSize(0),
      mToAll(nullptr),
      mReplicated(nullptr)
{
}

ReplicatableVector::ReplicatableVector(Vec vec)
    : mpData(nullptr),
      mSize(0),
      mToAll(nullptr),
      mReplicated(nullptr)
{
    ReplicatePetscVector(vec);
}

ReplicatableVector::ReplicatableVector(unsigned size)
    : mpData(nullptr),
      mSize(0),
      mToAll(nullptr),
      mReplicated(nullptr)
{
    Resize(size);
}

ReplicatableVector::~ReplicatableVector()
{
    RemovePetscContext();
}

// Vector interface methods

unsigned ReplicatableVector::GetSize()
{
    return mSize;
}

void ReplicatableVector::Resize(unsigned size)
{
    // PETSc stuff will be out of date
    RemovePetscContext();

    mSize = size;

    try
    {
        mpData = new double[mSize];
    }
// LCOV_EXCL_START
    catch(std::bad_alloc &badAlloc)
    {
        std::cout << "Failed to allocate a ReplicatableVector of size " << size  << std::endl;
        PetscTools::ReplicateException(true);
        throw badAlloc;
    }
// LCOV_EXCL_STOP

    PetscTools::ReplicateException(false);
}

double& ReplicatableVector::operator[](unsigned index)
{
    assert(index < mSize);
    return mpData[index];
}

// The workhorse methods

void ReplicatableVector::Replicate(unsigned lo, unsigned hi)
{
    // Create a PetSC vector with the array containing the distributed data
    Vec distributed_vec;

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
    //Extra argument is block size
    VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, hi-lo, this->GetSize(), &mpData[lo], &distributed_vec);
#else
    VecCreateMPIWithArray(PETSC_COMM_WORLD, hi-lo, this->GetSize(), &mpData[lo], &distributed_vec);
#endif
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
    VecSetOption(distributed_vec, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
    VecSetOption(distributed_vec, VEC_IGNORE_OFF_PROC_ENTRIES);
#endif
    // Now do the real replication
    ReplicatePetscVector(distributed_vec);

    // Clean up
    PetscTools::Destroy(distributed_vec);
}

void ReplicatableVector::ReplicatePetscVector(Vec vec)
{
    // If the size has changed then we'll need to make a new context
    PetscInt isize;
    VecGetSize(vec, &isize);
    unsigned size = isize;

    if (this->GetSize() != size)
    {
        Resize(size);
    }
    if (mReplicated == nullptr)
    {
        // This creates mToAll (the scatter context) and mReplicated (to store values)
        VecScatterCreateToAll(vec, &mToAll, &mReplicated);
    }

    // Replicate the data
//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(mToAll, vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd  (mToAll, vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD);
#else
//PETSc-2.3.2 or previous
    VecScatterBegin(vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
    VecScatterEnd  (vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
#endif

    // Information is now in mReplicated PETSc vector
    // Copy into mData
    double* p_replicated;
    VecGetArray(mReplicated, &p_replicated);
    for (unsigned i=0; i<size; i++)
    {
        mpData[i] = p_replicated[i];
    }
}
