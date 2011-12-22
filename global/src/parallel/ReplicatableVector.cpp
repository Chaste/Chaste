/*

Copyright (C) University of Oxford, 2005-2011

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

#include "ReplicatableVector.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

#include <cassert>
#include <iostream>

// Private methods

void ReplicatableVector::RemovePetscContext()
{
    if (mToAll != NULL)
    {
        VecScatterDestroy(mToAll);
        mToAll = NULL;
    }

    if (mReplicated != NULL)
    {
        VecDestroy(mReplicated);
        mReplicated = NULL;
    }

    if (mpData != NULL)
    {
        delete[] mpData;
        mpData = NULL;
    }
}

// Constructors & destructors

ReplicatableVector::ReplicatableVector()
    : mpData(NULL),
      mSize(0),
      mToAll(NULL),
      mReplicated(NULL)
{
}

ReplicatableVector::ReplicatableVector(Vec vec)
    : mpData(NULL),
      mSize(0),
      mToAll(NULL),
      mReplicated(NULL)
{
    ReplicatePetscVector(vec);
}

ReplicatableVector::ReplicatableVector(unsigned size)
    : mpData(NULL),
      mSize(0),
      mToAll(NULL),
      mReplicated(NULL)
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
    catch(std::bad_alloc &badAlloc)
    {
#define COVERAGE_IGNORE
        std::cout << "Failed to allocate a ReplicatableVector of size " << size  << std::endl;
        PetscTools::ReplicateException(true);
        throw badAlloc;
#undef COVERAGE_IGNORE
    }
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
    VecCreateMPIWithArray(PETSC_COMM_WORLD, hi-lo, this->GetSize(), &mpData[lo], &distributed_vec);

    // Now do the real replication
    ReplicatePetscVector(distributed_vec);

    // Clean up
    VecDestroy(distributed_vec);
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
    if (mReplicated == NULL)
    {
        // This creates mToAll (the scatter context) and mReplicated (to store values)
        VecScatterCreateToAll(vec, &mToAll, &mReplicated);
    }

    // Replicate the data
//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
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
