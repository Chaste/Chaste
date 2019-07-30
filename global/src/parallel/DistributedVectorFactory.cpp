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

#include <cassert>

#include "DistributedVectorFactory.hpp"
#include "PetscTools.hpp"

// Initialise static data
bool DistributedVectorFactory::msCheckNumberOfProcessesOnLoad = true;

void DistributedVectorFactory::CalculateOwnership(Vec vec)
{
#ifndef NDEBUG
    if (!mPetscStatusKnown)
    {
        CheckForPetsc();
    }
#endif

    // Calculate my range
    PetscInt petsc_lo, petsc_hi;
    VecGetOwnershipRange(vec, &petsc_lo, &petsc_hi);
    mGlobalLows.clear();
    mLo = (unsigned)petsc_lo;
    mHi = (unsigned)petsc_hi;
    // vector size
    PetscInt size;
    VecGetSize(vec, &size);
    mProblemSize = (unsigned) size;
    mNumProcs = PetscTools::GetNumProcs();
}

void DistributedVectorFactory::SetFromFactory(DistributedVectorFactory* pFactory)
{
    if (pFactory->GetNumProcs() != mNumProcs)
    {
        EXCEPTION("Cannot set from a factory for a different number of processes.");
    }
    if (pFactory->GetProblemSize() != mProblemSize)
    {
        EXCEPTION("Cannot set from a factory for a different problem size.");
    }
    mGlobalLows.clear();
    mLo = pFactory->GetLow();
    mHi = pFactory->GetHigh();
}

DistributedVectorFactory::DistributedVectorFactory(Vec vec)
    : mPetscStatusKnown(false),
      mpOriginalFactory(nullptr)
{
    CalculateOwnership(vec);
}

DistributedVectorFactory::DistributedVectorFactory(unsigned size, PetscInt local)
    : mPetscStatusKnown(false),
      mpOriginalFactory(nullptr)
{
#ifndef NDEBUG
    CheckForPetsc();
#endif
    Vec vec = PetscTools::CreateVec(size, local);
    CalculateOwnership(vec);
    PetscTools::Destroy(vec);
}

DistributedVectorFactory::DistributedVectorFactory(DistributedVectorFactory* pOriginalFactory)
    : mPetscStatusKnown(false),
      mpOriginalFactory(pOriginalFactory)
{
    assert(mpOriginalFactory != nullptr);

    /*
     * Normally called when mpOriginalFactory->GetNumProcs() != PetscTools::GetNumProcs()
     * so ignore mpOriginalFactory->GetLocalOwnership()
     */
    Vec vec = PetscTools::CreateVec(mpOriginalFactory->GetProblemSize());

    CalculateOwnership(vec);
    PetscTools::Destroy(vec);
}

DistributedVectorFactory::DistributedVectorFactory(unsigned lo, unsigned hi, unsigned size, unsigned numProcs)
    : mLo(lo),
      mHi(hi),
      mProblemSize(size),
      mNumProcs(numProcs),
      mPetscStatusKnown(false),
      mpOriginalFactory(nullptr)
{
#ifndef NDEBUG
    CheckForPetsc();
#endif
}

DistributedVectorFactory::~DistributedVectorFactory()
{
    delete mpOriginalFactory;
}

void DistributedVectorFactory::CheckForPetsc()
{
    assert(mPetscStatusKnown==false);
    PetscBool petsc_is_initialised;
    PetscInitialized(&petsc_is_initialised);

    /*
     * Tripping this assertion means that PETSc and MPI weren't intialised.
     * A unit test should include the global fixture:
     * #include "PetscSetupAndFinalize.hpp"
     */
    assert(petsc_is_initialised);
    mPetscStatusKnown = true;
}

bool DistributedVectorFactory::IsGlobalIndexLocal(unsigned globalIndex)
{
    return (mLo<=globalIndex && globalIndex<mHi);
}

Vec DistributedVectorFactory::CreateVec()
{
    Vec vec = PetscTools::CreateVec(mProblemSize, mHi-mLo);
    return vec;
}

Vec DistributedVectorFactory::CreateVec(unsigned stride)
{
    Vec vec;
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
    //In PETSc 3.3 we cannot alter the stride (block size) after it has been created
    //There is no VecCreateMPI with block size
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetBlockSize(vec, stride);
    VecSetSizes(vec, stride*(mHi-mLo), stride*mProblemSize);
    VecSetType(vec,VECMPI);
    //VecCreateMPIWithArray(PETSC_COMM_WORLD, stride, stride*(mHi-mLo), stride*mProblemSize, PETSC_NULL/*No array*/, &vec);
#else
    VecCreateMPI(PETSC_COMM_WORLD, stride*(mHi-mLo), stride*mProblemSize, &vec);
    VecSetBlockSize(vec, stride);
#endif
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        VecSetOption(vec, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
        VecSetOption(vec, VEC_IGNORE_OFF_PROC_ENTRIES);
#endif
    return vec;
}

DistributedVector DistributedVectorFactory::CreateDistributedVector(Vec vec, bool readOnly)
{
    DistributedVector dist_vector(vec, this, readOnly);
    return dist_vector;
}

std::vector<unsigned> &DistributedVectorFactory::rGetGlobalLows()
{
    if (mGlobalLows.size() != PetscTools::GetNumProcs())
    {
        assert( mGlobalLows.empty());
        mGlobalLows.resize(PetscTools::GetNumProcs());

        // Exchange data
        MPI_Allgather( &mLo, 1, MPI_UNSIGNED, &mGlobalLows[0], 1, MPI_UNSIGNED, PETSC_COMM_WORLD);
      }

    return  mGlobalLows;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DistributedVectorFactory)
