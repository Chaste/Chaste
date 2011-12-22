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

#include <cassert>

#include "DistributedVectorFactory.hpp"

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
      mpOriginalFactory(NULL)
{
    CalculateOwnership(vec);
}

DistributedVectorFactory::DistributedVectorFactory(unsigned size, PetscInt local)
    : mPetscStatusKnown(false),
      mpOriginalFactory(NULL)
{
#ifndef NDEBUG
    CheckForPetsc();
#endif
    Vec vec = PetscTools::CreateVec(size, local);
    CalculateOwnership(vec);
    VecDestroy(vec);
}

DistributedVectorFactory::DistributedVectorFactory(DistributedVectorFactory* pOriginalFactory)
    : mPetscStatusKnown(false),
      mpOriginalFactory(pOriginalFactory)
{
    assert(mpOriginalFactory != NULL);

    /*
     * Normally called when mpOriginalFactory->GetNumProcs() != PetscTools::GetNumProcs()
     * so ignore mpOriginalFactory->GetLocalOwnership()
     */
    Vec vec = PetscTools::CreateVec(mpOriginalFactory->GetProblemSize());

    CalculateOwnership(vec);
    VecDestroy(vec);
}

DistributedVectorFactory::DistributedVectorFactory(unsigned lo, unsigned hi, unsigned size, unsigned numProcs)
    : mLo(lo),
      mHi(hi),
      mProblemSize(size),
      mNumProcs(numProcs),
      mPetscStatusKnown(false),
      mpOriginalFactory(NULL)
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
    PetscTruth petsc_is_initialised;
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
    VecCreateMPI(PETSC_COMM_WORLD, stride*(mHi-mLo), stride*mProblemSize, &vec);
    return vec;
}

DistributedVector DistributedVectorFactory::CreateDistributedVector(Vec vec)
{
    DistributedVector dist_vector(vec, this);
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
