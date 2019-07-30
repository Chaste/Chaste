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

#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"

bool DistributedVector::IsGlobalIndexLocal(unsigned globalIndex)
{
    return (mLo<=globalIndex && globalIndex<mHi);
}

DistributedVector::DistributedVector(Vec vec, DistributedVectorFactory* pFactory, bool readOnly)
    : mVec(vec),
      mpFactory(pFactory),
      mReadOnly(readOnly)
{
    assert(pFactory != nullptr);

    // Set local copies of problem size, etc.
    mProblemSize = pFactory->GetProblemSize();
    mLo = pFactory->GetLow();
    mHi = pFactory->GetHigh();

    if (mReadOnly)
    {
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        // Request read-only access properly
        VecGetArrayRead(vec, (const PetscScalar**)&mpVec);
#else
        // Silently ignore the request and hope nothing breaks
        VecGetArray(vec, &mpVec);
#endif
    }
    else
    {
        VecGetArray(vec, &mpVec);
    }

    // Set mSizeMultiplier by reading the vec size.
    PetscInt size;
    VecGetSize(vec, &size);
    mSizeMultiplier = (unsigned) size / mProblemSize;
    assert ((mSizeMultiplier * mProblemSize) == (unsigned)size);
}

double& DistributedVector::operator[](unsigned globalIndex)
{
    assert(mSizeMultiplier == 1);
    if (mLo<=globalIndex && globalIndex<mHi)
    {
        return mpVec[globalIndex - mLo];
    }
    throw DistributedVectorException();
}

double& DistributedVector::operator[](Iterator index)
{
    assert(mSizeMultiplier==1);
    return mpVec[index.Local];
}

void DistributedVector::Restore()
{
    assert(mReadOnly == false);
    VecRestoreArray(mVec, &mpVec);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
    /**
     * mpVec is NULL after this function call
     */
    VecGetArray(mVec, &mpVec);
#endif

}

// Iterator class

bool DistributedVector::Iterator::operator!=(const Iterator& rOther)
{
   return(Global != rOther.Global);
}

DistributedVector::Iterator& DistributedVector::Iterator::operator++()
{
    Local++;
    Global++;
    return(*this);
}

// Iterator creation

DistributedVector::Iterator DistributedVector::Begin()
{
    Iterator index;
    index.Local = 0;
    index.Global = mLo;
    return index;
}

DistributedVector::Iterator DistributedVector::End()
{
    Iterator index;
    index.Local = mHi-mLo;
    index.Global = mHi;
    return index;
}
