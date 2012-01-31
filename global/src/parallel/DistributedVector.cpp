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

#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"

bool DistributedVector::IsGlobalIndexLocal(unsigned globalIndex)
{
    return (mLo<=globalIndex && globalIndex<mHi);
}

DistributedVector::DistributedVector(Vec vec, DistributedVectorFactory* pFactory)
    : mVec(vec),
      mpFactory(pFactory)
{
    assert(pFactory != NULL);

    // Set local copies of problem size, etc.
    mProblemSize = pFactory->GetProblemSize();
    mLo = pFactory->GetLow();
    mHi = pFactory->GetHigh();

    // Set mSizeMultiplier by reading the vec size.
    VecGetArray(vec, &mpVec);
    PetscInt size;
    VecGetSize(vec, &size);
    mSizeMultiplier = (unsigned) size / mProblemSize;
    assert ((mSizeMultiplier * mProblemSize) == (unsigned)size);
}

double& DistributedVector::operator[](unsigned globalIndex) throw (DistributedVectorException)
{
    assert(mSizeMultiplier == 1);
    if (mLo<=globalIndex && globalIndex<mHi)
    {
        return mpVec[globalIndex - mLo];
    }
    throw DistributedVectorException();
}

double& DistributedVector::operator[](Iterator index) throw (DistributedVectorException)
{
    assert(mSizeMultiplier==1);
    return mpVec[index.Local];
}

void DistributedVector::Restore()
{
    VecRestoreArray(mVec, &mpVec);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
    /** \todo #1994
     * mpVec is NULL after this function call
     * WARNING: VecRestoreArray function was called in many other places? should be checked - Arash
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
