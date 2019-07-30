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


#ifndef DISTRIBUTEDVECTOR_HPP_
#define DISTRIBUTEDVECTOR_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>
#include <cassert>

#include "DistributedVectorException.hpp"

class DistributedVectorFactory;

/**
 * Gives access to the local portion of a PETSc vector via an iterator.
 *
 * It also provides two nested classes for accessing vectors with particular
 * memory layouts: striped and chunked.
 */
class DistributedVector
{
private:
    friend class TestDistributedVector;

    // Data global to all vectors

    /** The first entry owned by the current processor. */
    unsigned mLo;

    /** One above the last entry owned by the current processor. */
    unsigned mHi;

    /** The problem size, i.e. the length of the vector of unknowns. */
    unsigned mProblemSize;

    // Data local to a single vector

    /** How much bigger this vector is than the problem size. */
    unsigned mSizeMultiplier;

    /** The underlying PETSc vector. */
    Vec mVec;

    /** The local part of the underlying PETSc vector. */
    double* mpVec;

    /**
     * Pointer to the factory that created this DistributedVector.
     * Gives access to local and global sizes
     */
    DistributedVectorFactory* mpFactory;

    /** A hint to PETSc that we only want to read from the vector (and never write).
     *  Note that it may still be possible to alter the content of the vector in an
     *  inconsistent way -- the user is responsible for ensuring that this doesn't happen.
     *  The only makes sense with PETSc 3.2 and above.
     */
    bool mReadOnly;


public:

    /**
     * Test if the given global index is owned by the current process, i.e. is local to it.
     *
     * @param globalIndex
     * @return true if the global index can be accessed by this process
     */
    bool IsGlobalIndexLocal(unsigned globalIndex);

    /**
     * Constructor.
     * This class represents the portion of a distributed PETSc vector on this process.
     *
     * Note that this class does NOT take over responsibility for destroying the Vec.
     *
     * @param vec PETSc vector of which this class shall be a portion
     * @param readOnly A suggestion to PETSc that we are not going to change the content of the Vec
     *        (Note that it may still be possible to alter the content of the vector in an
     *         inconsistent way -- the user is responsible for ensuring that this doesn't happen.)     *
     * @param pFactory pointer to the DistributedVectorFactory used to create this vector
     */
    DistributedVector(Vec vec, DistributedVectorFactory* pFactory, bool readOnly=false);

    /**
     * @return #mHi - The next index above the top one owned by the process.
     */
    unsigned GetHigh() const
    {
        return mHi;
    }

    /**
     * @return #mLo - The lowest index owned by the process.
     */
    unsigned GetLow() const
    {
        return mLo;
    }

    /**
     * @return  the factory used to create this vector.
     */
    DistributedVectorFactory* GetFactory()
    {
        return mpFactory;
    }

    /**
     * @param globalIndex
     * @return value of distributed vector at globalIndex
     * Do not use if stride>1.
     * For use in tests.
     * Will throw a DistributedVectorException if the specified element is not on this process.
     */
    double& operator[](unsigned globalIndex);

    /**
     * Store elements that have been written to
     * back into the PETSc vector. Call after you have finished writing.
     * It appears that you do not need to call this if you only read from the vector.
     *
     * Calling this method when the #mReadOnly flag is set results in error (as a
     * reminder that there should be no changes to the data in the original Vec).
     */
    void Restore();

    /**
     * Iterator class allows one to iterate over the elements of the distributed
     * vector on this process.
     */
    class Iterator
    {
    public:
        unsigned Local;  /**< Current index, local to this process. */
        unsigned Global; /**< Current index, global to the whole PETSc vector. */

        /**
         * Compare two indices for inequality.
         *
         * @param rOther
         * @return true when not equal
         */
        bool operator!=(const Iterator& rOther);

        /** Increment the iterator to the next index.
         * @return incremented operator
         */
        Iterator& operator++();
    };

    /**
     * Provide access to a particular stripe of a striped vector.
     *
     * A striped vector has multiple types of information encoded within a single
     * vector, with a layout like [x_1, y_1, z_1, x_2, y_2, z_2, ... x_n, y_n, z_n].
     * This class provides easy access to, for example, the x values.
     */
    class Stripe
    {
        unsigned mStride; /**< Number of types of information in the vector. */
        unsigned mStripe; /**< The number of this stripe within the vector starting from 0. */
        double* mpVec;    /**< The local part of the underlying PETSc vector. */
        unsigned mLo;     /**< The first entry owned by the current processor. */
        unsigned mHi;     /**< One above the last entry owned by the current processor. */
        DistributedVectorFactory* mpFactory; /**< The factory that created our parent vector. */

    public:
        /**
         * Constructor.
         *
         * @param parallelVec striped vector
         * @param stripe number of this stripe within the vector starting from 0
         */
        Stripe(DistributedVector parallelVec, unsigned stripe)
        {
            mStride = parallelVec.mSizeMultiplier;
            mStripe = stripe;
            assert(mStripe < mStride);
            mpVec = parallelVec.mpVec;
            mLo = parallelVec.GetLow();
            mHi = parallelVec.GetHigh();
            mpFactory = parallelVec.GetFactory();
        }

        /**
         * @return  the factory used to create this vector.
         */
        DistributedVectorFactory* GetFactory()
        {
            return mpFactory;
        }

        /**
         * Access a particular element of the stripe if on this processor.
         * For use in tests. Will throw a DistributedVectorException if
         * the specified element is not on this process.
         *
         * @param globalIndex index within the stripe
         * @return value of striped vector
         */
        double& operator[](unsigned globalIndex)
        {
            if (mLo <= globalIndex && globalIndex < mHi)
            {
                return mpVec[(globalIndex - mLo)*mStride + mStripe];
            }
            throw DistributedVectorException();
        }

        /**
         * @param index
         * @return value of striped distributed vector pointed to by index.
         */
        double& operator[](Iterator index)
        {
            return mpVec[index.Local*mStride + mStripe];
        }
    };

    /**
     * Provide access to a particular chunk of a chunked vector.
     *
     * A chunked vector has multiple types of information encoded within a single
     * vector, with a layout like [x_1, x_2, ..., x_n, y_1, y_2, ... y_n].
     * This class provides easy access to, for example, the x values.
     */
    class Chunk
    {
        unsigned mOffset; /**< The start of this chunk within the locally-owned part of the vector. */
        double* mpVec;    /**< The local part of the underlying PETSc vector. */
        unsigned mLo;     /**< The first entry owned by the current processor. */
        unsigned mHi;     /**< One above the last entry owned by the current processor. */

    public:
        /**
         * Constructor.
         *
         * @param parallelVec chunked vector
         * @param chunk number of this chunk within the vector starting from 0
         */
        Chunk(DistributedVector parallelVec, unsigned chunk)
        {
            assert(chunk < parallelVec.mSizeMultiplier);
            mLo = parallelVec.GetLow();
            mHi = parallelVec.GetHigh();
            mOffset = chunk * (mHi - mLo);
            mpVec = parallelVec.mpVec;
        }

        /**
         * Access a particular element of the chunk if on this processor.
         * For use in tests. Will throw a DistributedVectorException if
         * the specified element is not on this process.
         *
         * @param globalIndex index within the chunk
         * @return value of striped vector
         */
        double& operator[](unsigned globalIndex)
        {
            if (mLo <= globalIndex && globalIndex < mHi)
            {
                //localIndex = globalIndex - mLo
                return mpVec[mOffset + globalIndex - mLo];
            }
            throw DistributedVectorException();
         }

        /**
         * @param index
         * @return value of striped distributed vector pointed to by index.
         */
        double& operator[](Iterator index)
        {
            return mpVec[mOffset + index.Local];
        }
    };

    /**
     * @return iterator pointing to the first element of the distributed
     * vector on this process
     */
    Iterator Begin();

    /**
     * @return iterator pointing to one past the last element of the distributed
     * vector on this process
     */
    Iterator End();

    /**
     * @param index
     * @return value of distributed vector pointed to by index.
     * Do not use if stride>1.
     */
    double& operator[](Iterator index);
};

#endif /*DISTRIBUTEDVECTOR_HPP_*/
