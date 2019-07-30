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

#ifndef REPLICATABLEVECTOR_HPP_
#define REPLICATABLEVECTOR_HPP_

#include <vector>
#include <petscvec.h>

/**
 * Helper class for replicating a PETSc vector.
 */
class ReplicatableVector
{
private:

    double* mpData;     /**< The wrapped PeTSc vector. */
    unsigned mSize;     /**< The length of the vector. */
    VecScatter mToAll;   /**< Variable holding information for replicating a PETSc vector. */
    Vec mReplicated;     /**< Vector to hold concentrated copy of replicated vector. */

    /**
     * Clear data. Used in resize method and destructor.
     */
    void RemovePetscContext();

public:

    /**
     * Default constructor.
     * Note that the vector will need to be resized before it can be used.
     */
    ReplicatableVector();

    /**
     * Constructor taking in PETSc vector, which is immediately
     * replicated into the internal data
     *
     * @param vec a PETSc vector
     */
    ReplicatableVector(Vec vec);

    /**
     * Constructor to make a vector of given size.
     *
     * @param size the size of the vector
     */
    ReplicatableVector(unsigned size);

    /**
     * Default destructor.
     * Remove PETSc context.
     */
    ~ReplicatableVector();

    /**
     * @return the size of the vector.
     */
    unsigned GetSize();

    /**
     * Resize the vector.
     *
     * @param size  The number of elements to allocate memory for.
     */
    void Resize(unsigned size);

    /**
     * Access the vector.
     *
     * @param index the index of the vector to return
     * @return reference to component of the vector
     */
    double& operator[](unsigned index);

    /**
     * Replicate this vector over all processes.
     *
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes.
     *
     * @param lo  The start of our ownership range
     * @param hi  One past the end of our ownership range
     */
    void Replicate(unsigned lo, unsigned hi);

    /**
     * Replicate the given PETSc vector over all processes.
     *
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes, storing it in this object.
     *
     * Our data vector will automatically be resized to fit the whole PETSc vector.
     *
     * @param vec  The PETSc vector to replicate.
     */
    void ReplicatePetscVector(Vec vec);
};

#endif /*REPLICATABLEVECTOR_HPP_*/
