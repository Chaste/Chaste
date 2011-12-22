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
     * Return the size of the vector.
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
