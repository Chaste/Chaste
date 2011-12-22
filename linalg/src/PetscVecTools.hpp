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

#ifndef _PETSCVECTOOLS_HPP_
#define _PETSCVECTOOLS_HPP_

#include "UblasVectorInclude.hpp" // needs to be 'first'

#include <petscvec.h>

/**
 * A collection of static methods for working with PETSc vectors.
 */
class PetscVecTools
{
public:

    /**
     * Do parallel communication required to get the vector in a good state for further operations.
     * This is a wrapper to PETSc functions like VecAssemblyBegin.
     * @param vector  the vector to assemble
     */
    static void Finalise(Vec vector);

    /**
     * Display the given vector.
     *
     * @param vector  the vector to display
     */
    static void Display(Vec vector);

    /**
     * Zero all entries of a given vector.
     *
     * @param vector  the vector to be zero
     */
    static void Zero(Vec vector);

    /**
     * Set an element of a vector to a given value.
     *
     * @param vector  the vector to modify
     * @param row  the row index
     * @param value  the value to set this entry
     */
    static void SetElement(Vec vector, PetscInt row, double value);

    /**
     * Add a value to an element of a vector.
     *
     * @param vector  the vector to modify
     * @param row  the row index
     * @param value  the value to set this entry
     */
    static void AddToElement(Vec vector, PetscInt row, double value);

    /**
     * Get the size of a vector.
     *
     * @param vector  the vector
     */
    static unsigned GetSize(Vec vector);

    /**
     * Get this process's ownership range of the contents of the vector.
     *
     * @param vector  the vector
     * @param lo  lowest index owned by this process
     * @param hi  highest index owned by this process
     */
    static void GetOwnershipRange(Vec vector, PetscInt& lo, PetscInt& hi);

    /**
     * Return an element of a vector.
     * May only be called for elements you own.
     *
     * @param vector  the vector
     * @param row  the row index
     */
    static double GetElement(Vec vector, PetscInt row);

    /**
     * Computes y += ax, using the PETSc method VecAXPY (with appropriate arguments for the PETSc version).
     *
     * @param y the vector which is added to
     * @param x the vector which is scaled and added to y
     * @param scaleFactor the value 'a' above, the factor x is multiplied by.
     */
    static void AddScaledVector(Vec y, Vec x, double scaleFactor);

    /**
     * Scale the given vector. Calls VecScale (using the appropriate arguments for the PETSc version).
     * @param vector the vector
     * @param scaleFactor the scale factor
     */
    static void Scale(Vec vector, double scaleFactor);

    /**
     * Calls the PETSc function VecWAXPY (using the appropriate arguments for the PETSc version),
     * which does w = ax+y, where x,y,w are distinct vectors and a is scalar.
     * @param w the result vector
     * @param a the scale factor
     * @param x the scale vector
     * @param y the other vector
     */
    static void WAXPY(Vec w, double a, Vec x, Vec y);

    /**
     * Add multiple values to a vector.
     *
     * @param vector  the vector to modify
     * @param vectorIndices mapping from index of the ublas vector (see param below)
     *  to index of the vector of this linear system
     * @param smallVector Ublas vector containing the values to be added
     *
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
    template<size_t VECTOR_SIZE>
    static void AddMultipleValues(Vec vector, unsigned* vectorIndices, c_vector<double, VECTOR_SIZE>& smallVector)
    {
        PetscInt indices_owned[VECTOR_SIZE];
        PetscInt num_indices_owned = 0;
        PetscInt global_row;
        PetscInt lo, hi;
        GetOwnershipRange(vector, lo, hi);

        for (unsigned row = 0; row<VECTOR_SIZE; row++)
        {
            global_row = vectorIndices[row];
            if (global_row >=lo && global_row <hi)
            {
                indices_owned[num_indices_owned++] = global_row;
            }
        }

        if (num_indices_owned == VECTOR_SIZE)
        {
            VecSetValues(vector,
                         num_indices_owned,
                         indices_owned,
                         smallVector.data(),
                         ADD_VALUES);
        }
        else
        {
            /*
             * We need continuous data, if some of the rows do not belong
             * to the processor their values are not passed to VecSetValues.
             */
            double values[VECTOR_SIZE];
            unsigned num_values_owned = 0;

            for (unsigned row = 0; row<VECTOR_SIZE; row++)
            {
                global_row = vectorIndices[row];
                if (global_row >= lo && global_row < hi)
                {
                    values[num_values_owned++] = smallVector(row);
                }
            }

            VecSetValues(vector,
                         num_indices_owned,
                         indices_owned,
                         values,
                         ADD_VALUES);
        }
    }

    /**
     * Set up scatter/gather PETSc context for splitting a bidomain-like vector with interleaved values for
     * two variables into two separate PETSc Vec containing each of them.
     *
     * @param interleavedVec Source vector with interleaved values.
     * @param rFirstVariableScatterContext Context for scattering/gathering first variable
     * @param rSecondVariableScatterContext Context for scattering/gathering second variable
     */
    static void SetupInterleavedVectorScatterGather(Vec interleavedVec, VecScatter& rFirstVariableScatterContext, VecScatter& rSecondVariableScatterContext);

    /**
     * Performs scatter operation from a bidomain-like vector of interleaved values into two separate PETSc Vec.
     *
     * @param interleavedVec Source vector with interleaved values.
     * @param firstVariableScatterContext Context for scattering/gathering first variable
     * @param firstVariableVec Destination vector for first variable
     * @param secondVariableScatterContext Context for scattering/gathering second variable
     * @param secondVariableVec Destination vector for second variable
     */
    static void DoInterleavedVecScatter(Vec interleavedVec, VecScatter firstVariableScatterContext, Vec firstVariableVec, VecScatter secondVariableScatterContext, Vec secondVariableVec);

    /**
     * Performs scatter operation from a bidomain-like vector of interleaved values into two separate PETSc Vec.
     *
     * @param interleavedVec Destination vector for interleaved values.
     * @param firstVariableScatterContext Context for scattering/gathering first variable
     * @param firstVariableVec Source vector with first variable
     * @param secondVariableScatterContext Context for scattering/gathering second variable
     * @param secondVariableVec Source vector with second variable
     */
    static void DoInterleavedVecGather(Vec interleavedVec, VecScatter firstVariableScatterContext, Vec firstVariableVec, VecScatter secondVariableScatterContext, Vec secondVariableVec);
};

#endif //_PETSCVECTOOLS_HPP_
