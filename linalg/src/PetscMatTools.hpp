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

#ifndef _PETSCMATTOOLS_HPP_
#define _PETSCMATTOOLS_HPP_

#include "UblasMatrixInclude.hpp" // needs to be 'first'
#include "PetscTools.hpp"
#include <vector>
#include <petscvec.h>
#include <petscmat.h>

/**
 * A collection of static methods for working with PETSc matrices.
 */
class PetscMatTools
{
public:

    /**
     * Change one of the entries of a matrix to the specified value.
     *
     * @param matrix  the matrix
     * @param row  the row index
     * @param col  the column index
     * @param value  the value for this entry
     */
    static void SetElement(Mat matrix, PetscInt row, PetscInt col, double value);

    /**
     * Add the specified value to an entry of a matrix.
     *
     * @param matrix  the matrix
     * @param row  the row index
     * @param col  the column index
     * @param value  the value for this entry
     */
    static void AddToElement(Mat matrix, PetscInt row, PetscInt col, double value);

    /**
     * Sets up the matrix ready for use, e.g. in solving a linear system.
     * This is a wrapper to PETSc functions like MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY).
     * @param matrix  the matrix
     */
    static void Finalise(Mat matrix);

    /**
     * This must be called if switching between inserting or adding values to the matrix, to ensure all
     * processes are in sync.
     * This is a wrapper to PETSc functions like MatAssemblyBegin(matrix, MAT_FLUSH_ASSEMBLY).
     * @param matrix  the matrix
     */
    static void SwitchWriteMode(Mat matrix);

    /**
     * Display a matrix.
     *
     * @param matrix  the matrix
     */
    static void Display(Mat matrix);

    /**
     * Set all entries in a given row of a matrix to a certain value.
     * This must be called by the process who owns the row, (but other
     * processors will treat it as a null-op
     *
     * @param matrix  the matrix
     * @param row  the row index
     * @param value  the value to set each entry in this row
     */
    static void SetRow(Mat matrix, PetscInt row, double value);

    /**
     *  Zero several rows of a matrix, putting a given value in the diagonal entries.
     *
     *  *Massively* less expensive than zeroing each matrix row individually
     *
     * @param matrix  the matrix
     * @param rRows std::vector of rows to be zeroed
     * @param diagonalValue value to put in the diagonal entries (of the zeroed rows)
     */
    static void ZeroRowsWithValueOnDiagonal(Mat matrix, std::vector<unsigned>& rRows, double diagonalValue);

    /**
     * Zero several rows and columns of a matrix, putting a given value on the diagonal.
     *
     * @param matrix  the matrix
     * @param rRowColIndices A list of indices. All the rows with these indices, and all the columns
     * with these indices, will be zeroed.
     * @param diagonalValue value to put in the diagonal entries (of the zeroed rows)
     */
    static void ZeroRowsAndColumnsWithValueOnDiagonal(Mat matrix, std::vector<unsigned>& rRowColIndices, double diagonalValue);

    /**
     * Zero a column of a matrix.
     *
     * Unfortunately there is no equivalent method in Petsc, so this has to be
     * done carefully to ensure that the sparsity structure of the matrix
     * is not broken. Only owned entries which are non-zero are zeroed.
     *
     * @param matrix  the matrix
     * @param col  the column index
     */
    static void ZeroColumn(Mat matrix, PetscInt col);

    /**
     * Zero all entries of a matrix.
     * @param matrix  the matrix
     */
    static void Zero(Mat matrix);

    /**
     * Get the size of a matrix
     * @param matrix  the matrix
     */
    static unsigned GetSize(Mat matrix);

    /**
     * Get this process's ownership range of the contents of the system.
     *
     * @param matrix  the matrix
     * @param lo  lowest index owned by this process
     * @param hi  highest index owned by this process
     */
    static void GetOwnershipRange(Mat matrix, PetscInt& lo, PetscInt& hi);

    /**
     * Return an element of a matrix.
     * May only be called for elements you own.
     *
     * @param matrix  the matrix
     * @param row  the row index
     * @param col  the column index
     */
    static double GetElement(Mat matrix, PetscInt row, PetscInt col);

    /**
     * Set a PETSc matrix option to be true, using the PETSc method MatSetOption.
     *
     * @param matrix  the matrix for which to set the option
     * @param option  the option to set
     */
    static void SetOption(Mat matrix, MatOption option);

    /**
     * Returns the i-th row of the LHS matrix as a distributed PETSc Vec
     *
     * @param matrix  the matrix
     * @param rowIndex the row index
     * @return rowIndex-th row of the matrix in distributed format
     */
    static Vec GetMatrixRowDistributed(Mat matrix, unsigned rowIndex);

    /**
     * Check whether two matrices are equal to within a given tolerance.
     *
     * @param mat1  the first matrix
     * @param mat2  the second matrix
     * @param tol  the tolerance
     */
    static bool CheckEquality(const Mat mat1, const Mat mat2, double tol=1e-10);

    /**
     * Check whether a matrix is symmetrix, to within a given tolerance, by
     * checking if it is (approximately) equal to its transpose.
     *
     * Note that while there is a PETSc method MatIsSymmetric, it won't work in
     * parallel on some PETSc versions:
     * "Matrix of type <mpiaij> does not support checking for symmetric!"
     *
     * Also checking for exact equality of the transpose can break on 32bit systems.
     *
     * @param matrix  the matrix to check
     * @param tol  the tolerance
     */
    static bool CheckSymmetry(const Mat matrix, double tol=1e-10);

    /**
     * Add multiple values to a matrix.
     *
     * @param matrix  the matrix
     * @param matrixRowAndColIndices mapping from index of the ublas matrix (see param below)
     *  to index of the PETSc matrix of this linear system
     * @param rSmallMatrix Ublas matrix containing the values to be added
     *
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
    template<size_t MATRIX_SIZE>
    static void AddMultipleValues(Mat matrix, unsigned* matrixRowAndColIndices, c_matrix<double, MATRIX_SIZE, MATRIX_SIZE>& rSmallMatrix)
    {
        PetscInt matrix_row_indices[MATRIX_SIZE];
        PetscInt num_rows_owned = 0;
        PetscInt global_row;
        PetscInt lo, hi;
        GetOwnershipRange(matrix, lo, hi);

        for (unsigned row = 0; row<MATRIX_SIZE; row++)
        {
            global_row = matrixRowAndColIndices[row];
            if (global_row >= lo && global_row < hi)
            {
                matrix_row_indices[num_rows_owned++] = global_row;
            }
        }

        if ( num_rows_owned == MATRIX_SIZE)
        {
            MatSetValues(matrix,
                         num_rows_owned,
                         matrix_row_indices,
                         MATRIX_SIZE,
                         (PetscInt*) matrixRowAndColIndices,
                         rSmallMatrix.data(),
                         ADD_VALUES);
        }
        else
        {
            // We need continuous data, if some of the rows do not belong to the processor their values
            // are not passed to MatSetValues
            double values[MATRIX_SIZE*MATRIX_SIZE];
            unsigned num_values_owned = 0;
            for (unsigned row = 0; row < MATRIX_SIZE; row++)
            {
                global_row = matrixRowAndColIndices[row];
                if (global_row >= lo && global_row < hi)
                {
                    for (unsigned col = 0; col < MATRIX_SIZE; col++)
                    {
                        values[num_values_owned++] = rSmallMatrix(row, col);
                    }
                }
            }

            MatSetValues(matrix,
                         num_rows_owned,
                         matrix_row_indices,
                         MATRIX_SIZE,
                         (PetscInt*) matrixRowAndColIndices,
                         values,
                         ADD_VALUES);
        }
    }
};

#endif //_PETSCMATTOOLS_HPP_
