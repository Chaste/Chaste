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

#include "PetscMatTools.hpp"
#include <algorithm>
#include <cassert>


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

void PetscMatTools::SetElement(Mat matrix, PetscInt row, PetscInt col, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    if (row >= lo && row < hi)
    {
        MatSetValue(matrix, row, col, value, INSERT_VALUES);
    }
}

void PetscMatTools::AddToElement(Mat matrix, PetscInt row, PetscInt col, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    if (row >= lo && row < hi)
    {
        MatSetValue(matrix, row, col, value, ADD_VALUES);
    }
}

void PetscMatTools::Finalise(Mat matrix)
{
    MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
}

void PetscMatTools::SwitchWriteMode(Mat matrix)
{
    MatAssemblyBegin(matrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FLUSH_ASSEMBLY);
}

void PetscMatTools::Display(Mat matrix)
{
    //Give full precision, scientific notation
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7) // PETSc 3.7+
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
#else
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
#endif
    MatView(matrix,PETSC_VIEWER_STDOUT_WORLD);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7) // PETSc 3.7+
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
#endif
}

void PetscMatTools::SetRow(Mat matrix, PetscInt row, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    if (row >= lo && row < hi)
    {
        PetscInt rows, cols;
        MatGetSize(matrix, &rows, &cols);
        for (PetscInt i=0; i<cols; i++)
        {
            SetElement(matrix, row, i, value);
        }
    }
}

void PetscMatTools::ZeroRowsWithValueOnDiagonal(Mat matrix, std::vector<unsigned>& rRows, double diagonalValue)
{
    Finalise(matrix);

    /*
     * Important! PETSc by default will destroy the sparsity structure for this row and deallocate memory
     * when the row is zeroed, and if there is a next timestep, the memory will have to reallocated when
     * assembly to done again. This can kill performance. The following makes sure the zeroed rows are kept.
     *
     * Note: if the following lines are called, then once MatZeroRows() is called below, there will be an
     * error if some rows have no entries at all. Hence for problems with dummy variables, Stokes flow for
     * example, the identity block needs to be added before dirichlet BCs.
     */
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 1) //PETSc 3.1 or later
    MatSetOption(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0) //PETSc 3.0
    MatSetOption(matrix, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
#else //PETSc 2.x.x
    MatSetOption(matrix, MAT_KEEP_ZEROED_ROWS);
#endif

    PetscInt* rows = new PetscInt[rRows.size()];
    for (unsigned i=0; i<rRows.size(); i++)
    {
        rows[i] = rRows[i];
    }
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    IS is;
    ISCreateGeneral(PETSC_COMM_WORLD, rRows.size(), rows, &is);
    MatZeroRows(matrix, is, &diagonalValue);
    ISDestroy(PETSC_DESTROY_PARAM(is));
    /*
     *
[2]PETSC ERROR: MatMissingDiagonal_SeqAIJ() line 1011 in src/mat/impls/aij/seq/aij.c
[2]PETSC ERROR: PETSc has generated inconsistent data!
[2]PETSC ERROR: Matrix is missing diagonal number 15!
[2]PETSC ERROR: MatILUFactorSymbolic_SeqAIJ() line 906 in src/mat/impls/aij/seq/aijfact.c
     *
     *
    // There appears to be a problem with MatZeroRows not setting diagonals correctly
    // While we are supporting PETSc 2.2, we have to do this the slow way

    AssembleFinal(matrix);
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);
    PetscInt size=GetSize(matrix);
    ///assert(rRows.size() == 1);
    for (unsigned index=0; index<rRows.size(); index++)
    {
        PetscInt row = rRows[index];
        if (row >= lo && row < hi)
        {
            std::vector<unsigned> non_zero_cols;
            // This row is local, so zero it.
            for (PetscInt column = 0; column < size; column++)
            {
                if (GetElement(matrix, row, column) != 0.0)
                {
                    non_zero_cols.push_back(column);
                }
            }
            // Separate "gets" from "sets" so that we don't have to keep going into "assembled" mode
            for (unsigned i=0; i<non_zero_cols.size();i++)
            {
                SetElement(matrix, row, non_zero_cols[i], 0.0);
            }
            // Set the diagonal
            SetElement(matrix, row, row, diagonalValue);
        }
        // Everyone communicate after row is finished
        AssembleFinal(matrix);
    }
    */
#elif (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
    MatZeroRows(matrix, rRows.size(), rows, diagonalValue , nullptr, nullptr);
#else
    MatZeroRows(matrix, rRows.size(), rows, diagonalValue);
#endif
    delete [] rows;
}


void PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(Mat matrix, std::vector<unsigned> rowColIndices, double diagonalValue)
{
    Finalise(matrix);

    // sort the vector as we will be repeatedly searching for entries in it
    std::sort(rowColIndices.begin(), rowColIndices.end());

    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);
    unsigned size = hi-lo;

    std::vector<unsigned>* cols_to_zero_per_row = new std::vector<unsigned>[size];

    // collect the columns to be zeroed, for each row. We don't zero yet as we would
    // then have to repeatedly call Finalise before each MatGetRow
    for (PetscInt row = lo; row < hi; row++)
    {
        // get all the non-zero cols for this row
        PetscInt num_cols;
        const PetscInt* cols;
        MatGetRow(matrix, row, &num_cols, &cols, PETSC_NULL);

        // see which of these cols are in the list of cols to be zeroed
        for (PetscInt i=0; i<num_cols; i++)
        {
            if (std::binary_search(rowColIndices.begin(), rowColIndices.end(), cols[i]))
            {
                cols_to_zero_per_row[row-lo].push_back(cols[i]);
            }
        }

        // this must be called for each MatGetRow
        MatRestoreRow(matrix, row, &num_cols, &cols, PETSC_NULL);
    }

    // Now zero columns of each row
    for (PetscInt row = lo; row < hi; row++)
    {
        unsigned num_cols_to_zero_this_row = cols_to_zero_per_row[row-lo].size();

        if (num_cols_to_zero_this_row>0)
        {
            PetscInt* cols_to_zero = new PetscInt[num_cols_to_zero_this_row];
            double* zeros = new double[num_cols_to_zero_this_row];
            for (unsigned i=0; i<num_cols_to_zero_this_row; i++)
            {
                cols_to_zero[i] = cols_to_zero_per_row[row-lo][i];
                zeros[i] = 0.0;
            }

            PetscInt rows[1];
            rows[0] = row;
            MatSetValues(matrix, 1, rows, num_cols_to_zero_this_row, cols_to_zero, zeros, INSERT_VALUES);
            delete [] cols_to_zero;
            delete [] zeros;
        }
    }

    delete [] cols_to_zero_per_row;

    // Now zero the rows and add the diagonal entries
    ZeroRowsWithValueOnDiagonal(matrix, rowColIndices, diagonalValue);
}

void PetscMatTools::ZeroColumn(Mat matrix, PetscInt col)
{
    Finalise(matrix);

    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    // Determine which rows in this column are non-zero (and therefore need to be zeroed)
    std::vector<unsigned> nonzero_rows;
    for (PetscInt row = lo; row < hi; row++)
    {
        if (GetElement(matrix, row, col) != 0.0)
        {
            nonzero_rows.push_back(row);
        }
    }

    // Set those rows to be zero by calling MatSetValues
    unsigned size = nonzero_rows.size();
    PetscInt* rows = new PetscInt[size];
    PetscInt cols[1];
    double* zeros = new double[size];

    cols[0] = col;

    for (unsigned i=0; i<size; i++)
    {
        rows[i] = nonzero_rows[i];
        zeros[i] = 0.0;
    }

    MatSetValues(matrix, size, rows, 1, cols, zeros, INSERT_VALUES);
    delete [] rows;
    delete [] zeros;
}

void PetscMatTools::Zero(Mat matrix)
{
    MatZeroEntries(matrix);
}

unsigned PetscMatTools::GetSize(Mat matrix)
{
    PetscInt rows, cols;

    MatGetSize(matrix, &rows, &cols);
    assert(rows == cols);
    return (unsigned) rows;
}

void PetscMatTools::GetOwnershipRange(Mat matrix, PetscInt& lo, PetscInt& hi)
{
    MatGetOwnershipRange(matrix, &lo, &hi);
}

double PetscMatTools::GetElement(Mat matrix, PetscInt row, PetscInt col)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    assert(lo <= row && row < hi);
    PetscInt row_as_array[1];
    row_as_array[0] = row;
    PetscInt col_as_array[1];
    col_as_array[0] = col;

    double ret_array[1];

    MatGetValues(matrix, 1, row_as_array, 1, col_as_array, ret_array);

    return ret_array[0];
}

void PetscMatTools::SetOption(Mat matrix, MatOption option)
{
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
    MatSetOption(matrix, option, PETSC_TRUE);
#else
    MatSetOption(matrix, option);
#endif
}

Vec PetscMatTools::GetMatrixRowDistributed(Mat matrix, unsigned rowIndex)
{
    /*
     * We need to make sure that lhs_ith_row doesn't ignore off processor entries when assembling,
     * otherwise the VecSetValues call a few lines below will not work as expected.
     */

    PetscInt lo, hi;
    PetscMatTools::GetOwnershipRange(matrix, lo, hi);
    unsigned size = PetscMatTools::GetSize(matrix);

    Vec mat_ith_row = PetscTools::CreateVec(size, hi-lo, false);

    PetscInt num_entries;
    const PetscInt* column_indices;
    const PetscScalar* values;

    bool am_row_owner = (PetscInt)rowIndex >= lo && (PetscInt)rowIndex < hi;

    /*
     * Am I the owner of the row? If so get the non-zero entries and add them lhs_ith_row.
     * In parallel, VecAssembly{Begin,End} will send values to the rest of processors.
     */
    if (am_row_owner)
    {
        MatGetRow(matrix, rowIndex, &num_entries, &column_indices, &values);
        VecSetValues(mat_ith_row, num_entries, column_indices, values, INSERT_VALUES);
    }

    VecAssemblyBegin(mat_ith_row);
    VecAssemblyEnd(mat_ith_row);

    if (am_row_owner)
    {
        MatRestoreRow(matrix, rowIndex, &num_entries, &column_indices, &values);
    }

    return mat_ith_row;
}


bool PetscMatTools::CheckEquality(const Mat mat1, const Mat mat2, double tol)
{
    Mat y;
    MatDuplicate(mat2, MAT_COPY_VALUES, &y);

    double minus_one = -1.0;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    // MatAYPX(*a, X, Y) does  Y = X + a*Y.
    MatAYPX(&minus_one, mat1, y);
#elif (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 1) //PETSc 2.3.1
    // MatAYPX( Y, a, X) does Y = a*Y + X.
    MatAYPX(y, minus_one, mat1);
#else
    // MatAYPX( Y, a, X, structure) does Y = a*Y + X.
    MatAYPX(y, minus_one, mat1, DIFFERENT_NONZERO_PATTERN);
#endif
    PetscReal norm;
    MatNorm(y, NORM_INFINITY, &norm);
    PetscTools::Destroy(y);

    return (norm < tol);
}

bool PetscMatTools::CheckSymmetry(const Mat matrix, double tol)
{
    Mat trans;
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
    MatTranspose(matrix, MAT_INITIAL_MATRIX, &trans);
#else
    MatTranspose(matrix, &trans);
#endif
    bool is_symmetric = PetscMatTools::CheckEquality(matrix, trans, tol);
    PetscTools::Destroy(trans);
    return is_symmetric;
}

void PetscMatTools::TurnOffVariableAllocationError(Mat matrix)
{
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
    MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif
}

