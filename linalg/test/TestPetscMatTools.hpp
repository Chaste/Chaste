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

#ifndef TESTPETSCMATTOOLS_HPP_
#define TESTPETSCMATTOOLS_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscMatTools.hpp" // Includes Ublas so must come before PETSc

#include "PetscSetupAndFinalize.hpp"

/**
 * Tests methods in the PETSc helper class PetscMatTools that are not directly related to LinearSystem.
 */
class TestPetscMatTools : public CxxTest::TestSuite
{
public:
    void TestEqualityCheck()
    {
        Mat matrix1;
        Mat matrix2;
        const unsigned size = 10u;
        PetscTools::SetupMat(matrix1, size, size, size);
        PetscTools::SetupMat(matrix2, size, size, size);
        TS_ASSERT_EQUALS(PetscMatTools::GetSize(matrix1), size);
        TS_ASSERT_EQUALS(PetscMatTools::GetSize(matrix2), size);

        // Note: not efficient, but checks SetElement guards for ownership
        for (unsigned row=0; row<size; row++)
        {
            for (unsigned col=0; col<size; col++)
            {
                PetscMatTools::SetElement(matrix1, row, col, (double) size*row+col+1);
                PetscMatTools::SetElement(matrix2, row, col, (double) size*row+col+1);
            }
        }
        PetscMatTools::Finalise(matrix1);
        PetscMatTools::Finalise(matrix2);

        TS_ASSERT(PetscMatTools::CheckEquality(matrix1, matrix2));

        PetscMatTools::AddToElement(matrix2, 1, 1, 2e-10);
        PetscMatTools::Finalise(matrix2);
        TS_ASSERT(!PetscMatTools::CheckEquality(matrix1, matrix2));
        TS_ASSERT(PetscMatTools::CheckEquality(matrix1, matrix2, 1e-9));

        PetscTools::Destroy(matrix1);
        PetscTools::Destroy(matrix2);
    }

    void TestSymmetryCheck()
    {
        Mat matrix;
        const unsigned size = 10u;
        PetscTools::SetupMat(matrix, size, size, size);
        // Note: not efficient, but checks SetElement guards for ownership
        for (unsigned row=0; row<size; row++)
        {
            for (unsigned col=0; col<row; col++)
            {
                PetscMatTools::SetElement(matrix, row, col, (double) size*row+col+1);
                PetscMatTools::SetElement(matrix, col, row, (double) size*row+col+1);
            }
        }
        PetscMatTools::Finalise(matrix);
        TS_ASSERT(PetscMatTools::CheckSymmetry(matrix));

        PetscMatTools::AddToElement(matrix, 1, 2, 2e-10);
        PetscMatTools::Finalise(matrix);
        TS_ASSERT(!PetscMatTools::CheckSymmetry(matrix));
        TS_ASSERT(PetscMatTools::CheckSymmetry(matrix, 1e-9));

        PetscTools::Destroy(matrix);
    }

    void TestZeroRowsAndColumnsWithValueOnDiagonal()
    {
        Mat matrix;
        const unsigned size = 5u;
        PetscTools::SetupMat(matrix, size, size, size);
        for (unsigned row=0; row<size; row++)
        {
            for (unsigned col=0; col<size; col++)
            {
                PetscMatTools::SetElement(matrix, row, col, 2.78);
            }
        }
        PetscMatTools::Finalise(matrix);

        std::vector<unsigned> rows;
        rows.push_back(1);
        rows.push_back(3);
        rows.push_back(4);

        PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(matrix, rows, 3.41);

        double correct_mat[5][5] = { {2.78,    0, 2.78,    0,     0},
                                     {   0, 3.41,    0,    0,     0},
                                     {2.78,    0, 2.78,    0,     0},
                                     {   0,    0,    0, 3.41,     0},
                                     {   0,    0,    0,    0,  3.41}    };

        PetscInt lo, hi;
        PetscMatTools::GetOwnershipRange(matrix, lo, hi);

        for (int i=lo; i<hi; i++)
        {
            for (unsigned j=0; j<size; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(matrix, i, j), correct_mat[i][j], 1e-12);
            }
        }

        PetscTools::Destroy(matrix);
    }


    void TestTurnOffVariableAllocationError()
    {
        Mat matrix;
        const unsigned size = 5u;
        PetscTools::SetupMat(matrix, size, size, 1);

        //Shows that TurnOffVariableAllocationError doesn't *have* to be called straight after SetupMat,
        //but must be called before an allocation problem occurs.
        PetscMatTools::SetElement(matrix, 0, 0, 4.0);
        PetscMatTools::Finalise(matrix);

        PetscMatTools::TurnOffVariableAllocationError(matrix);

        for (unsigned row=0; row<size; row++)
        {
            for (unsigned col=0; col<size; col++)
            {
                PetscMatTools::SetElement(matrix, row, col, 2.78);
            }
        }

        PetscMatTools::Finalise(matrix);

        PetscInt lo, hi;
        PetscMatTools::GetOwnershipRange(matrix, lo, hi);
        if (lo<=3 && 3<hi)
        {
            TS_ASSERT_DELTA(PetscMatTools::GetElement(matrix, 3, 3), 2.78, 1e-6);
        }

        PetscTools::Destroy(matrix);
    }
};

#endif /*TESTPETSCMATTOOLS_HPP_*/
