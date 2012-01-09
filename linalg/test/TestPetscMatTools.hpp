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

#ifndef TESTPETSCMATTOOLS_HPP_
#define TESTPETSCMATTOOLS_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscMatTools.hpp" // Includes Ublas so must come before PETSc

#include "PetscSetupAndFinalize.hpp"

/**
 * Tests methods in the PETSc helper class PetscMatTools that are not directly related to LinearSystem.
 */
class TestPetscVecTools : public CxxTest::TestSuite
{
public:
    void TestEqualityCheck() throw (Exception)
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

        MatDestroy(matrix1);
        MatDestroy(matrix2);
    }

    void TestSymmetryCheck() throw (Exception)
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

        MatDestroy(matrix);
    }
};

#endif /*TESTPETSCMATTOOLS_HPP_*/
