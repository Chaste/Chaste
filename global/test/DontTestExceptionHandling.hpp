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

#ifndef _TESTEXCEPTIONHANDLING_HPP_
#define _TESTEXCEPTIONHANDLING_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * This tests that cxxtest handles exceptions thrown in tests gracefully.
 *
 * The tests are supposed to fail, so aren't run normally.
 */
class TestExceptionHandling : public CxxTest::TestSuite
{
public:

    void TestThrowingWithPetsc()
    {
        Vec test_vecc = PetscTool::CreateVec(20);

        VecAssemblyBegin(test_vec);
        VecAssemblyEnd(test_vec);

        EXCEPTION("Will cxxtest be nice if we do PETSc things?");
    }

    void ThrowExceptionMethod()
    {
        EXCEPTION("Exception thrown from method.");
    }

    void TestThrowingAnExceptionInATest()
    {
        EXCEPTION("Will cxxtest be nice I wonder?");
    }

    void TestCatchingExceptionWithCxxtest()
    {
        TS_ASSERT_THROWS_THIS(EXCEPTION("Will cxxtest be nice I wonder?"),"Will cxxtest be nice I wonder?");
        TS_ASSERT_THROWS_NOTHING(EXCEPTION("Will cxxtest be nice I wonder?"));
    }
};

#endif //_TESTEXCEPTIONHANDLING_HPP_
