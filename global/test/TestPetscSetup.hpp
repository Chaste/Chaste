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

#ifndef _TESTPETSCSETUP_HPP_
#define _TESTPETSCSETUP_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"
#include "Warnings.hpp"
#include "IsNan.hpp"

/**
 * This tests that the initialisation of PETSc does something.
 */
class TestPetscSetup : public CxxTest::TestSuite
{
public:

    void TestPetscIsThere()
    {
        PetscTruth is_there;
        PetscInitialized(&is_there);
        TS_ASSERT( is_there == PETSC_TRUE );
    }

    void TestPetscExceptions()
    {
        /*
         * Note we could test with TS_ASSERT_THROWS_THIS() but PetscException
         * includes line numbers so it isn't very robust.
         */
        int err = 0;
        TS_ASSERT_THROWS_NOTHING(PETSCEXCEPT(err));

        Vec v;
        err = VecCreateMPI(PETSC_COMM_WORLD, 2, 1, &v);
        VecDestroy(v);
        //#define PETSC_ERR_ARG_WRONGSTATE   73   /* object in argument is in wrong */
        //#define PETSC_ERR_ARG_INCOMP       75   /* two arguments are incompatible */
        TS_ASSERT_EQUALS(err, PETSC_ERR_ARG_INCOMP);
        TS_ASSERT_THROWS(PETSCEXCEPT(err), Exception);

        err=PETSC_ERR_FILE_OPEN;
        //#define PETSC_ERR_FILE_OPEN        65   /* unable to open file */
        TS_ASSERT_EQUALS(err, PETSC_ERR_FILE_OPEN);
        TS_ASSERT_THROWS(PETSCEXCEPT(err), Exception);

        // See if we can do it without a temporary
        TS_ASSERT_THROWS(
            PETSCEXCEPT(VecCreateMPI(PETSC_COMM_WORLD, 2, 1, &v)), Exception);
        VecDestroy(v);

        // This test give back an "unknown error" message
        TS_ASSERT_THROWS( PETSCEXCEPT(-3), Exception);
    }

    void TestKspExceptionsForCoverage()
    {
        TS_ASSERT_THROWS_NOTHING( KSPEXCEPT(2));

        /*
         * These next few lines are designed to force the coverage test to pass.
         * Some are hard to throw in normal circumstances --
         * "Unknown KSP error code" ought never to be thrown.
         */
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_ITS), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_DTOL), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN_BICG), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_NONSYMMETRIC), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_INDEFINITE_PC), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(-735827), Exception );


        KSPWARNIFFAILED(KSP_DIVERGED_ITS);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::QuietDestroy();
    }

    void TestDivideOneByZero() throw(Exception)
    {
        double one = 1.0;
        double zero = 0.0;
        double ans;
#ifdef TEST_FOR_FPE
// If we are testing for divide-by-zero, then this will throw an exception
        //TS_ASSERT_THROWS_ANYTHING(ans = one / zero);
        ans=zero*one;//otherwise compiler would complain
#else
// If we aren't testing for it, then there will be no exception
        TS_ASSERT_THROWS_NOTHING(ans = one / zero);
        double negative_infinity=std::numeric_limits<double>::infinity();
        TS_ASSERT_EQUALS(ans, negative_infinity);
#endif
    }

    void TestDivideZeroByZero() throw(Exception)
    {
        double zero = 0.0;
        double ans;
#ifdef TEST_FOR_FPE
// If we are testing for divide-by-zero, then this will throw an exception
        //TS_ASSERT_THROWS_ANYTHING(ans = zero / zero);
        ans=zero;//otherwise compiler would complain
#else
// If we aren't testing for it, then there will be no exception
        TS_ASSERT_THROWS_NOTHING(ans = zero / zero);
        TS_ASSERT(std::isnan(ans));
#endif
    }
};

#endif // _TESTPETSCSETUP_HPP_
