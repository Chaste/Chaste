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

#ifndef _TESTPETSCSETUP_HPP_
#define _TESTPETSCSETUP_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"
#include "PetscTools.hpp"
#include "Warnings.hpp"
#include "IsNan.hpp"

class TestPetscSetup : public CxxTest::TestSuite
{
public:

    void TestPetscIsThere()
    {
        PetscBool is_there;
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
        PetscTools::Destroy(v);
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
        PetscTools::Destroy(v);

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
        TS_ASSERT_THROWS_CONTAINS( KSPEXCEPT(KSP_DIVERGED_ITS), "DIVERGED_ITS in function \'TestKspExceptionsForCoverage\' on line");
        // The next one is deliberately fragile because it contains the line number in this test suite (to check that the line number is output correctly).
        TS_ASSERT_THROWS_THIS( KSPEXCEPT(KSP_DIVERGED_DTOL),  "DIVERGED_DTOL in function \'TestKspExceptionsForCoverage\' on line 102 of file ./global/test/TestPetscSetup.hpp");
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN_BICG), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_NONSYMMETRIC), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(KSP_DIVERGED_INDEFINITE_PC), Exception );
        TS_ASSERT_THROWS( KSPEXCEPT(-735827), Exception );


        KSPWARNIFFAILED(KSP_DIVERGED_ITS);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::QuietDestroy();
    }

    void TestDivideOneByZero()
    {
        double one = 1.0;
        double zero = 0.0;
        double ans;
#ifdef TEST_FOR_FPE
// If we are testing for divide-by-zero, then this will throw an exception
        //TS_ASSERT_THROWS_ANYTHING(ans = one / zero);
        ans = zero*one;//otherwise compiler would complain
        TS_ASSERT_EQUALS(ans, zero);
        ans=ans*zero;//otherwise compiler would complain
#else
// If we aren't testing for it, then there will be no exception
        TS_ASSERT_THROWS_NOTHING(ans = one / zero);
        double negative_infinity=std::numeric_limits<double>::infinity();
        TS_ASSERT_EQUALS(ans, negative_infinity);
#endif
    }

    void TestDivideZeroByZero()
    {
        double zero = 0.0;
        double ans;
#ifdef TEST_FOR_FPE
// If we are testing for divide-by-zero, then this will throw an exception
        //TS_ASSERT_THROWS_ANYTHING(ans = zero / zero);
        ans = zero;//otherwise compiler would complain
        TS_ASSERT_EQUALS(ans, zero);
        ans=ans*zero;//otherwise compiler would complain
#else
// If we aren't testing for it, then there will be no exception
        TS_ASSERT_THROWS_NOTHING(ans = zero / zero);
        TS_ASSERT(std::isnan(ans));
#endif
    }
};

#endif // _TESTPETSCSETUP_HPP_
