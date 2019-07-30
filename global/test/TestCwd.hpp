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

#ifndef _TESTCWD_HPP_
#define _TESTCWD_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdlib>
#include <climits>

#include "FakePetscSetup.hpp"
#include "ChasteBuildRoot.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "BoostFilesystem.hpp"
#include "IsNan.hpp"

/**
 * Test for a strange 'feature' of Debian sarge systems, where the
 * current working directory changes on PETSc initialisation.  There
 * is now code in PetscSetupUtils::CommonSetup() to work around this.
 */
class TestCwd : public CxxTest::TestSuite
{
public:
    void TestShowCwd()
    {
        std::string chaste_source_root(ChasteSourceRootDir());
        fs::path chaste_source_root_path(chaste_source_root);
        fs::path cwd(GetCurrentWorkingDirectory()+"/");
        TS_ASSERT(fs::equivalent(cwd,chaste_source_root));
    }

    void TestDivideOneByZero()
    {
        double one = 1.0;
        double zero = 0.0;
        double ans;
#ifdef TEST_FOR_FPE
// If we are testing for divide-by-zero, then this will throw an exception
        //TS_ASSERT_THROWS_ANYTHING(ans = one / zero);
        ans=zero*one;//otherwise compiler would complain
        ans=ans*zero;//otherwise compiler would complain
#else
// If we aren't testing for it, then there will be no exception
        TS_ASSERT_THROWS_NOTHING(ans = one / zero);
        double negative_infinity = std::numeric_limits<double>::infinity();
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
        ans=zero;//otherwise compiler would complain
        ans=ans*zero;//otherwise compiler would complain
#else
// If we aren't testing for it, then there will be no exception
        TS_ASSERT_THROWS_NOTHING(ans = zero / zero);
        TS_ASSERT(std::isnan(ans));
#endif
    }
};

#endif /*_TESTCWD_HPP_*/
