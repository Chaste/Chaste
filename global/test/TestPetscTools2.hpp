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

#ifndef TESTPETSCTOOLS2_HPP_
#define TESTPETSCTOOLS2_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscTools.hpp"

/**
 * Testing some methods is kind of tricky, since we really want to check
 * if they work when PETSc isn't set up as well.  Hence this file.
 *
 * \todo #2228 - should we include FakePetscSetup.hpp here?
 */
class TestPetscTools2 : public CxxTest::TestSuite
{
public:
    void TestSequentialBehaviour()
    {
        TS_ASSERT(!PetscTools::IsInitialised());
        TS_ASSERT(PetscTools::IsSequential());
        TS_ASSERT_EQUALS(PetscTools::GetMyRank(), 0u);
        TS_ASSERT_EQUALS(PetscTools::GetNumProcs(), 1u);
        TS_ASSERT(PetscTools::AmMaster());

        // This should be a no-op.  We can't check that, but just make sure it runs.
        PetscTools::Barrier();
    }
};

#endif /*TESTPETSCTOOLS2_HPP_*/
