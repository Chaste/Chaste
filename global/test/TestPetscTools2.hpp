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

#ifndef TESTPETSCTOOLS2_HPP_
#define TESTPETSCTOOLS2_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscTools.hpp"

/**
 * Testing some methods is kind of tricky, since we really want to check
 * if they work when PETSc isn't set up as well.  Hence this file.
 */
class TestPetscTools2 : public CxxTest::TestSuite
{
public:
    void TestSequentialBehaviour()
    {
        TS_ASSERT(PetscTools::IsSequential());
        TS_ASSERT_EQUALS(PetscTools::GetMyRank(), 0u);
        TS_ASSERT_EQUALS(PetscTools::GetNumProcs(), 1u);
        TS_ASSERT(PetscTools::AmMaster());

        // This should be a no-op.  We can't check that, but just make sure it runs.
        PetscTools::Barrier();
    }
};

#endif /*TESTPETSCTOOLS2_HPP_*/
