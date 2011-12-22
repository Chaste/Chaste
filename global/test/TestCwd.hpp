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

#ifndef _TESTCWD_HPP_
#define _TESTCWD_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdlib>

#include "PetscSetupAndFinalize.hpp"
#include "ChasteBuildRoot.hpp"
#include "GetCurrentWorkingDirectory.hpp"

/**
 * Test for a strange 'feature' of Debian sarge systems, where the
 * current working directory changes on PETSc initialisation.  There
 * is now code in PetscSetupAndFinalize.hpp to work around this.
 */
class TestCwd : public CxxTest::TestSuite
{
public:
    void TestShowCwd()
    {
        TS_ASSERT_EQUALS(system("pwd"), 0);
        TS_ASSERT_EQUALS(system("ls -l io/test/data"), 0);
        std::string chaste_build_root(ChasteBuildRootDir());
        TS_ASSERT_EQUALS(GetCurrentWorkingDirectory() + "/", chaste_build_root);
    }
};

#endif /*_TESTCWD_HPP_*/
