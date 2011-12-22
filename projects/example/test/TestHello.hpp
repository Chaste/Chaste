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

#ifndef TESTHELLO_HPP_
#define TESTHELLO_HPP_

#include <cxxtest/TestSuite.h>

#include "Hello.hpp"

/**
 * @file
 *
 * This is an example of a CxxTest test suite, used to test the source
 * code, and also used to run simulations (as it provides a handy
 * shortcut to compile and link against the correct libraries using scons).
 *
 * You can #include any of the files in the project 'src' folder.
 * For example here we #include "Hello.hpp"
 *
 * You can utilise any of the code in the main the Chaste trunk
 * in exactly the same way.
 * NOTE: you will have to alter the project SConscript file lines 41-44
 * to enable #including of code from the 'heart', 'cell_based' or 'crypt'
 * components of Chaste.
 */

class TestHello : public CxxTest::TestSuite
{
public:
    void testHello()
    {
        // Create an object called 'world' of class 'Hello',
        // (Hello.hpp is #included from the 'src' folder.)
        Hello world("Hello world!");

        // The TS_ASSERT macros are used to test that the object performs as expected
        TS_ASSERT_EQUALS(world.GetMessage(), "Hello world!");
        TS_ASSERT_THROWS_THIS(world.Complain("I don't like you"),
                              "I don't like you");
    }
};

#endif /*TESTHELLO_HPP_*/
