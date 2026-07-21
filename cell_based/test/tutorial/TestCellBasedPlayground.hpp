/*

Copyright (c) 2005-2026, University of Oxford.
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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTCELLBASEDPLAYGROUND_HPP_
#define TESTCELLBASEDPLAYGROUND_HPP_

/*
 * This is a blank playground for users to experiment with the cell-based functionality
 * offered by Chaste with a minimum of effort. The user must define void Test...()
 * member methods  * of the TestPlayground class, which will be run and will compile
 * out of the box, without any need to play with CMake etc.
 *
 * We include a healthy selection of header files, with documentation, in the below
 * setup file.
 */


#include "CellBasedIncludesAndDocs.hpp"
#include "AbstractCellBasedTestSuite.hpp"


// There is no need for an int main() {} function as we use the CXX test framework


class TestPlayground : public AbstractCellBasedTestSuite
{
public:
    static void TestFirstDummyExample()
    {
        /* Insert code you want to be executed here
         * you may define multiple independent void TestName() {} functions
         */

        std::cout << "Hello World!" << std::endl;
    }

    static void DoNotTest()
    {
        /* This is a function which one may call from within the other tests, but
         * it will not be directly executed itself
         */

        std::cout << "You can't see this unless you call me!" << std::endl;
    }
};

#endif // TESTCELLBASEDPLAYGROUND_HPP_
