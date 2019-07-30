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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTWRITINGTESTTUTORIAL_HPP_
#define TESTWRITINGTESTTUTORIAL_HPP_

/*
 * = Writing tests =
 * We do not use `int main()` methods in Chaste. Instead, we write ''tests'', which are run using !CxxTest.
 * Tests are used both as:
 * (i) part of the testing environment - every class in the source code has an equivalent test file which tests each aspect of its functionality, making use of the `TS_ASSERT`s as described below; and
 * (ii) for experimental work, which involve writing a 'test' as below but generally without `TS_ASSERT`s.
 *
 * This tutorial shows how to write a test using !CxxTest. Note that the full code is given at the bottom of the page.
 */

/* First, the following header file needs to be included.*/
#include <cxxtest/TestSuite.h>

/*
 * Now we have to define a class containing the tests. It is sensible to name the class with the same name as the file name. The class should inherit from {{{CxxTest::TestSuite}}}.
 */
class TestWritingTestsTutorial: public CxxTest::TestSuite
{
/*
 * Now we define some tests, which must be '''public''', begin with the word 'Test', return {{{void}}}, and take in no parameters.
 */
public:
    void TestOnePlusOneEqualsTwo()
    {
/*
 * To test whether two integers are equal, we can use the macro {{{TS_ASSERT_EQUALS}}}.
 */
        int some_number = 1 + 1;
        TS_ASSERT_EQUALS(some_number, 2);
/*
 * To test whether two numbers are equal to within a certain (absolute) tolerance we can use {{{TS_ASSERT_DELTA}}}.
 * This should almost always be used when comparing two {{{double}}}s.  (See also class:CompareDoubles for more
 * advanced comparisons.)
 */
        double another_number = 1.000001 + 1.0001;
        TS_ASSERT_DELTA(another_number, 2.0, 1e-2);
    }
/*
 * This second test shows some of the other {{{TS_ASSERT}}} macros that are available.
 * The {{{}} part of the signature is there to make sure that full details of any
 * uncaught exceptions are reported.
 */
    void TestSomeOtherStuff()
    {
        TS_ASSERT(1==1); // however, it is better to use TS_ASSERT_EQUALS, below
        TS_ASSERT_EQUALS((true||false), true);
        TS_ASSERT_DIFFERS(1.348329534564385643543957436, 1.348329534564395643543957436);
        TS_ASSERT_LESS_THAN(2.71828183, 3.14159265); // Note: to test if x is greater than y, use TS_ASSERT_LESS_THAN(y,x)
        TS_ASSERT_LESS_THAN_EQUALS(-1e100, 1e100);
        TS_ASSERT_THROWS_ANYTHING(throw 0;); // normally you would put a function call inside the brackets

        unsigned x;
        // The following TS_ASSERT_THROWS_NOTHING may be useful if you want to be certain that there are no uncaught exceptions
        TS_ASSERT_THROWS_NOTHING(x=1;);  // normally you would put a function call inside the brackets
        TS_ASSERT_EQUALS(x, 1u); //Note that x and 1u are of the same type: unsigned integer
    }

/* Other useful macros include `TS_ASSERT_THROWS_THIS` and `TS_ASSERT_THROWS_CONTAINS` for testing exception
 * messages.
 *
 * Note that methods that don't start with 'Test' are compiled but not run. So, if you want to stop a single
 * test running, just put an 'x' or a 'donot' (for instance) before its name.
 */
    void donotTestThis()
    {
        TS_ASSERT_EQUALS(1,   2);
        TS_ASSERT_EQUALS(1u,  2u);
        TS_ASSERT_EQUALS(1.0, 2.0);
    }
};
/*
 *
 * To run this code, first copy it into a file, say, called `TestWritingTests.hpp` in the directory `global/test/`.
 * Second, add the full name of your new file to the relevant continuous test pack, say `[path/to/Chaste]/global/test/ContinuousTestPack.txt`.
 * Third, from the command line, run
{{{
#!sh
cd [path/to/ChasteBuild]
ccmake [path/to/Chaste]
}}}
 * Then press `c` to configure, `e` to exit, and `g` to generate. Finally, run
{{{
#!sh
make global
ctest -V -R TestWritingTests
}}}
 */
#endif /*TESTWRITINGTESTTUTORIAL_HPP_*/
