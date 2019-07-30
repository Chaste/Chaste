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

#ifndef TESTMATHSCUSTOMFUNCTIONS_HPP_
#define TESTMATHSCUSTOMFUNCTIONS_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include "MathsCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestMathsCustomFunctions : public CxxTest::TestSuite
{
public:

    void TestSmallPow()
    {
        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_DELTA(SmallPow(0.0, i), pow(0.0, double(i)), DBL_EPSILON);
            //Integer version
            TS_ASSERT_EQUALS(SmallPow(5u, i), (unsigned)floor(0.5 + pow(5.0, double(i))));
            TS_ASSERT_DELTA(SmallPow(-1.67e3, i), pow(-1.67e3, double(i)), 1e14);//(1e3)^10/DBL_EPSILON
            //Note the following should be within DBL_EPSILON on a 64-bit machine.  However,
            //32-bit floating point accumulated errors push both SmallPow and pow away from integer values.
            //There also seems to be 32-bit issue within TS_ASSERT_DELTA...
            double difference = SmallPow(75.0, i) - pow(75.0, double(i));
            TS_ASSERT_EQUALS(difference, 0.0);
        }
    }

    void TestDivides()
    {
        TS_ASSERT_EQUALS(Divides(0.7, 0.1),  false);
        TS_ASSERT_EQUALS(Divides(0.07, 0.1),  false);
        TS_ASSERT_EQUALS(Divides(0.1, 0.1),  true);
        TS_ASSERT_EQUALS(Divides(1e10, 1e10),  true);
        TS_ASSERT_EQUALS(Divides(5.7e10, 5.7e20),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 0.1),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 1.0),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 10.0),  true);
        TS_ASSERT_EQUALS(Divides(0.01, 100.01),  true);

        // Note that Divides() returns false if you attempt to divide zero by a non-zero number
        TS_ASSERT_EQUALS(Divides(0.01, 0.00),  false);
    }

    void TestCeilDivide()
    {
        TS_ASSERT_EQUALS(CeilDivide(28u, 7u), 4u);
        TS_ASSERT_EQUALS(CeilDivide(27u, 7u), 4u);
        TS_ASSERT_EQUALS(CeilDivide(29u, 7u), 5u);
        TS_ASSERT_EQUALS(CeilDivide(6u, 3u), 2u);
        TS_ASSERT_EQUALS(CeilDivide(5u, 3u), 2u);
        TS_ASSERT_EQUALS(CeilDivide(0u, 3u), 0u);
        TS_ASSERT_EQUALS(CeilDivide(2342368901u, UINT_MAX), 1u);
    }

    void TestSignum()
    {
        TS_ASSERT_EQUALS(Signum(0.0), 0.0);
        TS_ASSERT_EQUALS(Signum(0.01), 1.0);
        TS_ASSERT_EQUALS(Signum(-0.01), -1.0);
        TS_ASSERT_EQUALS(Signum(DBL_MAX), 1.0);
        TS_ASSERT_EQUALS(Signum(-DBL_MAX), -1.0);
        TS_ASSERT_EQUALS(Signum(DBL_EPSILON), 1.0);
    }

    void TestCompareDoubles()
    {
        TS_ASSERT(CompareDoubles::IsNearZero(DBL_EPSILON, 2*DBL_EPSILON));
        TS_ASSERT(CompareDoubles::IsNearZero(-0.2, 0.200001));
        TS_ASSERT(CompareDoubles::IsNearZero(0.1, 0.1));
        TS_ASSERT(CompareDoubles::IsNearZero(-0.1, 0.1));
        TS_ASSERT(!CompareDoubles::IsNearZero(1/3.0, 0.33333));
        TS_ASSERT(!CompareDoubles::IsNearZero(-1/3.0, 0.33333));

        TS_ASSERT(CompareDoubles::WithinRelativeTolerance(10, 10.1, 0.01)); // 1%
        TS_ASSERT(CompareDoubles::WithinRelativeTolerance(10.1, 10, 0.01)); // 1%
        TS_ASSERT(!CompareDoubles::WithinRelativeTolerance(10.1000000001, 10, 0.01)); // 1%
        TS_ASSERT(!CompareDoubles::WithinRelativeTolerance(10, 0.9999999999, 0.01)); // 1%
        TS_ASSERT(CompareDoubles::WithinAbsoluteTolerance(99.99999, 100, 0.000011));
        TS_ASSERT(!CompareDoubles::WithinAbsoluteTolerance(100, 99.9999, 0.0001)); // Equality is interesting...
        TS_ASSERT(CompareDoubles::WithinAbsoluteTolerance(2.0, 3.0, 1.0));
        TS_ASSERT(!CompareDoubles::WithinAbsoluteTolerance(99.99999, 100, 0.000009));

        TS_ASSERT(CompareDoubles::WithinAnyTolerance(10, 10.1, 0.01, 0.001));
        TS_ASSERT(CompareDoubles::WithinAnyTolerance(10, 10.1, 0.01)); // Relative only
        TS_ASSERT(!CompareDoubles::WithinAnyTolerance(10, 10.1, 0.001, 0.01));
        TS_ASSERT(CompareDoubles::WithinAnyTolerance(10, 10.1, 0.001, 0.2));
        TS_ASSERT(CompareDoubles::WithinAnyTolerance(10, 10)); // Exact equality...
        TS_ASSERT(!CompareDoubles::WithinAnyTolerance(10, 10.1, 0.001, 0.01, true)); // Print error

        TS_ASSERT(CompareDoubles::WithinTolerance(0.001, 0.002, 0.0015, true)); // Absolute tol
        TS_ASSERT(!CompareDoubles::WithinTolerance(0.001, 0.002, 0.0015, false)); // Relative tol
        TS_ASSERT(!CompareDoubles::WithinTolerance(0.001, 0.002, 0.0005, true)); // Absolute tol

        TS_ASSERT_DELTA(CompareDoubles::Difference(0.001, 0.002, true), 0.001, 1e-12); // Absolute diff
        TS_ASSERT_DELTA(CompareDoubles::Difference(0.001, 0.002, false), 1.0, 1e-12); // Max. relative diff
    }

    void TestSafeDivide()
    {
        TS_ASSERT_EQUALS(SafeDivide(DBL_MAX, 0.5), DBL_MAX);
        TS_ASSERT_DELTA(SafeDivide(0.0, 1.0), 0.0, 1e-6);
        TS_ASSERT_DELTA(SafeDivide(1.0, 2.0), 0.5, 1e-6);
    }
};

#endif /*TESTMATHSCUSTOMFUNCTIONS_HPP_*/
