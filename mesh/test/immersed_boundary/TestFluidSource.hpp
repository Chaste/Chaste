/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTFLUIDSOURCE_HPP_
#define TESTFLUIDSOURCE_HPP_

#include <cxxtest/TestSuite.h>

#include "FluidSource.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestFluidSource : public CxxTest::TestSuite
{
public:

    void TestFluidSourceConstructors()
    {
        ChastePoint<2> point(0.1, -5.0);
        FluidSource<2> source_by_point(1, point);
        TS_ASSERT_EQUALS(source_by_point.GetIndex(), 1u);
        TS_ASSERT_DELTA(source_by_point.rGetLocation()[0], 0.1, 1e-6);
        TS_ASSERT_DELTA(source_by_point.rGetLocation()[1], -5.0, 1e-6);

        c_vector<double, 3> vector;
        vector[0] = 3.2;
        vector[1] = 0.1;
        vector[2] = -5.5;
        FluidSource<3> source_by_vector(0, vector);
        TS_ASSERT_EQUALS(source_by_vector.GetIndex(), 0u);
        TS_ASSERT_DELTA(source_by_vector.rGetLocation()[0], 3.2, 1e-6);
        TS_ASSERT_DELTA(source_by_vector.rGetLocation()[1], 0.1, 1e-6);
        TS_ASSERT_DELTA(source_by_vector.rGetLocation()[2], -5.5, 1e-6);

        FluidSource<2> source_2d(5, 1.3, 6.7);
        TS_ASSERT_EQUALS(source_2d.GetIndex(), 5u);
        TS_ASSERT_DELTA(source_2d.rGetLocation()[0], 1.3, 1e-6);
        TS_ASSERT_DELTA(source_2d.rGetLocation()[1], 6.7, 1e-6);

        FluidSource<3> source_3d(8, 8.0, 9.7, 7.1);
        TS_ASSERT_EQUALS(source_3d.GetIndex(), 8u);
        TS_ASSERT_DELTA(source_3d.rGetLocation()[0], 8.0, 1e-6);
        TS_ASSERT_DELTA(source_3d.rGetLocation()[1], 9.7, 1e-6);
        TS_ASSERT_DELTA(source_3d.rGetLocation()[2], 7.1, 1e-6);
    }

    void TestFluidSourceIndexMethods()
    {
        FluidSource<2> source(0, 1.3, 6.7);
        TS_ASSERT_EQUALS(source.GetIndex(), 0u);

        source.SetIndex(132);
        auto res = source.GetIndex();
        TS_ASSERT_EQUALS(res, 132u);
    }

    void TestStrengthMethods()
    {
        FluidSource<2> source(1, 0.1, -5.0);
        TS_ASSERT_DELTA(source.GetStrength(), 0.0, 1e-6);

        source.SetStrength(18.4);
        TS_ASSERT_DELTA(source.GetStrength(), 18.4, 1e-6);
    }

    void TestLocationMethods()
    {
        FluidSource<3> source(0, 8.0, 9.7, 7.1);

        ChastePoint<3> point = source.GetPoint();
        TS_ASSERT_DELTA(point[0], 8.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 9.7, 1e-6);
        TS_ASSERT_DELTA(point[2], 7.1, 1e-6);

        c_vector<double, 3> location = source.rGetLocation();
        TS_ASSERT_DELTA(location[0], 8.0, 1e-6);
        TS_ASSERT_DELTA(location[1], 9.7, 1e-6);
        TS_ASSERT_DELTA(location[2], 7.1, 1e-6);

        source.rGetModifiableLocation()[1] += 100.0;
        c_vector<double, 3> new_location = source.rGetLocation();
        TS_ASSERT_DELTA(new_location[0], 8.0, 1e-6);
        TS_ASSERT_DELTA(new_location[1], 109.7, 1e-6);
        TS_ASSERT_DELTA(new_location[2], 7.1, 1e-6);
    }

    void TestElementMethods()
    {
        FluidSource<2> source(0, 1.3, 6.7);

        TS_ASSERT_EQUALS(source.IsSourceAssociatedWithElement(), false);
        TS_ASSERT_EQUALS(source.GetAssociatedElementIndex(), UINT_MAX);

        source.SetIfSourceIsAssociatedWithElement(true);
        TS_ASSERT_EQUALS(source.IsSourceAssociatedWithElement(), true);
        TS_ASSERT_EQUALS(source.GetAssociatedElementIndex(), UINT_MAX);

        source.SetAssociatedElementIndex(15);
        TS_ASSERT_EQUALS(source.IsSourceAssociatedWithElement(), true);
        TS_ASSERT_EQUALS(source.GetAssociatedElementIndex(), 15u);
    }
};

#endif /*TESTFLUIDSOURCE_HPP_*/