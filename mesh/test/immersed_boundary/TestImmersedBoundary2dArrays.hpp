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

#ifndef TESTIMMERSEDBOUNDARY2DARRAYS_HPP_
#define TESTIMMERSEDBOUNDARY2DARRAYS_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundary2dArrays.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundary2dArrays : public CxxTest::TestSuite
{
public:

    void TestMethods()
    {
        // Create an ImmersedBoundary2dArrays object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();
        ImmersedBoundary2dArrays<2> arrays(p_mesh, 0.123, 0.246, true);

        // The default array size is set in the ImmersedBoundaryMesh class to 128
        TS_ASSERT_EQUALS(p_mesh->GetNumGridPtsX(), 128u);
        TS_ASSERT_EQUALS(p_mesh->GetNumGridPtsY(), 128u);

        // Test member variables are correctly initialised
        TS_ASSERT_EQUALS(arrays.HasActiveSources(), true);
        TS_ASSERT_EQUALS(arrays.GetMesh()->GetNumNodes(), 583u);
        TS_ASSERT_EQUALS(arrays.GetMesh()->GetNumElements(), 5u);

        // Test that each multi_array member has the correct shape
        TS_ASSERT_EQUALS(arrays.rGetModifiableForceGrids().shape()[0], 2u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableForceGrids().shape()[1], 128u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableForceGrids().shape()[2], 128u);

        TS_ASSERT_EQUALS(arrays.rGetModifiableRightHandSideGrids().shape()[0], 3u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableRightHandSideGrids().shape()[1], 128u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableRightHandSideGrids().shape()[2], 128u);

        TS_ASSERT_EQUALS(arrays.rGetOperator1().shape()[0], 128u);
        TS_ASSERT_EQUALS(arrays.rGetOperator1().shape()[1], 65u);
        TS_ASSERT_EQUALS(arrays.rGetOperator2().shape()[0], 128u);
        TS_ASSERT_EQUALS(arrays.rGetOperator2().shape()[1], 65u);

        TS_ASSERT_EQUALS(arrays.rGetModifiableFourierGrids().shape()[0], 3u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableFourierGrids().shape()[1], 128u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableFourierGrids().shape()[2], 65u);
        TS_ASSERT_EQUALS(arrays.rGetModifiablePressureGrid().shape()[0], 128u);
        TS_ASSERT_EQUALS(arrays.rGetModifiablePressureGrid().shape()[1], 65u);

        TS_ASSERT_EQUALS(arrays.rGetSin2x().size(), 128u);
        TS_ASSERT_EQUALS(arrays.rGetSin2y().size(), 65u);

        // Expected values of sin(x) and sin(2x), used below
        std::vector<double> sin_x;
        std::vector<double> sin_2x;
        for (unsigned i = 0; i < 128; ++i)
        {
            sin_x.push_back(sin(M_PI * (double) i / 128));
            sin_2x.push_back(sin(2 * M_PI * (double) i / 128));
        }

        // Test that mSin2x and mSin2y are correctly initialised
        for (unsigned i = 0; i < 128; ++i)
        {
            TS_ASSERT_DELTA(arrays.rGetSin2x()[i], sin_2x[i], 1e-4);
        }
        for (unsigned j = 0; j < 65; ++j)
        {
            TS_ASSERT_DELTA(arrays.rGetSin2y()[j], sin_2x[j], 1e-4);
        }

        for (unsigned i = 0; i < 128; ++i)
        {
            for (unsigned j = 0; j < 65; ++j)
            {
                TS_ASSERT_DELTA(arrays.rGetOperator1()[i][j], 0.5*128*128*(sin_2x[i]*sin_2x[i] + sin_2x[j]*sin_2x[j]), 1e-4);
                TS_ASSERT_DELTA(arrays.rGetOperator2()[i][j], 1.0 + 2.0*128*128*(sin_x[i]*sin_x[i] + sin_x[j]*sin_x[j]), 1e-4);
            }
        }

        // Test that these methods can be used to modify the respective members
        arrays.rGetModifiableForceGrids().resize(extents[23][1][1]);
        TS_ASSERT_EQUALS(arrays.rGetModifiableForceGrids().shape()[0], 23u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableForceGrids().shape()[1], 1u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableForceGrids().shape()[2], 1u);

        arrays.rGetModifiableRightHandSideGrids().resize(extents[17][1][1]);
        TS_ASSERT_EQUALS(arrays.rGetModifiableRightHandSideGrids().shape()[0], 17u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableRightHandSideGrids().shape()[1], 1u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableRightHandSideGrids().shape()[2], 1u);

        arrays.rGetModifiableFourierGrids().resize(extents[13][2][1]);
        TS_ASSERT_EQUALS(arrays.rGetModifiableFourierGrids().shape()[0], 13u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableFourierGrids().shape()[1], 2u);
        TS_ASSERT_EQUALS(arrays.rGetModifiableFourierGrids().shape()[2], 1u);

        arrays.rGetModifiablePressureGrid().resize(extents[7][4]);
        TS_ASSERT_EQUALS(arrays.rGetModifiablePressureGrid().shape()[0], 7u);
        TS_ASSERT_EQUALS(arrays.rGetModifiablePressureGrid().shape()[1], 4u);
    }
};

#endif /*TESTIMMERSEDBOUNDARY2DARRAYS_HPP_*/