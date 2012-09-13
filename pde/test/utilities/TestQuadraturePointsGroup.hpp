/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTQUADRATUREPOINTSGROUP_HPP_
#define TESTQUADRATUREPOINTSGROUP_HPP_

#include <cxxtest/TestSuite.h>

#include "QuadraturePointsGroup.hpp"
#include "TrianglesMeshReader.hpp"

class TestQuadraturePointsGroup : public CxxTest::TestSuite
{
public:

    void TestGetQuadPointLocations1d() throw(Exception)
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        GaussianQuadratureRule<1> quad_rule(2, 3);

        QuadraturePointsGroup<1> group(mesh,quad_rule);

        assert(quad_rule.GetNumQuadPoints()==2);
        TS_ASSERT_EQUALS(group.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(group.GetNumQuadPointsPerElement(), 2u);

        ChastePoint<1> local_quad_point_0 = quad_rule.rGetQuadPoint(0);
        ChastePoint<1> local_quad_point_1 = quad_rule.rGetQuadPoint(1);

        c_vector<double,2> X = group.Get(0,0);
        TS_ASSERT_DELTA(X(0), local_quad_point_0[0]/10, 1e-9);

        X = group.Get(0,1);
        TS_ASSERT_DELTA(X(0), local_quad_point_1[0]/10, 1e-9);

        X = group.Get(8,0);
        TS_ASSERT_DELTA(X(0), 0.8+local_quad_point_0[0]/10, 1e-9);

        X = group.Get(8,1);
        TS_ASSERT_DELTA(X(0), 0.8+local_quad_point_1[0]/10, 1e-9);

        X = group.Get(17);
        TS_ASSERT_DELTA(X(0), 0.8+local_quad_point_1[0]/10, 1e-9);
    }

    void TestGetQuadPointLocations2d() throw(Exception)
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        GaussianQuadratureRule<2> quad_rule(2,2);

        QuadraturePointsGroup<2> group(mesh,quad_rule);

        assert(quad_rule.GetNumQuadPoints()==4);
        TS_ASSERT_EQUALS(group.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(group.GetNumQuadPointsPerElement(), 4u);
        for (unsigned quad_index=0; quad_index<4; quad_index++)
        {
            c_vector<double,2> X = group.Get(0, quad_index);
            TS_ASSERT_LESS_THAN(X(0)+X(1), 1.0); // quad point in elem 0, so x+y<1

            X = group.Get(1, quad_index);
            TS_ASSERT_LESS_THAN(1.0, X(0)+X(1)); // quad point in elem 0, so x+y>1
        }

        TS_ASSERT_EQUALS(group.Size(), 8u);
        for (unsigned index=0; index<group.Size(); index++)
        {
            TS_ASSERT_LESS_THAN(group.Get(index)[0], 0.8);
            TS_ASSERT_LESS_THAN(0.2, group.Get(index)[0]);
        }
    }
};

#endif /*TESTQUADRATUREPOINTSGROUP_HPP_*/
