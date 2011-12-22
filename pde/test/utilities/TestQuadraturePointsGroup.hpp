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

        GaussianQuadratureRule<1> quad_rule(2);

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

        GaussianQuadratureRule<2> quad_rule(2);

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
