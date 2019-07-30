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

#ifndef TESTQUADRATUREPOINTSGROUP_HPP_
#define TESTQUADRATUREPOINTSGROUP_HPP_

#include <cxxtest/TestSuite.h>

#include "QuadraturePointsGroup.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestQuadraturePointsGroup : public CxxTest::TestSuite
{
public:

    void TestGetQuadPointLocations1d()
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        GaussianQuadratureRule<1> quad_rule(3);

        QuadraturePointsGroup<1> group(mesh,quad_rule);

        TS_ASSERT_EQUALS(quad_rule.GetNumQuadPoints(), 2u);
        TS_ASSERT_EQUALS(group.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(group.GetNumQuadPointsPerElement(), 2u);

        ChastePoint<1> local_quad_point_0 = quad_rule.rGetQuadPoint(0);
        ChastePoint<1> local_quad_point_1 = quad_rule.rGetQuadPoint(1);

        c_vector<double,2> quad_point = group.rGet(0, 0);
        TS_ASSERT_DELTA(quad_point(0), local_quad_point_0[0]/10, 1e-9);

        quad_point = group.rGet(0, 1);
        TS_ASSERT_DELTA(quad_point(0), local_quad_point_1[0]/10, 1e-9);

        quad_point = group.rGet(8, 0);
        TS_ASSERT_DELTA(quad_point(0), 0.8+local_quad_point_0[0]/10, 1e-9);

        quad_point = group.rGet(8, 1);
        TS_ASSERT_DELTA(quad_point(0), 0.8+local_quad_point_1[0]/10, 1e-9);

        quad_point = group.rGet(17);
        TS_ASSERT_DELTA(quad_point(0), 0.8+local_quad_point_1[0]/10, 1e-9);
    }

    void TestGetQuadPointLocations2d()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        GaussianQuadratureRule<2> quad_rule(2);

        QuadraturePointsGroup<2> group(mesh,quad_rule);

        TS_ASSERT_EQUALS(quad_rule.GetNumQuadPoints(), 3u);
        TS_ASSERT_EQUALS(group.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(group.GetNumQuadPointsPerElement(), 3u);

        // Element 0
        for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
        {
            c_vector<double,2> quad_point = group.rGet(0, quad_index);
            TS_ASSERT_LESS_THAN(quad_point(0)+quad_point(1), 1.0); // quad point in elem 0, so x+y<1
        }
        // Element 1
        for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
        {
            c_vector<double,2> quad_point = group.rGet(1, quad_index);
            TS_ASSERT_LESS_THAN(1.0, quad_point(0)+quad_point(1)); // quad point in elem 1, so x+y>1
        }

        TS_ASSERT_EQUALS(group.Size(), 6u);
        for (unsigned index=0; index<group.Size(); index++)
        {
            TS_ASSERT_LESS_THAN_EQUALS(group.rGet(index)[0], 5.0/6.0);
            TS_ASSERT_LESS_THAN_EQUALS(1.0/6.0, group.rGet(index)[0]);
        }
    }

    void TestGetQuadPointLocations2dDistributed()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements");
        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);

        mesh.ConstructFromMeshReader(reader);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);

        GaussianQuadratureRule<2> quad_rule(2);
        TS_ASSERT_EQUALS(quad_rule.GetNumQuadPoints(), 3u);

        QuadraturePointsGroup<2> group(mesh, quad_rule);
        TS_ASSERT_EQUALS(group.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(group.GetNumQuadPointsPerElement(), 3u);

        // Element 0
        try
        {
            mesh.GetElement(0); //Throws if not owned
            for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
            {
                c_vector<double,2> quad_point = group.rGet(0, quad_index);
                TS_ASSERT_LESS_THAN(quad_point(0)+quad_point(1), 1.0); // quad point in elem 0, so x+y<1
            }
        }
        catch (Exception&)
        {
            //If not, then we know nothing about these quad points
            for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
            {
                c_vector<double,2> quad_point = group.rGet(0, quad_index);
                TS_ASSERT_EQUALS(quad_point(0), DOUBLE_UNSET);
                TS_ASSERT_EQUALS(quad_point(1), DOUBLE_UNSET);
            }

        }
        // Element 1
        try
        {
            mesh.GetElement(0); //Throws if not owned
            for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
            {
                c_vector<double,2> quad_point = group.rGet(1, quad_index);
                TS_ASSERT_LESS_THAN(1.0, quad_point(0)+quad_point(1)); // quad point in elem 1, so x+y>1
            }
        }
        catch (Exception&)
        {
            //If not, then we know nothing about these quad points
            for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
            {
                c_vector<double,2> quad_point = group.rGet(0, quad_index);
                TS_ASSERT_EQUALS(quad_point(0), DOUBLE_UNSET);
                TS_ASSERT_EQUALS(quad_point(1), DOUBLE_UNSET);
            }
        }

        TS_ASSERT_EQUALS(group.Size(), 6u);
        for (unsigned index=0; index<group.Size(); index++)
        {
            double quad_x = group.rGet(index)[0];
            if (quad_x != DOUBLE_UNSET)
            {
                TS_ASSERT_LESS_THAN_EQUALS(quad_x, 5.0/6.0);
                TS_ASSERT_LESS_THAN_EQUALS(1.0/6.0, quad_x);
            }
        }
    }
};

#endif /*TESTQUADRATUREPOINTSGROUP_HPP_*/
