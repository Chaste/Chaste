/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYMESH_HPP_
#define TESTIMMERSEDBOUNDARYMESH_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryMesh : public CxxTest::TestSuite
{
public:
    void TestSolveNodeAndElementMapping() throw(Exception)
    {
    }

    void TestClear() throw(Exception)
    {
    }

    void TestSetupFluidVelocityGrids() throw(Exception)
    {
    }

    void TestArchiving() throw(Exception)
    {
    }

    void TestElementIterator() throw(Exception)
    {
    }

    void TestSetAndGetMethods() throw(Exception)
    {
    }

    void TestGetVectorFromAtoB() throw(Exception)
    {
        // Create a small mesh
        ImmersedBoundaryHoneycombMeshGenerator gen(1, 1, 3, 0.1, 0.3);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

        // Two helper vectors
        c_vector<double, 2> x_unit = unit_vector<double>(2, 0);
        c_vector<double, 2> y_unit = unit_vector<double>(2, 1);

        c_vector<double, 2> point_a;
        c_vector<double, 2> point_b;
        c_vector<double, 2> vec_a2b;

        // Immersed boundary meshes are always doubly-periodic on the square [0, 1] x [0, 1]

        // Regular cases where periodicity plays no part
        point_a = 0.0 * x_unit;
        point_b = 0.3 * x_unit + 0.4 * y_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], point_b[0], 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], point_b[1], 1e-6);
        TS_ASSERT_DELTA(norm_2(vec_a2b), 0.5, 1e-6);

        // x-periodicity
        point_a = 0.2 * x_unit;
        point_b = 0.8 * x_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], -0.4, 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], 0.0, 1e-6);

        // y-periodicity
        point_a = 0.2 * y_unit;
        point_b = 0.8 * y_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], -0.4, 1e-6);

        // x and y-periodicity
        point_a = 0.1 * x_unit + 0.1 * y_unit;
        point_b = 0.9 * x_unit + 0.9 * y_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], -0.2, 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], -0.2, 1e-6);
    }

    void TestGetSkewnessOfElementMassDistributionAboutAxis() throw(Exception)
    {
        // A square should have no skewness about any axis
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            for (unsigned i = 0; i < 16; i++)
            {
                double theta = 2.0 * M_PI * (double)i / 17.0;

                c_vector<double, 2> axis;
                axis[0] = cos(theta);
                axis[1] = sin(theta);

                TS_ASSERT_DELTA(p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(0, axis), 0.0, 1e-12);
            }

            delete (p_mesh);
        }

        // A triangle should have skewness
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            c_vector<double, 2> axis;
            axis[0] = 0.0;
            axis[1] = 1.0;

            double hand_calculated_skewness = -0.5656854249;

            // Test that the skewness is equal to the hand calculated value
            TS_ASSERT_DELTA(p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(0, axis) - hand_calculated_skewness, 0.0, 1e-9);

            // If we flip the axis, the skewness should be minus what it was before
            axis[1] = -1.0;
            TS_ASSERT_DELTA(p_mesh->GetSkewnessOfElementMassDistributionAboutAxis(0, axis) + hand_calculated_skewness, 0.0, 1e-9);
        }
    }

    void TestReMeshElement() throw(Exception)
    {
        // ReMeshElement where nothing should change, no regions involved
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            unsigned num_nodes = p_mesh->GetElement(0)->GetNumNodes();

            // Get locations before ReMesh
            std::vector<double> old_pos;
            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                old_pos.push_back(norm_2(p_mesh->GetElement(0)->GetNode(node_idx)->rGetLocation()));
            }

            p_mesh->ReMesh();

            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                TS_ASSERT_DELTA(old_pos[node_idx], norm_2(p_mesh->GetElement(0)->GetNode(node_idx)->rGetLocation()), 1e-12);
            }

            delete p_mesh;
        }

        // ReMeshElement where nothing should change, with regions involved
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            unsigned num_nodes = p_mesh->GetElement(0)->GetNumNodes();

            p_mesh->GetElement(0)->GetNode(0)->SetRegion(1);
            p_mesh->GetElement(0)->GetNode(1)->SetRegion(1);
            p_mesh->GetElement(0)->GetNode(2)->SetRegion(2);
            p_mesh->GetElement(0)->GetNode(3)->SetRegion(2);

            // Get locations before ReMesh
            std::vector<double> old_pos;
            std::vector<unsigned> old_regions;
            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                old_pos.push_back(norm_2(p_mesh->GetElement(0)->GetNode(node_idx)->rGetLocation()));
                old_regions.push_back(p_mesh->GetElement(0)->GetNode(node_idx)->GetRegion());
            }

            p_mesh->ReMesh();

            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                TS_ASSERT_DELTA(old_pos[node_idx], norm_2(p_mesh->GetElement(0)->GetNode(node_idx)->rGetLocation()), 1e-12);
                TS_ASSERT_EQUALS(old_regions[node_idx], p_mesh->GetElement(0)->GetNode(node_idx)->GetRegion());
            }

            delete p_mesh;
        }

        // ReMeshElement where positions should change, no regions involved
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            p_mesh->ReMesh();

            double cumulative_dist = 0.2 + sqrt(0.02);
            double node_spacing = cumulative_dist / 3.0;
            double epsilon = node_spacing - 0.1;

            // First node should not move
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(0)->rGetLocation()[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(0)->rGetLocation()[1], 0.0, 1e-12);

            // Second node up back edge slightly
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(1)->rGetLocation()[0], 0.1, 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(1)->rGetLocation()[1], epsilon, 1e-12);

            // Third node down diagonal towards (0,0)
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(2)->rGetLocation()[0], 0.1 - 0.2 * epsilon / sqrt(0.02), 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(2)->rGetLocation()[1], 0.1 - 0.2 * epsilon / sqrt(0.02), 1e-12);

            delete p_mesh;
        }

        // ReMeshElement where positions should change, with regions involved
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));
            nodes.push_back(new Node<2>(4, true, 0.0, 0.03));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2, 2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            p_mesh->GetElement(0)->GetNode(0)->SetRegion(1);
            p_mesh->GetElement(0)->GetNode(1)->SetRegion(1);
            p_mesh->GetElement(0)->GetNode(2)->SetRegion(2);
            p_mesh->GetElement(0)->GetNode(3)->SetRegion(2);
            p_mesh->GetElement(0)->GetNode(4)->SetRegion(1);

            p_mesh->ReMesh();

            // Third node should not move
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(2)->rGetLocation()[0], 0.1, 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(2)->rGetLocation()[1], 0.1, 1e-12);

            // Fourth moves towards third
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(3)->rGetLocation()[0], 0.02, 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(3)->rGetLocation()[1], 0.1, 1e-12);

            // Fifth moves towards fourth
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(4)->rGetLocation()[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(4)->rGetLocation()[1], 0.04, 1e-12);

            // First moves towards second
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(0)->rGetLocation()[0], 0.04, 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(0)->rGetLocation()[1], 0.0, 1e-12);

            // Second moves towards third
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(1)->rGetLocation()[0], 0.1, 1e-12);
            TS_ASSERT_DELTA(p_mesh->GetElement(0)->GetNode(1)->rGetLocation()[1], 0.02, 1e-12);

            // Due to reposition of nodes, node 4 should now be part of region 2
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNode(0)->GetRegion(), 1u)
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNode(1)->GetRegion(), 1u)
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNode(2)->GetRegion(), 2u)
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNode(3)->GetRegion(), 2u)
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNode(4)->GetRegion(), 2u)

            delete p_mesh;
        }
    }
};

#endif /*TESTIMMERSEDBOUNDARYMESH_HPP_*/
