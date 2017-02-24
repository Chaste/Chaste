/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef TESTMONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_
#define TESTMONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "Node.hpp"
#include "VertexElement.hpp"
#include "MonolayerVertexMeshGenerator.hpp"

#include <boost/assign/list_of.hpp>
using boost::assign::list_of;

#include "FakePetscSetup.hpp"

#include <boost/lexical_cast.hpp>

#include "Debug.hpp"

class TestMonolayerVertexMeshCustomFunctions : public CxxTest::TestSuite
{
public:
    void TestElementHasNode()
    {
        Node<3> n18(18);
        Node<3> n27(27);
        std::vector<Node<3>*> nodes;
        nodes.push_back(&n27);
        VertexElement<2, 3> f99(99, nodes);
        nodes[0] = &n18;
        VertexElement<3, 3> e66(66, nodes);

        TS_ASSERT_EQUALS(ElementHasNode(&f99, 18u), false);
        TS_ASSERT_EQUALS(ElementHasNode(&f99, 27u), true);

        TS_ASSERT_EQUALS(ElementHasNode(&e66, 18u), true);
        TS_ASSERT_EQUALS(ElementHasNode(&e66, 27u), false);
    }

    void TestGetSharedElementnFaceIndicesnBoundaryFace()
    {
        /*
         * Create a mesh comprising 3 elements, as shown below.
         *      ___ __
         *     |   |__|
         *  ___|___|
         * |   |   |
         * |___|___|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, true, 1.0, 1.0));
        nodes.push_back(new Node<3>(5, true, 2.0, 1.0));
        nodes.push_back(new Node<3>(6, true, 2.0, 1.5));
        nodes.push_back(new Node<3>(7, true, 2.5, 1.5));
        nodes.push_back(new Node<3>(8, true, 1.0, 2.0));
        nodes.push_back(new Node<3>(9, true, 2.0, 2.0));
        nodes.push_back(new Node<3>(10, true, 2.5, 2.0));
        const unsigned node_indices_elem_0[4] = { 0, 1, 4, 3 };
        const unsigned node_indices_elem_1[4] = { 1, 2, 5, 4 };
        const unsigned node_indices_elem_2[5] = { 4, 5, 6, 9, 8 };
        const unsigned node_indices_elem_3[4] = { 6, 7, 10, 9 };

        MonolayerVertexMeshGenerator builder(nodes, "TestGetSharedElementnFaceIndices");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(5, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder("MonolayerCustomFunctions", "Before");

        // Test some random node pairs
        std::set<unsigned> elems_tmp, faces_tmp;
        unsigned node_1, node_2;
        node_1 = 0;
        node_2 = 1;
        elems_tmp = list_of(0);
        faces_tmp = list_of(0)(2);
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);
        node_1 = 6;
        node_2 = 9;
        elems_tmp = list_of(2)(3);
        faces_tmp = list_of(11)(17)(14);
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);
        node_1 = 5;
        node_2 = 8;
        elems_tmp = list_of(2);
        faces_tmp = list_of(11);
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);
        node_1 = 18;
        node_2 = 21;
        elems_tmp = list_of(3);
        faces_tmp = list_of(18)(20);
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);
        node_1 = 12;
        node_2 = 13;
        elems_tmp = list_of(1);
        faces_tmp = list_of(7)(8);
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);
        node_1 = 4;
        node_2 = 16;
        elems_tmp = list_of(1)(2);
        faces_tmp = list_of(10);
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);
        node_1 = 7;
        node_2 = 15;
        elems_tmp.clear();
        faces_tmp.clear();
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);
        node_1 = 4;
        node_2 = 15;
        elems_tmp = list_of(0)(1)(2);
        faces_tmp = list_of(3)(4)(10)(16);
        TS_ASSERT_EQUALS(GetSharedElementIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), elems_tmp);
        TS_ASSERT_EQUALS(GetSharedFaceIndices(vertex_mesh.GetNode(node_1), vertex_mesh.GetNode(node_2)), faces_tmp);

        // Test for boundary faces
        const std::set<unsigned> non_boundary_indices = list_of(3)(10)(14);
        for (unsigned face_index = 0; face_index < vertex_mesh.GetNumFaces(); ++face_index)
        {
            const VertexElement<2, 3>* p_face = vertex_mesh.GetFace(face_index);
            TS_ASSERT_EQUALS(IsFaceOnBoundary(p_face), non_boundary_indices.count(face_index) == 0);
        }
    }

    void SomeForce(MutableVertexMesh<3, 3>& vertex_mesh, double dt = 0.1)
    {
        c_vector<double, 3> force = zero_vector<double>(3);

        double mVolumeParameter = 1, mTargetVolume = 0.5,
               mApicalAreaParameter = 0, mTargetApicalArea = 0, mApicalEdgeParameter = 0,
               mBasalAreaParameter = 0, mTargetBasalArea = 0, mBasalEdgeParameter = 0,
               mLateralEdgeParameter = 0;

        // Define some helper variables
        const unsigned num_nodes = vertex_mesh.GetNumNodes();
        const unsigned num_elements = vertex_mesh.GetNumElements();

        // Begin by computing the volumes of each element in the mesh, to avoid having to do this multiple times
        std::vector<double> element_volumes(num_elements);
        std::vector<double> apical_areas(num_elements);
        std::vector<double> basal_areas(num_elements);
        for (unsigned elem_index = 0; elem_index < num_elements; ++elem_index)
        {
            const VertexElement<3, 3>* p_elem = vertex_mesh.GetElement(elem_index);
            assert(elem_index == p_elem->GetIndex());
            element_volumes[elem_index] = vertex_mesh.GetVolumeOfElement(elem_index);
            apical_areas[elem_index] = vertex_mesh.CalculateAreaOfFace(GetApicalFace(p_elem));
            basal_areas[elem_index] = vertex_mesh.CalculateAreaOfFace(GetBasalFace(p_elem));
        }

        // Iterate over nodes in the cell population
        for (unsigned global_index = 0; global_index < num_nodes; ++global_index)
        {
            Node<3>* p_this_node = vertex_mesh.GetNode(global_index);

            assert(global_index == p_this_node->GetIndex());

            // Get the type of node. 1=basal; 2=apical
            const Monolayer::v_type node_type = GetNodeType(p_this_node);

            c_vector<double, 3> basal_face_contribution = zero_vector<double>(3);
            c_vector<double, 3> ab_edge_contribution = zero_vector<double>(3);
            c_vector<double, 3> apical_face_contribution = zero_vector<double>(3);
            c_vector<double, 3> lateral_edge_contribution = zero_vector<double>(3);
            c_vector<double, 3> volume_contribution = zero_vector<double>(3);

            // A variable to store such that the apical/basal edge forces are not counted twice for non-boundary edges.
            std::set<unsigned> neighbour_node_indices;

            // Find the indices of the elements owned by this node
            const std::set<unsigned> containing_elem_indices = p_this_node->rGetContainingElementIndices();

            // Iterate over these elements
            for (std::set<unsigned>::const_iterator iter = containing_elem_indices.begin();
                 iter != containing_elem_indices.end();
                 ++iter)
            {
                // Get this element, its index and its number of nodes
                VertexElement<3, 3>* p_element = vertex_mesh.GetElement(*iter);

                const unsigned elem_index = p_element->GetIndex();

                // Calculating volume contribution
                c_vector<double, 3> element_volume_gradient = vertex_mesh.GetVolumeGradientofElementAtNode(p_element, global_index);

                // Add the force contribution from this cell's volume compressibility (note the minus sign)
                volume_contribution -= mVolumeParameter * element_volume_gradient * (element_volumes[elem_index] - mTargetVolume);

                // {
                // const double k(std::abs(element_volume_gradient[0]));
                // PRINT_4_VARIABLES(element_volumes[elem_index], mTargetVolume, p_this_node->GetIndex(), k)
                // PRINT_CONTAINER(element_volume_gradient/k);
                // PRINT_CONTAINER(volume_contribution)
                // }
            }

            c_vector<double, 3> force_on_node = basal_face_contribution + ab_edge_contribution + apical_face_contribution
                + lateral_edge_contribution + volume_contribution;
            p_this_node->AddAppliedForceContribution(force_on_node);
        }

        for (unsigned global_index = 0; global_index < num_nodes; ++global_index)
        {
            Node<3>* p_this_node = vertex_mesh.GetNode(global_index);
            p_this_node->rGetModifiableLocation() += p_this_node->rGetAppliedForce() * dt;
            p_this_node->ClearAppliedForce();
        }
    }

    void TestFaceRearrangeNode()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 1.0, 1.0));

        VertexElement<2, 3>* face = new VertexElement<2, 3>(0, nodes);
        std::vector<VertexElement<2, 3>*> v_f;
        v_f.push_back(face);
        PrintElement(face);
        MutableVertexMesh<2, 3> meh(nodes, v_f);

        {
            VertexMeshWriter<2, 3> writer(OUTPUT_NAME + std::string("/FaceRearrangeNode"), "blabla");
            writer.WriteVtkUsingMesh(meh, "before");
        }

        c_vector<double, 3> vv = face->GetCentroid();
        vv[2] += 3;
        face->FaceRearrangeNodes(vv);
        PrintElement(face);
        {
            VertexMeshWriter<2, 3> writer(OUTPUT_NAME + std::string("/FaceRearrangeNode"), "blabla", false);
            writer.WriteVtkUsingMesh(meh, "after");
        }
    }

    void TestCheckingGradientDeviation()
    {
        /*
         * Create a mesh with a square element, as shown below.
         *  _____
         * |     |
         * |     |
         * |_____|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 1.0, 1.0));
        const unsigned node_indices_elem_0[4] = { 0, 1, 3, 2 };

        MonolayerVertexMeshGenerator builder(nodes, "TestGradientDeviation");
        builder.BuildElementWith(4, node_indices_elem_0);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        PrintElement(vertex_mesh.GetElement(0));

        for (unsigned i = 0; i < vertex_mesh.GetNumNodes(); ++i)
            vertex_mesh.GetNode(i)->rGetModifiableLocation()[2] += 5;

        builder.WriteVtkWithSubfolder("MonolayerCustomFunctions", "Before");

        for (unsigned t = 0; t < 50; ++t)
        {
            SomeForce(vertex_mesh);
            builder.WriteVtkWithSubfolder("MonolayerCustomFunctions", "time_" + boost::lexical_cast<std::string>(t));

            PRINT_2_VARIABLES(t, vertex_mesh.GetVolumeOfElement(0));
        }
    }
};

#endif /* TESTMONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_ */
