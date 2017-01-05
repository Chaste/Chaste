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
        const unsigned node_indices_elem_0[4] = {0, 1, 4, 3};
        const unsigned node_indices_elem_1[4] = {1, 2, 5, 4};
        const unsigned node_indices_elem_2[5] = {4, 5, 6, 9, 8};
        const unsigned node_indices_elem_3[4] = {6, 7, 10, 9};

        MonolayerVertexMeshGenerator builder(nodes, "TestGetSharedElementnFaceIndices");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(5, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder("MonolayerCustomFunctions","Before");

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
        for (unsigned face_index=0; face_index<vertex_mesh.GetNumFaces(); ++face_index)
        {
            const VertexElement<2,3>* p_face = vertex_mesh.GetFace(face_index);
            TS_ASSERT_EQUALS(IsFaceOnBoundary(p_face), non_boundary_indices.count(face_index)==0);
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
        const unsigned node_indices_elem_0[4] = {0, 1, 3, 2};

        MonolayerVertexMeshGenerator builder(nodes, "TestGradientDeviation");
        builder.BuildElementWith(4, node_indices_elem_0);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder("MonolayerCustomFunctions","Before");
    }
};

#endif /* TESTMONOLAYERVERTEXMESHCUSTOMFUNCTIONS_HPP_ */
