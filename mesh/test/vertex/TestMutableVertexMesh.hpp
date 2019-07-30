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

#ifndef TESTMUTABLEVERTEXMESH_HPP_
#define TESTMUTABLEVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"
#include "MutableVertexMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ArchiveOpener.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestMutableVertexMesh : public CxxTest::TestSuite
{
private:

    MutableVertexMesh<3,3>* ConstructCubeAndPyramidMesh()
    {
        // Make 8 nodes to assign to a cube and a pyramid element
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(8, false, 0.5, 0.5, 1.5));

        std::vector<std::vector<Node<3>*> > nodes_faces(10);

        // Make 6 square faces and 4 triangular faces out of these nodes
        unsigned node_indices_face_0[4] = {0, 2, 4, 1};
        unsigned node_indices_face_1[4] = {4, 7, 5, 2};
        unsigned node_indices_face_2[4] = {7, 6, 1, 4};
        unsigned node_indices_face_3[4] = {0, 3, 5, 2};
        unsigned node_indices_face_4[4] = {1, 6, 3, 0};
        unsigned node_indices_face_5[4] = {7, 6, 3, 5};
        unsigned node_indices_face_6[3] = {6, 7, 8};
        unsigned node_indices_face_7[3] = {6, 8, 3};
        unsigned node_indices_face_8[3] = {3, 8, 5};
        unsigned node_indices_face_9[3] = {5, 8, 7};
        for (unsigned i=0; i<4; i++)
        {
            nodes_faces[0].push_back(nodes[node_indices_face_0[i]]);
            nodes_faces[1].push_back(nodes[node_indices_face_1[i]]);
            nodes_faces[2].push_back(nodes[node_indices_face_2[i]]);
            nodes_faces[3].push_back(nodes[node_indices_face_3[i]]);
            nodes_faces[4].push_back(nodes[node_indices_face_4[i]]);
            nodes_faces[5].push_back(nodes[node_indices_face_5[i]]);
            if (i < 3)
            {
                nodes_faces[6].push_back(nodes[node_indices_face_6[i]]);
                nodes_faces[7].push_back(nodes[node_indices_face_7[i]]);
                nodes_faces[8].push_back(nodes[node_indices_face_8[i]]);
                nodes_faces[9].push_back(nodes[node_indices_face_9[i]]);
            }
        }

        // Make the faces
        std::vector<VertexElement<2,3>*> faces;

        for (unsigned i=0; i<10; i++)
        {
            faces.push_back(new VertexElement<2,3>(i, nodes_faces[i]));
        }

        // Make the elements
        std::vector<VertexElement<2,3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;

        // Cube element
        for (unsigned i=0; i<6; i++)
        {
            faces_element_0.push_back(faces[i]);
            orientations_0.push_back(true);
        }

        // Pyramid element
        for (unsigned i=6; i<10; i++)
        {
            faces_element_1.push_back(faces[i]);
            orientations_1.push_back(true);
        }
        faces_element_1.push_back(faces[5]);
        orientations_1.push_back(false);

        std::vector<VertexElement<3,3>*> elements;
        elements.push_back(new VertexElement<3,3>(0, faces_element_0, orientations_0));
        elements.push_back(new VertexElement<3,3>(1, faces_element_1, orientations_1));

        return new MutableVertexMesh<3,3>(nodes, elements);
    }

public:

    void TestMutableVertexElementIterator()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 9u);

        unsigned counter = 0;
        for (MutableVertexMesh<2,2>::VertexElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumElements(), counter);

        // For coverage, test with an empty mesh
        MutableVertexMesh<2,2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        MutableVertexMesh<2,2>::VertexElementIterator el_iter = empty_mesh.GetElementIteratorBegin();

        // Check that the iterator is now at the end (we need to check this as a double-negative,
        // as we only have a NOT-equals operator defined on the iterator).
        bool iter_is_not_at_end = (el_iter != empty_mesh.GetElementIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);

        // Delete an element from mesh and test the iterator
        mesh.DeleteElementPriorToReMesh(0);

        counter = 0;
        for (MutableVertexMesh<2,2>::VertexElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter+1, element_index); // assumes the iterator will give elements 1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumElements(), counter);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), counter+1);
        TS_ASSERT_EQUALS(mesh.IsMeshChanging(), true);
    }

    void TestBasic2dMutableVertexMesh()
    {
        // Make seven nodes to assign to two elements
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        basic_nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = {0, 1, 2, 3, 4};
        unsigned node_indices_elem_1[3] = {2, 5, 6};
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem_0.push_back(basic_nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(basic_nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 7u);

        TS_ASSERT_DELTA(basic_vertex_mesh.GetNode(2)->rGetLocation()[0], 1.5, 1e-3);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);

        // Check that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0);

        // Nodes 1 and 4 are only in element 0
        TS_ASSERT_EQUALS(basic_nodes[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(basic_nodes[4]->rGetContainingElementIndices(), temp_list1);

        // Node 2 is in elements 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[2]->rGetContainingElementIndices(), temp_list1);

        // Node 5 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[5]->rGetContainingElementIndices(), temp_list2);

        // Test Set and Get methods
        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.01, 1e-4); // Default value
        TS_ASSERT_DELTA(basic_vertex_mesh.GetT2Threshold(), 0.001, 1e-4); // Default value

        basic_vertex_mesh.SetCellRearrangementThreshold(0.03);
        basic_vertex_mesh.SetT2Threshold(0.003);

        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.03, 1e-4);
        TS_ASSERT_DELTA(basic_vertex_mesh.GetT2Threshold(), 0.003, 1e-4);

        // Coverage
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveBoundaryElementMapping(0), 0u);
    }

    void TestBasic3dMutableVertexMesh()
    {
        MutableVertexMesh<3,3>* p_mesh = ConstructCubeAndPyramidMesh();

        // Set/get parameter values
        p_mesh->SetCellRearrangementThreshold(0.1);
        p_mesh->SetT2Threshold(0.01);
        p_mesh->SetCellRearrangementRatio(1.6);

        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementThreshold(), 0.1, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetT2Threshold(), 0.01, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementRatio(), 1.6, 1e-6);

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 2u);

        // Test the location of one of the nodes
        Node<3>* p_node_2 = p_mesh->GetNode(2);
        TS_ASSERT_DELTA(p_node_2->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(p_node_2->rGetLocation()[1], 1.0, 1e-3);
        TS_ASSERT_DELTA(p_node_2->rGetLocation()[2], 0.0, 1e-3);

        // Test a couple of the elements
        VertexElement<3,3>* p_element_0 = p_mesh->GetElement(0);
        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetNumFaces(), 6u);

        VertexElement<3,3>* p_element_1 = p_mesh->GetElement(1);
        TS_ASSERT_EQUALS(p_element_1->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_element_1->GetNumFaces(), 5u);

        // Check that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0);

        // Nodes 0, 1, 2 and 4 are only in element 0
        TS_ASSERT_EQUALS(p_mesh->GetNode(0)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(1)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(2)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(4)->rGetContainingElementIndices(), temp_list1);

        // Node 3, 5, 6 and 7 are in elements 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(3)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(5)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(6)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(7)->rGetContainingElementIndices(), temp_list1);

        // Node 8 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(8)->rGetContainingElementIndices(), temp_list2);

        // Coverage
        TS_ASSERT_EQUALS(p_mesh->SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveBoundaryElementMapping(0), 0u);

        // Tidy up
        delete p_mesh;
    }

    void TestMeshConstructionFromMeshReader()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Test Get methods
        TS_ASSERT_DELTA(mesh.GetCellRearrangementThreshold(), 0.01, 1e-4); // Default value
        TS_ASSERT_DELTA(mesh.GetT2Threshold(), 0.001, 1e-4); // Default value

        // Check we have the right number of nodes and elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        // Check second element has the right nodes
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
    }

    void TestSetNode()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<2> point = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 2.0, 1e-6);

        // Nudge node
        point.SetCoordinate(0, 1.1);
        mesh.SetNode(3, point);

        ChastePoint<2> point2 = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point2[0], 1.1, 1e-6);
        TS_ASSERT_DELTA(point2[1], 2.0, 1e-6);

        // Nudge node again
        point.SetCoordinate(1, 1.9);
        mesh.SetNode(3, point);

        ChastePoint<2> point3 = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point3[0], 1.1, 1e-6);
        TS_ASSERT_DELTA(point3[1], 1.9, 1e-6);
    }

    void TestAddNodeAndReMesh()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Choose a node
        ChastePoint<2> point = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 2.0, 1e-6);

        // Create a new node close to this node
        point.SetCoordinate(0, 1.1);
        point.SetCoordinate(1, 2.1);
        Node<2>* p_node = new Node<2>(mesh.GetNumNodes(), point);

        unsigned old_num_nodes = mesh.GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = mesh.AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Remesh to update correspondences
        VertexElementMap map(mesh.GetNumElements());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[0], 1.1, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[1], 2.1, 1e-7);

        // Now test AddNode() when mDeletedNodeIndices is populated

        // Label node 3 as deleted
        mesh.mDeletedNodeIndices.push_back(3);

        // Create a new node close to this node
        ChastePoint<2> point2;
        point2.SetCoordinate(0, 0.9);
        point2.SetCoordinate(1, 1.9);
        Node<2>* p_node2 = new Node<2>(mesh.GetNumNodes(), point);

        // Add this new node to the mesh
        new_index = mesh.AddNode(p_node2);
        TS_ASSERT_EQUALS(new_index, 3u);
    }

    void TestAddElement()
    {
        // Make four nodes to assign to two elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = {0, 1, 2, 3, 4};
        unsigned node_indices_elem_1[3] = {2, 5, 6};
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<VertexElement<2,2>*> elements;
        VertexElement<2,2>* p_replaced_vertex_element = new VertexElement<2,2>(0, nodes_elem_0);
        elements.push_back(p_replaced_vertex_element);
        elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        std::vector<Node<2>*> nodes_elem_2, nodes_elem_3;

        // Make two triangular elements out of these nodes
        nodes_elem_2.push_back(nodes[6]);
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[2]);

        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[1]);
        nodes_elem_3.push_back(nodes[2]);
        nodes_elem_3.push_back(nodes[3]);

        // Add a new element to the mesh
        mesh.AddElement(new VertexElement<2,2>(2, nodes_elem_2));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);

        // Replace element 0 in the mesh
        mesh.AddElement(new VertexElement<2,2>(0, nodes_elem_3));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);

        // Tidy up
        delete p_replaced_vertex_element;
    }

    void TestDeletingNodes()
    {
        // Make a simple vertex mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = {0, 1, 2, 3, 4};
        unsigned node_indices_elem_1[3] = {2, 5, 6};
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        MutableVertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        mesh.DeleteElementPriorToReMesh(0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
    }

    void TestDivideVertexElementGivenNodes()
    {
        // Make four nodes
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        basic_nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        basic_nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, -2.0, -1.0));

        // Make one rectangular element out of these nodes. Ordering for coverage.
        std::vector<Node<2>*> nodes_elem;
        nodes_elem.push_back(basic_nodes[2]);
        nodes_elem.push_back(basic_nodes[3]);
        nodes_elem.push_back(basic_nodes[0]);
        nodes_elem.push_back(basic_nodes[1]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 4u);

        // Divide element using two given nodes, putting the original element below the new element
        unsigned new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(0), 2, 0, true);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 0u);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 1u);

        // For coverage, divide an element when mDeletedElementIndices is not empty
        basic_vertex_mesh.DeleteElementPriorToReMesh(0);
        new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(1), 2, 3, true);

        TS_ASSERT_EQUALS(new_element_index, 0u);
    }

    void TestDivideVertexElementGivenNodesForCoverage()
    {
        /*
         * Divide a square element like so
         *   ___      ___
         *  |   |    |  /|
         *  |   | -> | / |
         *  |___|    |/__|
         */

        // Make four nodes
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        basic_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 0.0, 0.0));

        // Make one rectangular element out of these nodes. Ordering for coverage.
        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem.push_back(basic_nodes[i]);
        }

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 4u);

        // Divide element using two given nodes, putting the original element above the new element (for coverage)
        unsigned new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(0), 1, 3, false);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNumNodes(), 3u);

        // Test elements have correct nodes
        unsigned expected_node_indices_element_0[3] = {1, 2, 3};
        unsigned expected_node_indices_element_1[3] = {0, 1, 3};
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), expected_node_indices_element_0[i]);
            TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), expected_node_indices_element_1[i]);
        }

        // For coverage, divide an element when mDeletedElementIndices is not empty
        basic_vertex_mesh.DeleteElementPriorToReMesh(0);
        new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(1), 2, 3, true);

        TS_ASSERT_EQUALS(new_element_index, 0u);
    }


    void TestDivideVertexElementAbove()
    {
        // Make four nodes
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        basic_nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        basic_nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, -2.0, -1.0));

        // Make one rectangular element out of these nodes. Ordering for coverage.
        std::vector<Node<2>*> nodes_elem;
        nodes_elem.push_back(basic_nodes[2]);
        nodes_elem.push_back(basic_nodes[3]);
        nodes_elem.push_back(basic_nodes[0]);
        nodes_elem.push_back(basic_nodes[1]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 4u);

        // Divide element using two given nodes, putting the original element above the new element
        unsigned new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(0), 2, 0, false);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 1u);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 0u);
    }

    // This also tests that boundary nodes are updated on element division
    void TestDivideVertexElementGivenAxisOfDivision()
    {
        // Make five nodes, 0, 1 and 2 are boundary nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 1.0, -2.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 2.0));
        nodes.push_back(new Node<2>(2, true, -1.0, 2.0));
        nodes.push_back(new Node<2>(3, false, -1.0, -2.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 3.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_0;
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[i]);
        }

        // Make a triangular element out of nodes 1,4,2
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[2]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);

        c_vector<double, 2> axis_of_division;
        axis_of_division(0) = 1.0;
        axis_of_division(1) = 0.0;

        // Divide element 0 along given axis
        unsigned new_element_index = vertex_mesh.DivideElementAlongGivenAxis(vertex_mesh.GetElement(0), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, vertex_mesh.GetNumElements()-1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Now test the position of new nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.0, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], -1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.0, 1e-8);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);

        // Now test the nodes in each element
        unsigned expected_node_indices_element_0[4] = {0, 5, 6, 3};
        unsigned expected_node_indices_element_1[3] = {1, 4, 2};
        unsigned expected_node_indices_element_2[4] = {5, 1, 2, 6};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), expected_node_indices_element_0[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), expected_node_indices_element_1[i]);
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), expected_node_indices_element_2[i]);
        }

        // Test boundary nodes updated
        for (unsigned i=0; i<7; i++)
        {
            bool expected_boundary_node = (i!=3 && i!=4 && i!=6);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->rGetContainingElementIndices(), expected_elements_containing_node_6);
    }

    void TestDivideVertexElementWithBoundaryNodes()
    {

        /*
         * This test checks that the new node created in the centre of the mesh is not a boundary node.
         *  _________       _________
         * |    |    |     |    |    |
         * |    |    | --> |____|    |
         * |    |    |     |    |    |
         * |____|____|     |____|____|
         *
         */

        // Make five nodes, all boundary nodes.
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(5, true, 2.0, 1.0));

        // Make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 2, 3};
        unsigned node_indices_elem_1[4] = {1, 4, 5, 2};
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        c_vector<double, 2> axis_of_division;
        axis_of_division(0) = 1.0;
        axis_of_division(1) = 0.0;

        // Divide element 0 along given axis
        unsigned new_element_index = vertex_mesh.DivideElementAlongGivenAxis(vertex_mesh.GetElement(0), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, vertex_mesh.GetNumElements()-1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Now test the position of new nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.5, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(7)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(7)->rGetLocation()[1], 0.5, 1e-8);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);

        // Now test the nodes in each element
        unsigned expected_node_indices_element_0[4] = {0, 1 , 6, 7};
        unsigned expected_node_indices_element_1[5] = {1, 4, 5, 2, 6};
        unsigned expected_node_indices_element_2[4] = {6, 2, 3, 7};
        for (unsigned i=0; i<5; i++)
        {
            if (i < 4)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), expected_node_indices_element_0[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), expected_node_indices_element_2[i]);
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), expected_node_indices_element_1[i]);
        }

        // Test boundary nodes updated
        for (unsigned i=0; i<8; i++)
        {
            bool expected_boundary_node = (i!=6);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(1);
        expected_elements_containing_node_6.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->rGetContainingElementIndices(), expected_elements_containing_node_6);

        std::set<unsigned> expected_elements_containing_node_7;
        expected_elements_containing_node_7.insert(0);
        expected_elements_containing_node_7.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(7)->rGetContainingElementIndices(), expected_elements_containing_node_7);
    }

    void TestDeleteElementWithBoundaryNodes()
    {

        /*
         * This test checks node 'boundaryness' when a boundary element is deleted.
         *  ____       ____
         * | /  |     |    |
         * |/   | --> |    |
         * |    |     |    |
         * |____|     |____|
         *
         */

        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.5, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.5, 1.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(6, true, 0.0, 0.5));
        nodes.push_back(new Node<2>(7, false, 0.4, 0.6));

        // Create elements
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[7] = {0, 1, 2, 3, 4, 7, 6};
        unsigned node_indices_elem_1[4] = {6, 7, 4, 5};
        for (unsigned i=0; i<7; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 4)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Create mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Test mesh has correct numbers of elements and nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Test correct nodes are boundary nodes
        for (unsigned i=0; i<8; i++)
        {
            bool expected_boundary_node = (i!=7);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

        // Delete smaller element
        vertex_mesh.DeleteElementPriorToReMesh(1);

        // Test mesh has correct numbers of elements and nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Test ownership
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 7u);
        unsigned expected_node_indices[7] = {0, 1, 2, 3, 4, 7, 6};
        for (unsigned i=0; i<7; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), expected_node_indices[i]);
        }

        // Test correct nodes are boundary nodes
        for (unsigned i=0; i<7; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    /**
     * Test that in the case where the given axis of division does not
     * cross two edges of the element, an exception is thrown.
     */
    void TestDivideVertexElementGivenAxisOfDivisionFailsForBadElement()
    {
        // Create a mesh consisting of a single non-convex element
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.4, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.4, 1.0));
        nodes.push_back(new Node<2>(3, false, 1.2, 1.0));
        nodes.push_back(new Node<2>(4, false, 1.2, 0.2));
        nodes.push_back(new Node<2>(5, false, 1.0, 0.2));
        nodes.push_back(new Node<2>(6, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(7, false, 0.0, 1.0));

        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            nodes_elem.push_back(nodes[i]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Provide an axis of division that does not cross two edges of the element (it crosses four)
        c_vector<double, 2> axis_of_division;
        axis_of_division(0) = 1.0;
        axis_of_division(1) = 0.0;

        // Divide element 0 along given axis
        TS_ASSERT_THROWS_THIS(vertex_mesh.DivideElementAlongGivenAxis(vertex_mesh.GetElement(0), axis_of_division, true),
                              "Cannot proceed with element division: the given axis of division does not cross two edges of the element");
    }

    void TestDivideVertexElementAlongShortAxis()
    {
        // Make five nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, false, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 2.0));

        // Make a rectangular element and a triangular element
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 2, 3};
        unsigned node_indices_elem_1[3] = {1, 4, 2};
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);

        // Divide element 0 along short axis
        unsigned new_element_index = vertex_mesh.DivideElementAlongShortAxis(vertex_mesh.GetElement(0), true);

        TS_ASSERT_EQUALS(new_element_index, vertex_mesh.GetNumElements()-1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Now test the position of new nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 1.0, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], -1.0, 1e-8);

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        unsigned expected_node_indices_element_0[4] = {0, 1, 5, 6};
        unsigned expected_node_indices_element_1[4] = {1, 4, 2, 5};
        unsigned expected_node_indices_element_2[4] = {5, 2, 3, 6};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), expected_node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), expected_node_indices_element_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), expected_node_indices_element_2[i]);
        }

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);
        expected_elements_containing_node_5.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->rGetContainingElementIndices(), expected_elements_containing_node_6);
    }

    void TestDivideVertexElementWithNonRegularElement()
    {
        // Make six nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, 3.0, 2.0));
        nodes.push_back(new Node<2>(3, false, 3.0, 3.0));
        nodes.push_back(new Node<2>(4, false, 1.0, 2.0));

        // Make one element out of these nodes
        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem.push_back(nodes[i]);
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);

        // Test that the centroid, moments and short axis of the element are calculated correctly
        // (i.e. agreee with Matlab and pen-and-paper calculations)
        c_vector<double,2> centroid = mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(centroid[0], 2.0000, 1e-4);
        TS_ASSERT_DELTA(centroid[1], 1.8666, 1e-4);

        c_vector<double,3> moments = mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(moments[0], 0.5388, 1e-4);
        TS_ASSERT_DELTA(moments[1], 0.7500, 1e-4);
        TS_ASSERT_DELTA(moments[2], 0.3750, 1e-4);

        c_vector<double,2> short_axis = mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(short_axis[0],  0.6037, 1e-4);
        TS_ASSERT_DELTA(short_axis[1], -0.7971, 1e-4);

        // Divide element using two given nodes
        unsigned new_element_index = mesh.DivideElementAlongShortAxis(mesh.GetElement(0), true);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 5u);
        unsigned expected_node_indices_element_0[5] = {0, 1, 5, 6, 4};
        for (unsigned i=0; i<5; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), expected_node_indices_element_0[i]);
        }

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(3), 6u);

        // Test locations of new nodes
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[0], 2.3735, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[1], 1.3735, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 1.6520, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 2.3260, 1e-4);
    }

    void TestDivideVertexElementWhereNewNodesAreCloseToOldNodes1()
    {
        // Make 6 nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, -1.0, 0.0));
        nodes.push_back(new Node<2>(1, false, -0.009, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 1.0));
        nodes.push_back(new Node<2>(5, false, -1.0, 1.0));

        // Make one rectangular element out of these nodes
        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<6; i++)
        {
            nodes_elem.push_back(nodes[i]);
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);

        // Divide element
        unsigned new_element_index = mesh.DivideElementAlongShortAxis(mesh.GetElement(0), true);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 8u);

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 5u);

        // Test elements have correct nodes
        unsigned node_indices_element_0[5] = {0, 1, 6, 7, 5};
        unsigned node_indices_element_1[5] = {6, 2, 3, 4, 7};
        for (unsigned i=0; i<5; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
        }

        // Test locations of new nodes
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], -0.009+1.5*mesh.GetCellRearrangementThreshold(), 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[1], 1.0, 1e-4);
    }

    void TestDivideVertexElementWhereNewNodesAreCloseToOldNodes2()
    {
        // Make 6 nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, -1.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.009, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 1.0));
        nodes.push_back(new Node<2>(5, false, -1.0, 1.0));

        // Make one rectangular element out of these nodes
        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<6; i++)
        {
            nodes_elem.push_back(nodes[i]);
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);

        // Divide element
        unsigned new_element_index = mesh.DivideElementAlongShortAxis(mesh.GetElement(0), true);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 8u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(1), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(3), 5u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(3), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(4), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(5), 7u);

        // Test locations of new nodes
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 0.009-1.5*mesh.GetCellRearrangementThreshold(), 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[1], 1.0, 1e-4);
    }

    void TestDivideVertexElementGivenAxisOfDivisionWithShortEdge()
    {
        // Make five nodes, 0, 1 and 2 are boundary nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 1.0, -0.1));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.12));
        nodes.push_back(new Node<2>(2, true, -1.0, 3.0));
        nodes.push_back(new Node<2>(3, true, -1.0, -3.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_0;
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[i]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.09);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        c_vector<double, 2> axis_of_division;
        axis_of_division(0) = 1.0;
        axis_of_division(1) = 0.0;

        // Divide element 0 along given axis
        unsigned new_element_index = vertex_mesh.DivideElementAlongGivenAxis(vertex_mesh.GetElement(0), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, vertex_mesh.GetNumElements()-1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Now test the position of new nodes 4 is at the midpoint of 0 and 1 as 0-1 is less than the threshold distance
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.01, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], -1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.00345123, 1e-8); // height of centroid of original element

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        unsigned node_indices_element_0[4] = {0, 4, 5, 3};
        unsigned node_indices_element_1[4] = {4, 1, 2, 5};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
        }

        // Test boundary nodes updated
        for (unsigned i=0; i<6; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_4;
        expected_elements_containing_node_4.insert(0);
        expected_elements_containing_node_4.insert(1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->rGetContainingElementIndices(), expected_elements_containing_node_4);

        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);
    }

    void TestArchive2dMutableVertexMesh()
    {
        // Set archiving location
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mutable_vertex_2d.arch";
        ArchiveLocationInfo::SetMeshFilename("mutable_vertex_2d");

        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_5_by_3");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set member variables
        mesh.SetCellRearrangementThreshold(0.54);
        mesh.SetT2Threshold(0.012);
        mesh.SetCellRearrangementRatio(1.6);
        mesh.SetDistanceForT3SwapChecking(7.3);

        AbstractMesh<2,2>* const p_mesh = &mesh;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects
         * changing during the save, and so object tracking leading to wrong results.
         *
         * E.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.
         */

        // Create an output archive
        {
            TS_ASSERT_EQUALS((static_cast<MutableVertexMesh<2,2>*>(p_mesh))->GetNumNodes(), 46u);
            TS_ASSERT_EQUALS((static_cast<MutableVertexMesh<2,2>*>(p_mesh))->GetNumElements(), 15u);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost
            (*p_arch) << p_mesh;
        }

        {
            // De-serialize and compare
            AbstractMesh<2,2>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            MutableVertexMesh<2,2>* p_mesh_original = static_cast<MutableVertexMesh<2,2>*>(p_mesh2);
            MutableVertexMesh<2,2>* p_mesh_loaded = static_cast<MutableVertexMesh<2,2>*>(p_mesh);

            // Test member variables were archived correctly
            TS_ASSERT_DELTA(p_mesh_original->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementRatio(), 1.6, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetDistanceForT3SwapChecking(), 7.3, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetDistanceForT3SwapChecking(), 7.3, 1e-6);

            // Compare the loaded mesh against the original
            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), p_mesh_loaded->GetNumNodes());

            for (unsigned node_index=0; node_index<p_mesh_original->GetNumNodes(); node_index++)
            {
                Node<2>* p_node = p_mesh_original->GetNode(node_index);
                Node<2>* p_node2 = p_mesh_loaded->GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned dimension=0; dimension<2; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), p_mesh_loaded->GetNumElements());

            for (unsigned elem_index=0; elem_index < p_mesh_original->GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumNodes(),
                                 p_mesh_loaded->GetElement(elem_index)->GetNumNodes());

                for (unsigned local_index=0; local_index<p_mesh_original->GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNodeGlobalIndex(local_index),
                                     p_mesh_loaded->GetElement(elem_index)->GetNodeGlobalIndex(local_index));
                }
            }

            // Tidy up
            delete p_mesh_original;
        }
    }

    void TestArchive3dMutableVertexMesh()
    {
        // Set archiving location
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mutable_vertex_3d.arch";
        ArchiveLocationInfo::SetMeshFilename("mutable_vertex_3d");

        // Create mesh
        MutableVertexMesh<3,3>* p_mutable_mesh = ConstructCubeAndPyramidMesh();

        // Set member variables
        p_mutable_mesh->SetCellRearrangementThreshold(0.54);
        p_mutable_mesh->SetT2Threshold(0.012);
        p_mutable_mesh->SetCellRearrangementRatio(1.6);

        AbstractMesh<3,3>* const p_mesh = p_mutable_mesh;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects
         * changing during the save, and so object tracking leading to wrong results.
         *
         * E.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.
         */

        // Create an output archive
        {
            TS_ASSERT_EQUALS((static_cast<MutableVertexMesh<3,3>*>(p_mesh))->GetNumNodes(), 9u);
            TS_ASSERT_EQUALS((static_cast<MutableVertexMesh<3,3>*>(p_mesh))->GetNumElements(), 2u);
            TS_ASSERT_EQUALS((static_cast<MutableVertexMesh<3,3>*>(p_mesh))->GetNumFaces(), 10u);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost
            (*p_arch) << p_mesh;
        }

        {
            // De-serialize and compare
            AbstractMesh<3,3>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            MutableVertexMesh<3,3>* p_mesh_original = static_cast<MutableVertexMesh<3,3>*>(p_mesh);
            MutableVertexMesh<3,3>* p_mesh_loaded = static_cast<MutableVertexMesh<3,3>*>(p_mesh2);

            // Compare the loaded mesh against the original
            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), p_mesh_loaded->GetNumNodes());

            for (unsigned node_index=0; node_index<p_mesh_original->GetNumNodes(); node_index++)
            {
                Node<3>* p_node = p_mesh_original->GetNode(node_index);
                Node<3>* p_node2 = p_mesh_loaded->GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned dimension=0; dimension<3; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), p_mesh_loaded->GetNumElements());

            for (unsigned elem_index=0; elem_index < p_mesh_original->GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumNodes(),
                                 p_mesh_loaded->GetElement(elem_index)->GetNumNodes());

                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumFaces(),
                                 p_mesh_loaded->GetElement(elem_index)->GetNumFaces());

                for (unsigned local_index=0; local_index<p_mesh_original->GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNodeGlobalIndex(local_index),
                                     p_mesh_loaded->GetElement(elem_index)->GetNodeGlobalIndex(local_index));
                }
            }

            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementRatio(), 1.6, 1e-6);

            TS_ASSERT_DELTA(p_mesh_original->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetCellRearrangementRatio(), 1.6, 1e-6);

            // Tidy up
            delete p_mesh_loaded;
        }

        delete p_mesh;
    }

    ///\todo include boundary nodes in the tests
    void TestDivideEdge()
    {
        /*
         *     Element
         *   0    2     1
         *
         *    3________2
         *    /|      |\
         * 4 / |      | \ 5
         *   \ |      | /
         *    \|______|/
         *    0        1
         */

        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(3, true, 1.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 0.5, 1.5));
        nodes.push_back(new Node<2>(5, true, 2.5, 1.5));

        // Create elements
        std::vector<Node<2>*> nodes_in_element0, nodes_in_element1, nodes_in_element2;
        unsigned node_indices_in_element0[3] = {0, 3, 4};
        unsigned node_indices_in_element1[3] = {1, 5, 2};
        unsigned node_indices_in_element2[4] = {0, 1, 2, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_in_element0.push_back(nodes[node_indices_in_element0[i]]);
                nodes_in_element1.push_back(nodes[node_indices_in_element1[i]]);
            }
            nodes_in_element2.push_back(nodes[node_indices_in_element2[i]]);
        }

        /*
         *  Create three elements, elements0 and 2 share nodes 0 and 3,
         *  and elements 1 and 2 share nodes 1 and 2
         */

        VertexElement<2,2>* p_element0 = new VertexElement<2,2>(0, nodes_in_element0);
        VertexElement<2,2>* p_element1 = new VertexElement<2,2>(1, nodes_in_element1);
        VertexElement<2,2>* p_element2 = new VertexElement<2,2>(2, nodes_in_element2);
        TS_ASSERT_EQUALS(p_element0->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_element1->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_element2->GetNumNodes(), 4u);

        // Create mesh

        std::vector<VertexElement<2,2>* > elements;
        elements.push_back(p_element0);
        elements.push_back(p_element1);
        elements.push_back(p_element2);

        MutableVertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 4u);

        // Divide the edge joining nodes 0 and 1
        mesh.DivideEdge(mesh.GetNode(0), mesh.GetNode(1));

        // Test edge is divided
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);

        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(2), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(2), 4.0, 1e-6);

        // Test other nodes are updated
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(1), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 4u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 2u);

        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(1), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(3), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(4), 3u);

        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[1], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[0], 2.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[1], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[1], 1.0, 1e-9);

        // Divide the edge joining nodes 3 and 0
        mesh.DivideEdge(mesh.GetNode(3), mesh.GetNode(0)); // Ordering for coverage

        // Divide the edge joining nodes 2 and 1
        mesh.DivideEdge(mesh.GetNode(1), mesh.GetNode(2)); // Ordering for coverage

        // Test edges are divided
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(2), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(2), 4.0, 1e-6);

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 7u);

        // Test other nodes are updated
        unsigned expected_node_indices_element_0[4] = {0, 7, 3, 4};
        unsigned expected_node_indices_element_1[4] = {1, 5, 2, 8};
        unsigned expected_node_indices_element_2[7] = {0, 6, 1, 8, 2, 3, 7};
        for (unsigned i=0; i<7; i++)
        {
            if (i < 4)
            {
                TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), expected_node_indices_element_0[i]);
                TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(i), expected_node_indices_element_1[i]);
            }
            TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(i), expected_node_indices_element_2[i]);
        }

        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[1], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[0], 2.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[1], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(7)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(7)->GetPoint()[1], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[1], 1.5, 1e-9);

        // Test boundary property of nodes. All are boundary nodes except nodes 7 and 8.
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = (i!=7 && i!=8);
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestElementIncludesPointAndGetLocalIndexForElementEdgeClosestToPoint()
    {
        /*
         * Create a simple mesh comprising four nodes contained in a square element.
         * We will test for inclusion in the element for several points.
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        MutableVertexMesh<2,2> mesh(nodes, elements);

        // Test a point far outside the element
        c_vector<double, 2> test_point1;
        test_point1[0] = -1.0;
        test_point1[1] = -1.0;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point1, 0), false);

        // Test a point far inside the element
        c_vector<double, 2> test_point2;
        test_point2[0] = 0.5;
        test_point2[1] = 0.5;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point2, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point2, 0), 0u);

        // Test a point on the left edge
        c_vector<double, 2> test_point3;
        test_point3[0] = 0.0;
        test_point3[1] = 0.5;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point3, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point3, 0), 3u);

        // Test a point on the right edge
        c_vector<double, 2> test_point4;
        test_point4[0] = 1.0;
        test_point4[1] = 0.5;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point4, 0), false);

        // A point on the bottom edge
        c_vector<double, 2> test_point5;
        test_point5[0] = 0.5;
        test_point5[1] = 0.0;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point5, 0), false); ///\todo #2387
//        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point5, 0), true);

        // Test a point on the top edge
        c_vector<double, 2> test_point6;
        test_point6[0] = 0.5;
        test_point6[1] = 1.0;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point6, 0), false);

        // Test a point just inside the element
        c_vector<double, 2> test_point7;
        test_point7[0] = 0.999;
        test_point7[1] = 0.998;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point7, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point7, 0), 1u);

        // Test a point just outside the element
        c_vector<double, 2> test_point8;
        test_point8[0] = 1.001;
        test_point8[1] = 0.5;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point8, 0), false);

        // Test a point coinciding with a vertex
        c_vector<double, 2> test_point9;
        test_point9[0] = 1.0;
        test_point9[1] = 1.0;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point9, 0), false);
        unsigned edge_closest_to_point_9 = mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point9,0);
        TS_ASSERT(edge_closest_to_point_9 == 1u || edge_closest_to_point_9 == 2u);
    }

    void TestGetLocalIndexForElementEdgeClosestToPointWithinDoubleEpsilon()
    {
        /*
         * Create a simple mesh comprising four nodes contained in a square element.
         * We will test the GetLocalIndexForElementEdgeClosestToPoint method for a point
         * which is within DBL_EPSILON distance of but not exactly at one of the nodes.
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.1, 0.1));
        nodes.push_back(new Node<2>(3, false, 0.0, 0.1));

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        MutableVertexMesh<2,2> mesh(nodes, elements);

        // Create a point within DBL_EPSILON of node 2
        c_vector<double, 2> test_point1;
        test_point1[0] = 0.1+2e-16;
        test_point1[1] = 0.1+1e-16;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point1, 0), false);

        // Since the point is basically at node 2 we are happy with the method returning either edge 1
        // or edge 2.  The closest edge is not well defined in this situation.
        unsigned edge_closest_to_point1 = mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point1,0);
        TS_ASSERT(edge_closest_to_point1 == 1u || edge_closest_to_point1 == 2u);
    }

    void TestSetAndGetDistanceForT3SwapChecking()
    {
        HoneycombVertexMeshGenerator mesh_generator(10,10); // is_flat_bottom, T1swaptthreshold, T2swapthreshold, elementArea
        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        double standard_distance = p_mesh->GetDistanceForT3SwapChecking();
        TS_ASSERT_EQUALS(standard_distance, 5.0);

        p_mesh->SetDistanceForT3SwapChecking( 10.0 );
        double new_distance = p_mesh->GetDistanceForT3SwapChecking();
        TS_ASSERT_EQUALS(new_distance, 10.0);
    }

    void TestHandleHighOrderJunctions()
    {
        //\todo need to re-implement this
    }
};

#endif /*TESTMUTABLEVERTEXMESH_HPP_*/
