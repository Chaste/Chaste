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

#ifndef TESTVERTEXMESH_HPP_
#define TESTVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexMesh.hpp"
#include "ArchiveOpener.hpp"
#include "MutableMesh.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestVertexMesh : public CxxTest::TestSuite
{
private:

    VertexMesh<3,3>* ConstructCubeAndPyramidMesh()
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

        // Make six square faces and four triangular faces out of these nodes
        std::vector<std::vector<Node<3>*> > nodes_faces(10);
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

        return new VertexMesh<3,3>(nodes, faces, elements);
    }

    VertexMesh<3,3>* ConstructPrismMesh()
    {
        // Create nodes
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 1.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 3.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 3.0));

        // Make five faces out of these nodes
        std::vector<std::vector<Node<3>*> > nodes_faces(5);

        nodes_faces[0].push_back(nodes[0]);
        nodes_faces[0].push_back(nodes[4]);
        nodes_faces[0].push_back(nodes[5]);
        nodes_faces[0].push_back(nodes[1]);

        nodes_faces[1].push_back(nodes[0]);
        nodes_faces[1].push_back(nodes[3]);
        nodes_faces[1].push_back(nodes[4]);

        nodes_faces[2].push_back(nodes[3]);
        nodes_faces[2].push_back(nodes[2]);
        nodes_faces[2].push_back(nodes[5]);
        nodes_faces[2].push_back(nodes[4]);

        nodes_faces[3].push_back(nodes[1]);
        nodes_faces[3].push_back(nodes[5]);
        nodes_faces[3].push_back(nodes[2]);

        nodes_faces[4].push_back(nodes[3]);
        nodes_faces[4].push_back(nodes[2]);
        nodes_faces[4].push_back(nodes[1]);
        nodes_faces[4].push_back(nodes[0]);

        std::vector<VertexElement<2,3>*> faces(5);
        std::vector<bool> orientations(5);
        for (unsigned i=0; i<5; i++)
        {
            faces[i] = new VertexElement<2,3>(i, nodes_faces[i]);
            orientations[i] = true;
        }

        // Create cuboidal element
        std::vector<VertexElement<3,3>*> elements;
        elements.push_back(new VertexElement<3,3>(0, faces, orientations));

        // Create mesh
        return new VertexMesh<3,3>(nodes, elements);
    }

public:

    void TestNodeIterator()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 30u);

        unsigned counter = 0;
        for (VertexMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            unsigned node_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, node_index); // assumes the iterator will give nodes 0,1..,N in that order
            counter++;
        }
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), counter);

        // Check that the node iterator correctly handles deleted nodes
        mesh.GetNode(0)->MarkAsDeleted();

        counter = 0;
        for (VertexMesh<2,2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            unsigned node_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter+1, node_index); // assumes the iterator will give nodes 1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), counter+1);

        // For coverage, test with an empty mesh
        VertexMesh<2,2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        VertexMesh<2,2>::NodeIterator iter = empty_mesh.GetNodeIteratorBegin();

        // Check that the iterator is now at the end (we need to check this as a double-negative,
        // as we only have a NOT-equals operator defined on the iterator).
        bool iter_is_not_at_end = (iter != empty_mesh.GetNodeIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);

        // Coverage of AbstractMesh::SetElementOwnerships()
        TS_ASSERT_THROWS_NOTHING(empty_mesh.SetElementOwnerships());
    }

    void TestVertexElementIterator()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 9u);

        unsigned counter = 0;
        for (VertexMesh<2,2>::VertexElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumElements(), counter);

        // For coverage, test with an empty mesh
        VertexMesh<2,2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        VertexMesh<2,2>::VertexElementIterator iter = empty_mesh.GetElementIteratorBegin();

        // Check that the iterator is now at the end (we need to check this as a double-negative,
        // as we only have a NOT-equals operator defined on the iterator).
        bool iter_is_not_at_end = (iter != empty_mesh.GetElementIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);

        // Check that the number of elements matches and that the mesh is not mutable
        TS_ASSERT_EQUALS(mesh.GetNumElements(), counter);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), counter);
        TS_ASSERT_EQUALS(mesh.IsMeshChanging(), false);
    }

    void TestBasic1dVertexMesh()
    {
        // Create a 1D mesh comprising four nodes and three elements
        std::vector<Node<1>*> nodes_1d;
        for (unsigned i=0; i<4; i++)
        {
            nodes_1d.push_back(new Node<1>(i, false, 0.5*(double)i));
        }

        std::vector<std::vector<Node<1>*> > nodes_elements_1d(3);
        std::vector<VertexElement<1,1>*> elements_1d;
        for (unsigned i=0; i<3; i++)
        {
            nodes_elements_1d[i].push_back(nodes_1d[i]);
            nodes_elements_1d[i].push_back(nodes_1d[i+1]);
            elements_1d.push_back(new VertexElement<1,1>(i, nodes_elements_1d[i]));
        }

        VertexMesh<1,1> mesh_1d(nodes_1d, elements_1d);

        // Test the mesh has the correct number of nodes and elements
        TS_ASSERT_EQUALS(mesh_1d.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh_1d.GetNumElements(), 3u);

        // Test the elements have the correct nodes
        TS_ASSERT_EQUALS(mesh_1d.GetElement(0)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh_1d.GetElement(0)->GetNodeLocation(0)[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetElement(0)->GetNodeLocation(1)[0], 0.5, 1e-6);

        TS_ASSERT_EQUALS(mesh_1d.GetElement(1)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh_1d.GetElement(1)->GetNodeLocation(0)[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetElement(1)->GetNodeLocation(1)[0], 1.0, 1e-6);

        TS_ASSERT_EQUALS(mesh_1d.GetElement(2)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh_1d.GetElement(2)->GetNodeLocation(0)[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetElement(2)->GetNodeLocation(1)[0], 1.5, 1e-6);
    }

    void TestBasic2dVertexMesh()
    {
        // Create a 2D mesh comprising seven nodes and two elements
        std::vector<Node<2>*> nodes_2d;
        nodes_2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes_2d.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes_2d.push_back(new Node<2>(2, false, 1.5, 1.0));
        nodes_2d.push_back(new Node<2>(3, false, 1.0, 2.0));
        nodes_2d.push_back(new Node<2>(4, false, 0.0, 1.0));
        nodes_2d.push_back(new Node<2>(5, false, 2.0, 0.0));
        nodes_2d.push_back(new Node<2>(6, false, 2.0, 3.0));

        std::vector<std::vector<Node<2>*> > nodes_elements_2d(2);
        for (unsigned i=0; i<5; i++)
        {
            nodes_elements_2d[0].push_back(nodes_2d[i]);
        }
        nodes_elements_2d[1].push_back(nodes_2d[2]);
        nodes_elements_2d[1].push_back(nodes_2d[5]);
        nodes_elements_2d[1].push_back(nodes_2d[6]);

        std::vector<VertexElement<2,2>*> elements_2d;
        elements_2d.push_back(new VertexElement<2,2>(0, nodes_elements_2d[0]));
        elements_2d.push_back(new VertexElement<2,2>(1, nodes_elements_2d[1]));

        VertexMesh<2,2> mesh_2d(nodes_2d, elements_2d);

        // Test the mesh has the correct numbers of nodes and elements
        TS_ASSERT_EQUALS(mesh_2d.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_2d.GetNumNodes(), 7u);

        // Further tests of the mesh
        TS_ASSERT_DELTA(mesh_2d.GetNode(2)->rGetLocation()[0], 1.5, 1e-3);
        TS_ASSERT_EQUALS(mesh_2d.GetElement(1)->GetNode(2)->GetIndex(), 6u);

        // Test that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);

        // Nodes 1 and 4 are only in element 0
        TS_ASSERT_EQUALS(nodes_2d[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes_2d[4]->rGetContainingElementIndices(), temp_list1);

        // Node 2 is in elements 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(nodes_2d[2]->rGetContainingElementIndices(), temp_list1);

        // Node 5 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(nodes_2d[5]->rGetContainingElementIndices(), temp_list2);

        // Coverage of some methods
        TS_ASSERT_EQUALS(mesh_2d.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(mesh_2d.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(mesh_2d.SolveBoundaryElementMapping(0), 0u);
    }

    void TestBasic3dVertexMesh()
    {
        // Create a 3D mesh comprising a cube and pyramid
        VertexMesh<3,3>* p_mesh = ConstructCubeAndPyramidMesh();

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

        // Coverage of some methods
        TS_ASSERT_EQUALS(p_mesh->SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveBoundaryElementMapping(0), 0u);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestMeshConstructionFromMeshReader()
    {
        // Test construction of a simple 2D mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));

        // Test construction of a mesh whose elements have attributes
        VertexMeshReader<2,2> mesh_reader2("mesh/test/data/TestVertexMesh/vertex_mesh_with_attributes");
        VertexMesh<2,2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);

        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 2u);

        TS_ASSERT_DELTA(mesh2.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh2.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh2.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNode(1), mesh2.GetNode(5));

        TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetUnsignedAttribute(), 76u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetUnsignedAttribute(), 89u);

        // Test the correct exception is thrown in other cases
        VertexMeshReader<1,1> mesh_reader_11("mesh/test/data/TestVertexMesh/vertex_mesh_with_attributes");
        VertexMesh<1,1> mesh_11;
        TS_ASSERT_THROWS_THIS(mesh_11.ConstructFromMeshReader(mesh_reader_11),
            "VertexMesh<1,1>::ConstructFromMeshReader() is not implemented");

        VertexMeshReader<1,2> mesh_reader_12("mesh/test/data/TestVertexMesh/vertex_mesh_with_attributes");
        VertexMesh<1,2> mesh_12;
        TS_ASSERT_THROWS_THIS(mesh_12.ConstructFromMeshReader(mesh_reader_12),
            "VertexMesh<1,2>::ConstructFromMeshReader() is not implemented");

        VertexMeshReader<1,3> mesh_reader_13("mesh/test/data/TestVertexMesh/vertex_mesh_with_attributes");
        VertexMesh<1,3> mesh_13;
        TS_ASSERT_THROWS_THIS(mesh_13.ConstructFromMeshReader(mesh_reader_13),
            "VertexMesh<1,3>::ConstructFromMeshReader() is not implemented");

        VertexMeshReader<2,3> mesh_reader_23("mesh/test/data/TestVertexMesh/vertex_mesh_with_attributes");
        VertexMesh<2,3> mesh_23;
        TS_ASSERT_THROWS_THIS(mesh_23.ConstructFromMeshReader(mesh_reader_23),
            "VertexMesh<2,3>::ConstructFromMeshReader() is not implemented");
    }

    void TestMeshConstructionFromMeshReaderIndexedFromOne()
    {
        // Test construction of a simple 2D mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh_elements_indexed_from_1");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
    }

    void TestMeshConstructionFromMeshReaderWithFaces()
    {
        // Test construction of a simple 3D mesh
        VertexMeshReader<3,3> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_3d_with_faces");
        VertexMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumFaces(), 4u);

        // Check Voronoi nodes are correct
        c_vector<double, 3> node_0_location = mesh.GetNode(0)->rGetLocation();
        TS_ASSERT_DELTA(node_0_location[0], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_0_location[1], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_0_location[2], -0.25, 1e-6);

        c_vector<double, 3> node_1_location = mesh.GetNode(1)->rGetLocation();
        TS_ASSERT_DELTA(node_1_location[0], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_1_location[1], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_1_location[2], 1.25, 1e-6);

        c_vector<double, 3> node_2_location = mesh.GetNode(2)->rGetLocation();
        TS_ASSERT_DELTA(node_2_location[0], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_2_location[1], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_2_location[2], 1.25, 1e-6);

        c_vector<double, 3> node_3_location = mesh.GetNode(3)->rGetLocation();
        TS_ASSERT_DELTA(node_3_location[0], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_3_location[1], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_3_location[2], -0.25, 1e-6);

        // Check Voronoi faces are correct
        VertexElement<2,3>* p_face_0 = mesh.GetFace(0);
        TS_ASSERT_EQUALS(p_face_0->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(mesh.CalculateAreaOfFace(p_face_0), 9.0*sqrt(3.0)/8.0, 1e-4);

        VertexElement<2,3>* p_face_1 = mesh.GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_DELTA(mesh.CalculateAreaOfFace(p_face_1), 9.0*sqrt(3.0)/8.0, 1e-4);

        VertexElement<2,3>* p_face_2 = mesh.GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(mesh.CalculateAreaOfFace(p_face_2), 9.0*sqrt(3.0)/8.0, 1e-4);

        VertexElement<2,3>* p_face_3 = mesh.GetFace(3);
        TS_ASSERT_EQUALS(p_face_3->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(mesh.CalculateAreaOfFace(p_face_3), 9.0*sqrt(3.0)/8.0, 1e-4);

        // Check Voronoi element is correct
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumFaces(), 4u);
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 1.125, 1e-4);
        TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(0), 9.0*sqrt(3.0)/2.0, 1e-4);

        // Create mesh in which elements have attributes
        VertexMeshReader<3,3> mesh_reader2("mesh/test/data/TestVertexMeshWriter/vertex_mesh_3d_with_attributes");
        VertexMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);

        // Check we have the right number of nodes, elements and faces
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumFaces(), 4u);

        // Check element attributes
        TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetAttribute(), 49u);
    }

    void TestArchive2dVertexMesh()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vertex_mesh_2d.arch";
        ArchiveLocationInfo::SetMeshFilename("vertex_mesh");

        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_5_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

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
            TS_ASSERT_EQUALS((static_cast<VertexMesh<2,2>*>(p_mesh))->GetNumNodes(), 46u);
            TS_ASSERT_EQUALS((static_cast<VertexMesh<2,2>*>(p_mesh))->GetNumElements(), 15u);

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

            VertexMesh<2,2>* p_mesh_original = static_cast<VertexMesh<2,2>*>(p_mesh);
            VertexMesh<2,2>* p_mesh_loaded = static_cast<VertexMesh<2,2>*>(p_mesh2);

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

            // Avoid memory leaks
            delete p_mesh_loaded;
        }
    }

    void TestArchive3dVertexMesh()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vertex_mesh_3d.arch";
        ArchiveLocationInfo::SetMeshFilename("vertex_mesh");

        AbstractMesh<3,3>* const p_mesh = ConstructCubeAndPyramidMesh();

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
            TS_ASSERT_EQUALS((static_cast<VertexMesh<3,3>*>(p_mesh))->GetNumNodes(), 9u);
            TS_ASSERT_EQUALS((static_cast<VertexMesh<3,3>*>(p_mesh))->GetNumElements(), 2u);
            TS_ASSERT_EQUALS((static_cast<VertexMesh<3,3>*>(p_mesh))->GetNumFaces(), 10u);

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

            VertexMesh<3,3>* p_mesh_original = static_cast<VertexMesh<3,3>*>(p_mesh);
            VertexMesh<3,3>* p_mesh_loaded = static_cast<VertexMesh<3,3>*>(p_mesh2);

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

            // Avoid memory leaks
            delete p_mesh_loaded;
        }

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestVertexElementMap()
    {
        VertexElementMap map(10);
        TS_ASSERT_EQUALS(map.Size(), 10u);

        map.ResetToIdentity();
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        map.SetNewIndex(0,1);
        map.SetNewIndex(1,0);

        TS_ASSERT_EQUALS(map.GetNewIndex(0), 1u);
        TS_ASSERT_EQUALS(map.GetNewIndex(1), 0u);
        TS_ASSERT_EQUALS(map.GetNewIndex(2), 2u);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);

        map.ResetToIdentity();
        map.SetDeleted(4);
        TS_ASSERT_THROWS_THIS(map.GetNewIndex(4), "VertexElement has been deleted");
        TS_ASSERT_EQUALS(map.IsDeleted(4), true);
        TS_ASSERT_EQUALS(map.IsDeleted(5), false);
        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);
    }

    void TestNeighbouringNodeAndElementMethods()
    {
        // Test methods with a small regular mesh comprising hexagonal elements
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_2_by_2");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 16u);

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(1), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(3), 8u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(4), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(5), 2u);

        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(0), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(1), 9u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(2), 12u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(3), 14u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(4), 11u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(5), 8u);

        // Check we have the correct neighbours for node 6
        std::set<unsigned> node_neighbours = mesh.GetNeighbouringNodeIndices(6);

        std::set<unsigned> expected_node_neighbours;
        expected_node_neighbours.insert(3);
        expected_node_neighbours.insert(8);
        expected_node_neighbours.insert(9);

        TS_ASSERT_EQUALS(node_neighbours, expected_node_neighbours);

        // Check that the only neighbour not also in element 2 is node 3
        std::set<unsigned> node_neighbours_not_in_elem2 = mesh.GetNeighbouringNodeNotAlsoInElement(6, 2);

        TS_ASSERT_EQUALS(node_neighbours_not_in_elem2.size(), 1u);
        TS_ASSERT_EQUALS(*(node_neighbours_not_in_elem2.begin()), 3u);

        // Check an exception is thrown if we use the index of a node not contained in this element
        TS_ASSERT_THROWS_THIS(mesh.GetNeighbouringNodeNotAlsoInElement(0, 2),
                              "The given node is not contained in the given element.");

        // Check element neighbours
        std::set<unsigned> element_neighbours = mesh.GetNeighbouringElementIndices(0);

        std::set<unsigned> expected_element_neighbours;
        expected_element_neighbours.insert(1);
        expected_element_neighbours.insert(2);

        TS_ASSERT_EQUALS(element_neighbours, expected_element_neighbours);
    }

    void TestGetRosetteRankOfElement()
    {
        // Test method on a honeycomb mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        // The rosette rank of each element should be 3
        for (unsigned elem_idx = 0 ; elem_idx < regular_mesh.GetNumElements() ; elem_idx++)
        {
            TS_ASSERT_EQUALS(regular_mesh.GetRosetteRankOfElement(elem_idx), 3u);
        }

        // Generate a five-element rosette
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.00000, 0.00000));
        nodes.push_back(new Node<2>(2, true, 0.80902, 0.58779));
        nodes.push_back(new Node<2>(3, true, 0.30902, 0.95106));
        nodes.push_back(new Node<2>(4, true, -0.30902, 0.95106));
        nodes.push_back(new Node<2>(5, true, -0.80902, 0.58779));
        nodes.push_back(new Node<2>(6, true, -1.00000, 0.00000));
        nodes.push_back(new Node<2>(7, true, -0.80902, -0.58779));
        nodes.push_back(new Node<2>(8, true, -0.30902, -0.95106));
        nodes.push_back(new Node<2>(9, true, 0.30902, -0.95106));
        nodes.push_back(new Node<2>(10, true, 0.80902, -0.58779));

        // Make 5 quadrangular cells
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);
        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[5]);
        nodes_elem_3.push_back(nodes[6]);
        nodes_elem_3.push_back(nodes[7]);
        std::vector<Node<2>*> nodes_elem_4;
        nodes_elem_4.push_back(nodes[0]);
        nodes_elem_4.push_back(nodes[7]);
        nodes_elem_4.push_back(nodes[8]);
        nodes_elem_4.push_back(nodes[9]);
        std::vector<Node<2>*> nodes_elem_5;
        nodes_elem_5.push_back(nodes[0]);
        nodes_elem_5.push_back(nodes[9]);
        nodes_elem_5.push_back(nodes[10]);
        nodes_elem_5.push_back(nodes[1]);

        // Make 5 vertex elements
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_4));
        vertex_elements.push_back(new VertexElement<2,2>(4, nodes_elem_5));

        VertexMesh<2,2> rosette_mesh(nodes, vertex_elements);

        // The rosette rank of each element should be five
        for (unsigned elem_idx = 0 ; elem_idx < rosette_mesh.GetNumElements() ; elem_idx++)
        {
            TS_ASSERT_EQUALS(rosette_mesh.GetRosetteRankOfElement(elem_idx), 5u);
        }
    }

    void TestGetCentroidOfElement()
    {
        // Test method with a 1D mesh
        std::vector<Node<1>*> nodes_1d;
        for (unsigned i=0; i<4; i++)
        {
            nodes_1d.push_back(new Node<1>(i, false, 0.5*(double)i));
        }
        std::vector<std::vector<Node<1>*> > nodes_elements_1d(3);
        std::vector<VertexElement<1,1>*> elements_1d;
        for (unsigned i=0; i<3; i++)
        {
            nodes_elements_1d[i].push_back(nodes_1d[i]);
            nodes_elements_1d[i].push_back(nodes_1d[i+1]);
            elements_1d.push_back(new VertexElement<1,1>(i, nodes_elements_1d[i]));
        }
        VertexMesh<1,1> mesh_1d(nodes_1d, elements_1d);

        TS_ASSERT_DELTA(mesh_1d.GetCentroidOfElement(0)[0], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetCentroidOfElement(1)[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetCentroidOfElement(2)[0], 1.25, 1e-6);

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> triangle_elements;
        triangle_elements.push_back(new VertexElement<2,2>(0, triangle_nodes));
        VertexMesh<2,2> triangle_mesh(triangle_nodes, triangle_elements);

        c_vector<double, 2> triangle_centroid = triangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(triangle_centroid[0], 2.0/3.0, 1e-4);
        TS_ASSERT_DELTA(triangle_centroid[1], 1.0/3.0, 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));
        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        c_vector<double, 2> square_centroid = square_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(square_centroid[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(square_centroid[1], 0.5, 1e-6);

        // Test method with a single rectangular element away from the origin
        std::vector<Node<2>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new Node<2>(0, false, 10.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(1, false, 11.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(2, false, 11.0, 14.0));
        far_rectangle_nodes.push_back(new Node<2>(3, false, 10.0, 14.0));
        std::vector<VertexElement<2,2>*> far_rectangle_elements;
        far_rectangle_elements.push_back(new VertexElement<2,2>(0, far_rectangle_nodes));
        VertexMesh<2,2> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_elements);

        c_vector<double, 2> far_rectangle_centroid = far_rectangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(far_rectangle_centroid[0], 10.5, 1e-4);
        TS_ASSERT_DELTA(far_rectangle_centroid[1], 12.0, 1e-4);

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<Node<2>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new Node<2>(0, false,  2.0*0.5*sqrt(3.0) - 1.0*0.5,  2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(1, false, -2.0*0.5*sqrt(3.0) - 1.0*0.5, -2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(2, false, -2.0*0.5*sqrt(3.0) + 1.0*0.5, -2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(3, false,  2.0*0.5*sqrt(3.0) + 1.0*0.5,  2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        std::vector<VertexElement<2,2>*> angled_rectangle_elements;
        angled_rectangle_elements.push_back(new VertexElement<2,2>(0, angled_rectangle_nodes));
        VertexMesh<2,2> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_elements);

        c_vector<double, 2> angled_rectangle_centroid = angled_rectangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(angled_rectangle_centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(angled_rectangle_centroid(1), 0.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> circle_elements;
        circle_elements.push_back(new VertexElement<2,2>(0, circle_nodes));
        VertexMesh<2,2> circle_mesh(circle_nodes, circle_elements);

        c_vector<double, 2> circle_centroid = circle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(circle_centroid[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(circle_centroid[1], 0.0, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> hexagon_elements;
        hexagon_elements.push_back(new VertexElement<2,2>(0, hexagon_nodes));
        VertexMesh<2,2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        c_vector<double, 2> hexagon_centroid = hexagon_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(hexagon_centroid[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(hexagon_centroid[1], 0.0, 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158,  1.5588));
        std::vector<VertexElement<2,2>*> irregular_elements;
        irregular_elements.push_back(new VertexElement<2,2>(0, irregular_nodes));
        VertexMesh<2,2> irregular_mesh(irregular_nodes, irregular_elements);

        c_vector<double, 2> irregular_centroid = irregular_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(irregular_centroid[0], 2.6269, 1e-4);
        TS_ASSERT_DELTA(irregular_centroid[1], 0.8930, 1e-4);

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        c_vector<double, 2> regular_centroid_5 = regular_mesh.GetCentroidOfElement(5);
        TS_ASSERT_DELTA(regular_centroid_5[0], 2.0,           1e-4);
        TS_ASSERT_DELTA(regular_centroid_5[1], 2.5/sqrt(3.0), 1e-4);

        c_vector<double, 2> regular_centroid_7 = regular_mesh.GetCentroidOfElement(7);
        TS_ASSERT_DELTA(regular_centroid_7[0], 4.0,           1e-4);
        TS_ASSERT_DELTA(regular_centroid_7[1], 2.5/sqrt(3.0), 1e-4);

        // Test method with a 3D mesh
        VertexMesh<3,3>* p_mesh = ConstructPrismMesh();

        // By symmetry, the centroid of the prism should lie in the plane x=0.5
        c_vector<double, 3> centroid = p_mesh->GetCentroidOfElement(0);
        TS_ASSERT_DELTA(centroid(0), 0.5, 1e-5);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestGetVolumeOfElement()
    {
        // Note that the method GetVolumeOfElement() only works in 2D and 3D

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> triangle_elements;
        triangle_elements.push_back(new VertexElement<2,2>(0, triangle_nodes));
        VertexMesh<2,2> triangle_mesh(triangle_nodes, triangle_elements);

        TS_ASSERT_DELTA(triangle_mesh.GetVolumeOfElement(0), 1.0, 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));
        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        TS_ASSERT_DELTA(square_mesh.GetVolumeOfElement(0), 1.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> circle_elements;
        circle_elements.push_back(new VertexElement<2,2>(0, circle_nodes));
        VertexMesh<2,2> circle_mesh(circle_nodes, circle_elements);

        TS_ASSERT_DELTA(circle_mesh.GetVolumeOfElement(0), M_PI, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> hexagon_elements;
        hexagon_elements.push_back(new VertexElement<2,2>(0, hexagon_nodes));
        VertexMesh<2,2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        TS_ASSERT_DELTA(hexagon_mesh.GetVolumeOfElement(0), 1.5*sqrt(3.0), 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158,  1.5588));
        std::vector<VertexElement<2,2>*> irregular_elements;
        irregular_elements.push_back(new VertexElement<2,2>(0, irregular_nodes));
        VertexMesh<2,2> irregular_mesh(irregular_nodes, irregular_elements);

        TS_ASSERT_DELTA(irregular_mesh.GetVolumeOfElement(0), 1.4684, 1e-3);

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        for (VertexMesh<2,2>::VertexElementIterator iter = regular_mesh.GetElementIteratorBegin();
             iter != regular_mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            TS_ASSERT_DELTA(regular_mesh.GetVolumeOfElement(elem_index), 0.5*sqrt(3.0), 1e-4);
        }

        // Test method with a 3D mesh
        VertexMesh<3,3>* p_mesh = ConstructPrismMesh();

        // The volume of the prism should be 0.5 * 3 * 2 * 1 = 3
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(0), 3.0, 1e-6);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestGetSurfaceAreaOfElement()
    {
        // Note that the method GetSurfaceAreaOfElement() only works in 2D and 3D

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> triangle_elements;
        triangle_elements.push_back(new VertexElement<2,2>(0, triangle_nodes));
        VertexMesh<2,2> triangle_mesh(triangle_nodes, triangle_elements);

        TS_ASSERT_DELTA(triangle_mesh.GetSurfaceAreaOfElement(0), 3.0 + sqrt(5.0), 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));
        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        TS_ASSERT_DELTA(square_mesh.GetSurfaceAreaOfElement(0), 4.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> circle_elements;
        circle_elements.push_back(new VertexElement<2,2>(0, circle_nodes));
        VertexMesh<2,2> circle_mesh(circle_nodes, circle_elements);

        TS_ASSERT_DELTA(circle_mesh.GetSurfaceAreaOfElement(0), 2.0*M_PI, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> hexagon_elements;
        hexagon_elements.push_back(new VertexElement<2,2>(0, hexagon_nodes));
        VertexMesh<2,2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        TS_ASSERT_DELTA(hexagon_mesh.GetSurfaceAreaOfElement(0), 6.0, 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158,  1.5588));
        std::vector<VertexElement<2,2>*> irregular_elements;
        irregular_elements.push_back(new VertexElement<2,2>(0, irregular_nodes));
        VertexMesh<2,2> irregular_mesh(irregular_nodes, irregular_elements);

        TS_ASSERT_DELTA(irregular_mesh.GetSurfaceAreaOfElement(0), 5.1263, 1e-3);

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        for (VertexMesh<2,2>::VertexElementIterator iter = regular_mesh.GetElementIteratorBegin();
             iter != regular_mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            TS_ASSERT_DELTA(regular_mesh.GetSurfaceAreaOfElement(elem_index), 2*sqrt(3.0), 1e-4);
        }

        // Test method with a 3D mesh
        VertexMesh<3,3>* p_mesh = ConstructPrismMesh();

        // The surface area of the prism should be the sum of the face areas
        TS_ASSERT_DELTA(p_mesh->GetSurfaceAreaOfElement(0), 11 + sqrt(13.0), 1e-6);

        // Avoid memory leaks
        delete p_mesh;
    }

    void Test3dMethodsWithPrism()
    {
        // Test method with a 3D mesh
        VertexMesh<3,3>* p_mesh = ConstructPrismMesh();

        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 5u);

        // Face 0 has four vertices, is perpendicular to the y axis, and has area 1*3 = 3
        VertexElement<2,3>* p_face_0 = p_mesh->GetFace(0);
        TS_ASSERT_EQUALS(p_face_0->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_0;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_0, unit_normal_0);
        TS_ASSERT_DELTA(unit_normal_0[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_0[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_0[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_0), 3.0, 1e-6);

        // Face 1 has three vertices, is perpendicular to the x axis, and has area 0.5*2*3 = 3
        VertexElement<2,3>* p_face_1 = p_mesh->GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 3u);
        c_vector<double, 3> unit_normal_1;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_1, unit_normal_1);
        TS_ASSERT_DELTA(unit_normal_1[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_1), 3.0, 1e-6);

        // Face 2 has four vertices, is at an angle theta to the y axis where tan(theta) = 2/3,
        // and has area 1*sqrt(2^2 + 3^2) = sqrt(13.0)
        VertexElement<2,3>* p_face_2 = p_mesh->GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_2;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_2, unit_normal_2);
        TS_ASSERT_DELTA(unit_normal_2[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[1], -sin(atan2(3.0,2.0)), 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[2], -cos(atan2(3.0,2.0)), 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_2), sqrt(13.0), 1e-6);

        // Face 1 has three vertices, is perpendicular to the x axis, and has area 0.5*2*3 = 3
        VertexElement<2,3>* p_face_3 = p_mesh->GetFace(3);
        TS_ASSERT_EQUALS(p_face_3->GetNumNodes(), 3u);
        c_vector<double, 3> unit_normal_3;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_3, unit_normal_3);
        TS_ASSERT_DELTA(unit_normal_3[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_3[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_3[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_3), 3.0, 1e-6);

        // Face 4 has four vertices, is perpendicular to the z axis, and has area 1*2 = 2
        VertexElement<2,3>* p_face_4 = p_mesh->GetFace(4);
        TS_ASSERT_EQUALS(p_face_4->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_4;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_4, unit_normal_4);
        TS_ASSERT_DELTA(unit_normal_4[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_4[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_4[2], -1.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_4), 2.0, 1e-6);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestGetAreaGradientOfElementAtNode()
    {
        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));

        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        VertexElement<2,2>* p_element = square_mesh.GetElement(0);

        c_vector<double, 2> element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 0);
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);

        element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 1);
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);

        element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 2);
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);

        element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 3);
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);
    }

    void TestGetPerimeterGradientAtNode()
    {
        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));

        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        VertexElement<2,2>* p_element = square_mesh.GetElement(0);

        c_vector<double, 2> element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 0);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 1);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 2);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);

        element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 3);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);
    }

    void TestMeshGetWidthAndBoundingBoxMethod()
    {
        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        // Test CalculateBoundingBox()
        ChasteCuboid<2> bounds = regular_mesh.CalculateBoundingBox();
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], 3.50,          1e-4);
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], 5.0/sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], 0.0,           1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], 0.0,           1e-4);

        // Test GetWidth()
        double width = regular_mesh.GetWidth(0);
        TS_ASSERT_DELTA(width, 3.5000, 1e-4);

        double height = regular_mesh.GetWidth(1);
        TS_ASSERT_DELTA(height, 5.0/sqrt(3.0), 1e-4);
    }

    void TestCalculateMomentsOfElement()
    {
        // Test method with a single triangular element
        std::vector<Node<2>*> isos_triangle_nodes;
        isos_triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        isos_triangle_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        isos_triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> isos_triangle_elements;
        isos_triangle_elements.push_back(new VertexElement<2,2>(0, isos_triangle_nodes));
        VertexMesh<2,2> isos_triangle_mesh(isos_triangle_nodes, isos_triangle_elements);

        c_vector<double, 3> isos_triangle_moments = isos_triangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(isos_triangle_moments[0], 0.0277, 1e-4);
        TS_ASSERT_DELTA(isos_triangle_moments[1], 0.0277, 1e-4);
        TS_ASSERT_DELTA(isos_triangle_moments[2], -0.0138, 1e-4);

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> triangle_elements;
        triangle_elements.push_back(new VertexElement<2,2>(0, triangle_nodes));
        VertexMesh<2,2> triangle_mesh(triangle_nodes, triangle_elements);

        c_vector<double, 3> triangle_moments = triangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(triangle_moments(0), 1.0/18.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(triangle_moments(1), 2.0/9.0, 1e-6);   // Iyy
        TS_ASSERT_DELTA(triangle_moments(2), -5.0/90.0, 1e-6); // Ixy

        // Test method with a single rectangular element parallel to the x-axis
        std::vector<Node<2>*> horizontal_rectangle_nodes;
        horizontal_rectangle_nodes.push_back(new Node<2>(0, false,  2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new Node<2>(1, false, -2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new Node<2>(2, false, -2.0, -1.0));
        horizontal_rectangle_nodes.push_back(new Node<2>(3, false,  2.0, -1.0));
        std::vector<VertexElement<2,2>*> horizontal_rectangle_elements;
        horizontal_rectangle_elements.push_back(new VertexElement<2,2>(0, horizontal_rectangle_nodes));
        VertexMesh<2,2> horizontal_rectangle_mesh(horizontal_rectangle_nodes, horizontal_rectangle_elements);

        c_vector<double, 3> horizontal_rectangle_moments = horizontal_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(horizontal_rectangle_moments(0), 8.0/3.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(horizontal_rectangle_moments(1), 32.0/3.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(horizontal_rectangle_moments(2), 0.0, 1e-6);      // Ixy = 0 by symmetry

        // Test method with the same shape, but supply nodes in clockwise manner
        std::vector<Node<2>*> clockwise_rectangle_nodes;
        clockwise_rectangle_nodes.push_back(new Node<2>(0, false, -2.0, -1.0));
        clockwise_rectangle_nodes.push_back(new Node<2>(1, false, -2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new Node<2>(2, false,  2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new Node<2>(3, false,  2.0, -1.0));
        std::vector<VertexElement<2,2>*> clockwise_rectangle_elements;
        clockwise_rectangle_elements.push_back(new VertexElement<2,2>(0, clockwise_rectangle_nodes));
        VertexMesh<2,2> clockwise_rectangle_mesh(clockwise_rectangle_nodes, clockwise_rectangle_elements);

        c_vector<double, 3> clockwise_rectangle_moments = clockwise_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(clockwise_rectangle_moments(0), 8.0/3.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(clockwise_rectangle_moments(1), 32.0/3.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(clockwise_rectangle_moments(2), 0.0, 1e-6);      // Ixy = 0 by symmetry

        // Test method with a single rectangular element parallel to the y-axis
        std::vector<Node<2>*> vertical_rectangle_nodes;
        vertical_rectangle_nodes.push_back(new Node<2>(0, false,  1.0,  2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(1, false, -1.0,  2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(2, false, -1.0, -2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(3, false,  1.0, -2.0));
        std::vector<VertexElement<2,2>*> vertical_rectangle_elements;
        vertical_rectangle_elements.push_back(new VertexElement<2,2>(0, vertical_rectangle_nodes));
        VertexMesh<2,2> vertical_rectangle_mesh(vertical_rectangle_nodes, vertical_rectangle_elements);

        c_vector<double, 3> vertical_rectangle_moments = vertical_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_moments(0), 32.0/3.0, 1e-6); // Ixx
        TS_ASSERT_DELTA(vertical_rectangle_moments(1), 8.0/3.0, 1e-6);  // Iyy
        TS_ASSERT_DELTA(vertical_rectangle_moments(2), 0.0, 1e-6);      // Ixy = 0 by symmetry

        // Test method with a single rectangular element away from the origin
        std::vector<Node<2>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new Node<2>(0, false, 10.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(1, false, 11.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(2, false, 11.0, 14.0));
        far_rectangle_nodes.push_back(new Node<2>(3, false, 10.0, 14.0));
        std::vector<VertexElement<2,2>*> far_rectangle_elements;
        far_rectangle_elements.push_back(new VertexElement<2,2>(0, far_rectangle_nodes));
        VertexMesh<2,2> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_elements);

        c_vector<double, 3> far_rectangle_moments = far_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(far_rectangle_moments[0], 16.0/3.0, 1e-4); // Ixx
        TS_ASSERT_DELTA(far_rectangle_moments[1], 1.0/3.0, 1e-4);  // Iyy
        TS_ASSERT_DELTA(far_rectangle_moments[2], 0.0, 1e-4);      // Ixy = 0 by symmetry

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<Node<2>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new Node<2>(0, false,  2.0*0.5*sqrt(3.0) - 1.0*0.5,  2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(1, false, -2.0*0.5*sqrt(3.0) - 1.0*0.5, -2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(2, false, -2.0*0.5*sqrt(3.0) + 1.0*0.5, -2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(3, false,  2.0*0.5*sqrt(3.0) + 1.0*0.5,  2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        std::vector<VertexElement<2,2>*> angled_rectangle_elements;
        angled_rectangle_elements.push_back(new VertexElement<2,2>(0, angled_rectangle_nodes));
        VertexMesh<2,2> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_elements);

        c_vector<double, 3> angled_rectangle_moments = angled_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(angled_rectangle_moments[0], 14.0/3.0, 1e-4);    // Ixx
        TS_ASSERT_DELTA(angled_rectangle_moments[1], 26.0/3.0, 1e-4);    // Iyy
        TS_ASSERT_DELTA(angled_rectangle_moments[2], 2.0*sqrt(3.0), 1e-4); // Ixy

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158,  1.5588));
        std::vector<VertexElement<2,2>*> irregular_elements;
        irregular_elements.push_back(new VertexElement<2,2>(0, irregular_nodes));
        VertexMesh<2,2> irregular_mesh(irregular_nodes, irregular_elements);

        c_vector<double, 3> irregular_moments = irregular_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(irregular_moments[0], 0.3521, 1e-4); // Ixx
        TS_ASSERT_DELTA(irregular_moments[1], 0.1264, 1e-4); // Iyy
        TS_ASSERT_DELTA(irregular_moments[2], 0.1162, 1e-4); // Ixy

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        for (unsigned i=0; i<regular_mesh.GetNumElements(); i++)
        {
            c_vector<double, 3> regular_moments = regular_mesh.CalculateMomentsOfElement(i);
            TS_ASSERT_DELTA(regular_moments(0), 5*sqrt(3.0)/16/9, 1e-6); // Ixx
            TS_ASSERT_DELTA(regular_moments(1), 5*sqrt(3.0)/16/9, 1e-6); // Iyy
            TS_ASSERT_DELTA(regular_moments(2), 0.0,            1e-6); // Ixy = 0 by symmetry
        }
    }

    void TestGetShortAxisOfElement()
    {
        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> triangle_elements;
        triangle_elements.push_back(new VertexElement<2,2>(0, triangle_nodes));
        VertexMesh<2,2> triangle_mesh(triangle_nodes, triangle_elements);

        c_vector<double, 2> triangle_short_axis = triangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(triangle_short_axis[0], 1/sqrt(2.0), 1e-4);
        TS_ASSERT_DELTA(triangle_short_axis[1], 1/sqrt(2.0), 1e-4);

        // Test method with a single rectangular element parallel to the x-axis
        std::vector<Node<2>*> horizontal_rectangle_nodes;
        horizontal_rectangle_nodes.push_back(new Node<2>(0, false,  2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new Node<2>(1, false, -2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new Node<2>(2, false, -2.0, -1.0));
        horizontal_rectangle_nodes.push_back(new Node<2>(3, false,  2.0, -1.0));
        std::vector<VertexElement<2,2>*> horizontal_rectangle_elements;
        horizontal_rectangle_elements.push_back(new VertexElement<2,2>(0, horizontal_rectangle_nodes));
        VertexMesh<2,2> horizontal_rectangle_mesh(horizontal_rectangle_nodes, horizontal_rectangle_elements);

        c_vector<double, 2> horizontal_rectangle_short_axis = horizontal_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(horizontal_rectangle_short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(horizontal_rectangle_short_axis(1), 1.0, 1e-6);

        // Test method with the same shape, but supply nodes in clockwise manner
        std::vector<Node<2>*> clockwise_rectangle_nodes;
        clockwise_rectangle_nodes.push_back(new Node<2>(0, false, -2.0, -1.0));
        clockwise_rectangle_nodes.push_back(new Node<2>(1, false, -2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new Node<2>(2, false,  2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new Node<2>(3, false,  2.0, -1.0));
        std::vector<VertexElement<2,2>*> clockwise_rectangle_elements;
        clockwise_rectangle_elements.push_back(new VertexElement<2,2>(0, clockwise_rectangle_nodes));
        VertexMesh<2,2> clockwise_rectangle_mesh(clockwise_rectangle_nodes, clockwise_rectangle_elements);

        c_vector<double, 3> clockwise_rectangle_short_axis = clockwise_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(clockwise_rectangle_short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(clockwise_rectangle_short_axis(1), 1.0, 1e-6);

        // Test method with a single rectangular element parallel to the y-axis
        std::vector<Node<2>*> vertical_rectangle_nodes;
        vertical_rectangle_nodes.push_back(new Node<2>(0, false,  1.0,  2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(1, false, -1.0,  2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(2, false, -1.0, -2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(3, false,  1.0, -2.0));
        std::vector<VertexElement<2,2>*> vertical_rectangle_elements;
        vertical_rectangle_elements.push_back(new VertexElement<2,2>(0, vertical_rectangle_nodes));
        VertexMesh<2,2> vertical_rectangle_mesh(vertical_rectangle_nodes, vertical_rectangle_elements);

        c_vector<double, 2> vertical_rectangle_short_axis = vertical_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_short_axis(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(vertical_rectangle_short_axis(1), 0.0, 1e-6);

        // Test method with a single rectangular element parallel to the y-axis, centred away from the origin
        std::vector<Node<2>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new Node<2>(0, false, 10.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(1, false, 11.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(2, false, 11.0, 14.0));
        far_rectangle_nodes.push_back(new Node<2>(3, false, 10.0, 14.0));
        std::vector<VertexElement<2,2>*> far_rectangle_elements;
        far_rectangle_elements.push_back(new VertexElement<2,2>(0, far_rectangle_nodes));
        VertexMesh<2,2> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_elements);

        c_vector<double, 2> far_rectangle_short_axis = far_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(far_rectangle_short_axis[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(far_rectangle_short_axis[1], 0.0, 1e-4);

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<Node<2>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new Node<2>(0, false,  2.0*0.5*sqrt(3.0) - 1.0*0.5,  2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(1, false, -2.0*0.5*sqrt(3.0) - 1.0*0.5, -2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(2, false, -2.0*0.5*sqrt(3.0) + 1.0*0.5, -2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(3, false,  2.0*0.5*sqrt(3.0) + 1.0*0.5,  2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        std::vector<VertexElement<2,2>*> angled_rectangle_elements;
        angled_rectangle_elements.push_back(new VertexElement<2,2>(0, angled_rectangle_nodes));
        VertexMesh<2,2> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_elements);

        c_vector<double, 2> angled_rectangle_short_axis = angled_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(angled_rectangle_short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(angled_rectangle_short_axis(1), -0.5*sqrt(3.0), 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158,  1.5588));
        std::vector<VertexElement<2,2>*> irregular_elements;
        irregular_elements.push_back(new VertexElement<2,2>(0, irregular_nodes));
        VertexMesh<2,2> irregular_mesh(irregular_nodes, irregular_elements);

        c_vector<double, 2> irregular_short_axis = irregular_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(irregular_short_axis[0], 0.9210, 1e-4);
        TS_ASSERT_DELTA(irregular_short_axis[1], -0.3894, 1e-4);

        // Test with a single trapezoidal element of width 1, lengths 3*sqrt(3.0) and sqrt(3.0), rotated 30 degrees anticlockwise
        std::vector<Node<2>*> trapezium_nodes;
        trapezium_nodes.push_back(new Node<2>(0, false,  1.0, 0.0));
        trapezium_nodes.push_back(new Node<2>(1, false,  2.0, sqrt(3.0)));
        trapezium_nodes.push_back(new Node<2>(2, false, -2.5, -sqrt(3.0)/2.0));
        trapezium_nodes.push_back(new Node<2>(3, false, -0.5, -sqrt(3.0)/2.0));
        std::vector<VertexElement<2,2>*> trapezium_elements;
        trapezium_elements.push_back(new VertexElement<2,2>(0, trapezium_nodes));
        VertexMesh<2,2> trapezium_mesh(trapezium_nodes, trapezium_elements);

        c_vector<double, 2> trapezium_short_axis = trapezium_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(trapezium_short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(trapezium_short_axis(1), -0.5*sqrt(3.0), 1e-6);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> hexagon_elements;
        hexagon_elements.push_back(new VertexElement<2,2>(0, hexagon_nodes));
        VertexMesh<2,2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        // Since the element is a regular polygon, the short axis is a random vector, so test the random seed
        c_vector<double, 2> hexagon_short_axis = hexagon_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(norm_2(hexagon_short_axis), 1.0, 1e-6);
        TS_ASSERT_DELTA(hexagon_short_axis(0), 0.5488, 1e-4);
        TS_ASSERT_DELTA(hexagon_short_axis(1), 0.8359, 1e-4);
    }

    void TestGetElongationShapeFactorOfElement()
    {
        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false,  1.0,  1.0));
        square_nodes.push_back(new Node<2>(1, false, -1.0,  1.0));
        square_nodes.push_back(new Node<2>(2, false, -1.0, -1.0));
        square_nodes.push_back(new Node<2>(3, false,  1.0, -1.0));
        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));
        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        // By symmetry, the two eigenvalues are equal, so the elongation shape factor equals one
        double square_elongation_shape_factor = square_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(square_elongation_shape_factor, 1.0, 1e-6);

        // We should obtain the same result regardless of the square's size
        square_mesh.Scale(5.0, 5.0);
        square_elongation_shape_factor = square_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(square_elongation_shape_factor, 1.0, 1e-6);

        // Test method with a single rectangular element
        std::vector<Node<2>*> vertical_rectangle_nodes;
        vertical_rectangle_nodes.push_back(new Node<2>(0, false,  1.0,  2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(1, false, -1.0,  2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(2, false, -1.0, -2.0));
        vertical_rectangle_nodes.push_back(new Node<2>(3, false,  1.0, -2.0));
        std::vector<VertexElement<2,2>*> vertical_rectangle_elements;
        vertical_rectangle_elements.push_back(new VertexElement<2,2>(0, vertical_rectangle_nodes));
        VertexMesh<2,2> vertical_rectangle_mesh(vertical_rectangle_nodes, vertical_rectangle_elements);

        /*
         * For a rectangle of width a (parallel to the x axis) and height b (parallel to the
         * y axis), the second moments of area are given by J_xx = a*b^3/12, J_yy = b*a^3/12
         * and J_xy = 0. Therefore the two eigenvalues are given by J_xx and J_yy and if b>a,
         * the elongation shape factor should be equal to sqrt(a*b^3/b*a^3) = b/a.
         */
        double vertical_rectangle_elongation_shape_factor = vertical_rectangle_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_elongation_shape_factor, 4.0/2.0, 1e-6);

        // We should obtain the same result regardless of the orientation of the rectangle
        vertical_rectangle_mesh.Scale(2.0, 0.5);
        vertical_rectangle_elongation_shape_factor = vertical_rectangle_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_elongation_shape_factor, 4.0/2.0, 1e-6);

        // Check the formula is correct for a different aspect ratio
        vertical_rectangle_mesh.Scale(3.0, 1.0);
        vertical_rectangle_elongation_shape_factor = vertical_rectangle_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_elongation_shape_factor, 12.0/2.0, 1e-6);
    }

    void TestScaleAndTranslate()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA(mesh.GetWidth(0), 3.5000, 1e-4);
        TS_ASSERT_DELTA(mesh.GetWidth(1), 2.8867, 1e-4);

        // Squash in the x direction by a factor of 2
        mesh.Scale(0.5);
        TS_ASSERT_DELTA(mesh.GetWidth(0), 1.7500, 1e-4);
        TS_ASSERT_DELTA(mesh.GetWidth(1), 2.8867, 1e-4);

        // Stretch in the x and y directions by a factor of 2
        mesh.Scale(2.0, 2.0);
        TS_ASSERT_DELTA(mesh.GetWidth(0), 3.5000, 1e-4);
        TS_ASSERT_DELTA(mesh.GetWidth(1), 5.7735, 1e-4);

        // Create 3D mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 1.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 2.0, 3.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 3.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 2.0, 3.0));
        nodes.push_back(new Node<3>(7, false, 0.0, 2.0, 3.0));

        std::vector<VertexElement<3,3>*> elements;
        elements.push_back(new VertexElement<3,3>(0, nodes));

        VertexMesh<3,3> mesh3d(nodes, elements);

        TS_ASSERT_DELTA(mesh3d.GetWidth(0), 1.0, 1e-4);
        TS_ASSERT_DELTA(mesh3d.GetWidth(1), 2.0, 1e-4);
        TS_ASSERT_DELTA(mesh3d.GetWidth(2), 3.0, 1e-4);

        // Stretch the mesh
        mesh3d.Scale(4.0, 2.0, 4.0/3.0);

        TS_ASSERT_DELTA(mesh3d.GetWidth(0), 4.0, 1e-4);
        TS_ASSERT_DELTA(mesh3d.GetWidth(1), 4.0, 1e-4);
        TS_ASSERT_DELTA(mesh3d.GetWidth(2), 4.0, 1e-4);

        // Pick a certain node and store spatial position
        Node<3>* p_node = mesh3d.GetNode(7);
        ChastePoint<3> original_coordinate = p_node->GetPoint();

        const double x_movement = 1.0;
        const double y_movement = 2.5;
        const double z_movement = 2.5;

        mesh3d.Translate(x_movement, y_movement, z_movement);

        // Test the translate method
        ChastePoint<3>  new_coordinate = p_node->GetPoint();
        TS_ASSERT_DELTA(original_coordinate[0], new_coordinate[0] - x_movement, 1e-6);
        TS_ASSERT_DELTA(original_coordinate[1], new_coordinate[1] - y_movement, 1e-6);
        TS_ASSERT_DELTA(original_coordinate[2], new_coordinate[2] - z_movement, 1e-6);
    }

    void TestBoundaryNodes()
    {
        // Test with a single square element with all boundary nodes
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));
        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        for (unsigned i=0; i<square_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(square_mesh.GetNode(i)->IsBoundaryNode(), false);
        }

        // Test with a small regular honeycomb mesh with some interior nodes
        VertexMeshReader<2,2> regular_mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_2_by_2");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(regular_mesh_reader);

        for (unsigned i=0; i<regular_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = (i==6 || i==9) ? false : true;
            TS_ASSERT_EQUALS(regular_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

        // Test with a larger regular honeycomb mesh with some interior nodes
        VertexMeshReader<2,2> larger_regular_mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> larger_regular_mesh;
        larger_regular_mesh.ConstructFromMeshReader(larger_regular_mesh_reader);

        for (unsigned i=0; i<larger_regular_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==8 || i==9 || i==12 || i==13 || i==16 || i==17 || i==20 || i==21)
            {
                expected_boundary_node = false;
            }

            TS_ASSERT_EQUALS(larger_regular_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestTranslation2DWithUblas()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        c_vector<double, 2> old_location1 = mesh.GetNode(4)->rGetLocation();
        c_vector<double, 2> old_location2 = mesh.GetNode(9)->rGetLocation();

        // Set translation vector
        c_vector<double, 2> trans_vec;
        trans_vec(0) = 2.0;
        trans_vec(1) = 3.0;

        // Translate
        mesh.Translate(trans_vec);
        c_vector<double, 2> new_location1 = mesh.GetNode(4)->rGetLocation();
        c_vector<double, 2> new_location2 = mesh.GetNode(9)->rGetLocation();

        // Spot check a couple of nodes
        TS_ASSERT_DELTA(new_location1[0], old_location1[0] + 2.0, 1e-6);
        TS_ASSERT_DELTA(new_location1[1], old_location1[1] + 3.0, 1e-6);
        TS_ASSERT_DELTA(new_location2[0], old_location2[0] + 2.0, 1e-6);
        TS_ASSERT_DELTA(new_location2[1], old_location2[1] + 3.0, 1e-6);
    }

    void TestTranslation2DMethod()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Pick a random node and store spatial position
        Node<2>* p_node = mesh.GetNode(10);
        ChastePoint<2> original_coordinate = p_node->GetPoint();

        const double x_movement = 1.0;
        const double y_movement = 2.5;
        mesh.Translate(x_movement, y_movement);

        ChastePoint<2> new_coordinate = p_node->GetPoint();

        TS_ASSERT_DELTA(original_coordinate[0], new_coordinate[0] - x_movement, 1e-6);
        TS_ASSERT_DELTA(original_coordinate[1], new_coordinate[1] - y_movement, 1e-6);
    }

    void TestGenerateVerticesFromElementCircumcentres()
    {
        // Create a simple 3D tetrahedral mesh, the Delaunay triangulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  1.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(1, true, -1.0, -1.0,  1.0));
        nodes.push_back(new Node<3>(2, true, -1.0,  1.0, -1.0));
        nodes.push_back(new Node<3>(3, true,  1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(4, false, 0.0,  0.0,  0.0));
        MutableMesh<3,3> delaunay_mesh(nodes);

        /*
         * The Voronoi tessellation is not unique for this mesh, since four points are co-spherical.
         * We need to check how the mesher is breaking ties.
         */
        Element<3,3>* p_element = delaunay_mesh.GetElement(0);

        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 4u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 0u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3), 2u);

        // Create Voronoi tessellation mesh
        VertexMesh<3,3> tessellation(delaunay_mesh);
        tessellation.GenerateVerticesFromElementCircumcentres(delaunay_mesh);

        TS_ASSERT_EQUALS(tessellation.GetNumNodes(), 8u);

        c_vector<double,3> this_vertex = tessellation.GetNode(2)->rGetLocation();
        TS_ASSERT_DELTA(this_vertex[0],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);

        this_vertex = tessellation.GetNode(3)->rGetLocation();
        TS_ASSERT_DELTA(this_vertex[0],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2],  1.5, 1e-7);

        this_vertex = tessellation.GetNode(0)->rGetLocation();
        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2],  1.5, 1e-7);

        this_vertex = tessellation.GetNode(1)->rGetLocation();
        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);
    }

    void TestTessellationConstructor2d()
    {
        // Create a simple 2D tetrahedral mesh, the Delaunay triangulation
        std::vector<Node<2> *> delaunay_nodes;
        delaunay_nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        delaunay_nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        delaunay_nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        delaunay_nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        delaunay_nodes.push_back(new Node<2>(4, false, 0.5, 0.5));
        MutableMesh<2,2> delaunay_mesh(delaunay_nodes);

        TS_ASSERT_EQUALS(delaunay_mesh.CheckIsVoronoi(), true);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumNodes(), 5u);

        // Create a vertex mesh, the Voronoi tessellation, using the tetrahedral mesh
        VertexMesh<2,2> voronoi_mesh(delaunay_mesh);

        // Test the Voronoi tessellation has the correct number of nodes and elements
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumNodes(), 4u);

        // Test the location of the Voronoi nodes
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(0)->rGetLocation()[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(1)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(1)->rGetLocation()[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(2)->rGetLocation()[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(2)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(3)->rGetLocation()[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(3)->rGetLocation()[1], 1.0, 1e-6);

        // Test the number of nodes owned by each Voronoi element
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(0)->GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(1)->GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(2)->GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(3)->GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(4)->GetNumNodes(), 4u);

        // Test element areas
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(1), 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(2), 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(3), 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(4), 0.5, 1e-6);
    }

    void TestTessellationConstructor2dWithGhostNodes()
    {
        // Create a simple 2D tetrahedral mesh, the Delaunay triangulation
        std::vector<Node<2> *> delaunay_nodes;
        delaunay_nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        delaunay_nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        delaunay_nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        delaunay_nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        delaunay_nodes.push_back(new Node<2>(4, false, 0.5, 0.5));
        MutableMesh<2,2> delaunay_mesh(delaunay_nodes);

        TS_ASSERT_EQUALS(delaunay_mesh.CheckIsVoronoi(), true);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumNodes(), 5u);

        // Create a vertex mesh, the Voronoi tessellation, using the tetrahedral mesh and the 'real' nodes
        VertexMesh<2,2> voronoi_mesh(delaunay_mesh);

        // Test the Voronoi tessellation has the correct number of nodes and elements
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumNodes(), 4u);

        // Test the location of the Voronoi nodes
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(0)->rGetLocation()[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(1)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(1)->rGetLocation()[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(2)->rGetLocation()[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(2)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(3)->rGetLocation()[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(3)->rGetLocation()[1], 1.0, 1e-6);

        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(4)->GetNumNodes(), 4u);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(4), 0.5, 1e-6);
    }

    void TestGetEdgeLengthWithSimpleMesh()
    {
        // Create a simple 2D tetrahedral mesh, the Delaunay triangulation
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0.0,  0));
        nodes.push_back(new Node<2>(0, true,  0.0,  1));
        nodes.push_back(new Node<2>(0, true, -1.0,  0));
        nodes.push_back(new Node<2>(0, true,  1.0,  0));
        nodes.push_back(new Node<2>(0, true,  0.5, -0.5*sqrt(3.0)));
        nodes.push_back(new Node<2>(0, true, -0.5, -0.5*sqrt(3.0)));
        MutableMesh<2,2> delaunay_mesh(nodes);

        TS_ASSERT_EQUALS(delaunay_mesh.CheckIsVoronoi(), true);

        // Create Voronoi tessellation
        VertexMesh<2,2> voronoi_mesh(delaunay_mesh);

        TS_ASSERT_EQUALS(voronoi_mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumNodes(), 5u);

        // Measure the length of the edge separating the centre element and each of its neighbours
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,1), 1.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,2), 0.5 + 1.0/(sqrt(3.0)*2.0), 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,3), 0.5 + 1.0/(sqrt(3.0)*2.0), 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,4), 1.0/sqrt(3.0), 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,5), 1.0/sqrt(3.0), 1e-6);

        // All other neighbouring elements share an infinite edge
        TS_ASSERT_THROWS_THIS(voronoi_mesh.GetEdgeLength(1,2), "Elements 1 and 2 share only one node.");

        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(0), 0.5 + 0.25*sqrt(3.0), 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetSurfaceAreaOfElement(0), 2.0 + sqrt(3.0), 1e-6);
    }

    void TestTessellationConstructor3dWithGhostNode()
    {
        // Create a simple 3D tetrahedral mesh, the Delaunay triangulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> delaunay_mesh(nodes);

        TS_ASSERT_EQUALS(delaunay_mesh.CheckIsVoronoi(), true);

        VertexMesh<3,3> voronoi_mesh(delaunay_mesh);

        // Check there are as many nodes in the Voronoi mesh as there are elements in the Delaunay mesh
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumNodes(), 4u);

        // Check there are as many elements in the Voronoi mesh as there are non-boundary nodes in the Delaunay mesh
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumElements(), 1u);

        // Check there are as many faces in the Voronoi mesh as there are boundary nodes in the Delaunay mesh
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumFaces(), 4u);

        // Check Voronoi nodes are correct
        c_vector<double, 3> node_0_location = voronoi_mesh.GetNode(0)->rGetLocation();
        TS_ASSERT_DELTA(node_0_location[0], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_0_location[1], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_0_location[2], -0.25, 1e-6);

        c_vector<double, 3> node_1_location = voronoi_mesh.GetNode(1)->rGetLocation();
        TS_ASSERT_DELTA(node_1_location[0], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_1_location[1], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_1_location[2], 1.25, 1e-6);

        c_vector<double, 3> node_2_location = voronoi_mesh.GetNode(2)->rGetLocation();
        TS_ASSERT_DELTA(node_2_location[0], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_2_location[1], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_2_location[2], 1.25, 1e-6);

        c_vector<double, 3> node_3_location = voronoi_mesh.GetNode(3)->rGetLocation();
        TS_ASSERT_DELTA(node_3_location[0], -0.25, 1e-6);
        TS_ASSERT_DELTA(node_3_location[1], 1.25, 1e-6);
        TS_ASSERT_DELTA(node_3_location[2], -0.25, 1e-6);

        /* The Voronoi element is a tetrahedron formed by
         * (5/4, -1/4, -1/4) (-1/4, 5/4, -1/4) (-1/4, -1/4, 5/4) (5/4, 5/4, 5/4)
         * This is equivalent to (1.5, 0, 0), (0, 0, 1.5) etc.
         */

        // Check Voronoi faces are correct
        VertexElement<2,3>* p_face_0 = voronoi_mesh.GetFace(0);
        TS_ASSERT_EQUALS(p_face_0->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(p_face_0), 9.0*sqrt(3.0)/8.0, 1e-4);

        VertexElement<2,3>* p_face_1 = voronoi_mesh.GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(p_face_1), 9.0*sqrt(3.0)/8.0, 1e-4);

        VertexElement<2,3>* p_face_2 = voronoi_mesh.GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(p_face_2), 9.0*sqrt(3.0)/8.0, 1e-4);

        VertexElement<2,3>* p_face_3 = voronoi_mesh.GetFace(3);
        TS_ASSERT_EQUALS(p_face_3->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(p_face_3), 9.0*sqrt(3.0)/8.0, 1e-4);

        // Check Voronoi element is correct
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(0)->GetNumFaces(), 4u);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(0), 1.125, 1e-4);
        TS_ASSERT_DELTA(voronoi_mesh.GetSurfaceAreaOfElement(0), 9.0*sqrt(3.0)/2.0, 1e-4);
    }

    void TestTessellationConstructor3dWithRepeatedCircumcentres()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0,  false,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1,  false,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2,  false,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3,  false,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4,  false,  0.5, 0.5, 0.5));
        nodes.push_back(new Node<3>(5,  true,  -1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(6,  true,   2.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(7,  true,   2.0,  2.0, -1.0));
        nodes.push_back(new Node<3>(8,  true,  -1.0,  2.0, -1.0));
        nodes.push_back(new Node<3>(9,  true,  -1.0, -1.0,  2.0));
        nodes.push_back(new Node<3>(10, true,   2.0, -1.0,  2.0));
        nodes.push_back(new Node<3>(11, true,   2.0,  2.0,  2.0));
        nodes.push_back(new Node<3>(12, true,  -1.0,  2.0,  2.0));
        MutableMesh<3,3> delaunay_mesh(nodes);

        TS_ASSERT_EQUALS(delaunay_mesh.CheckIsVoronoi(), true);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumElements(), 32u);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumNodes(), 13u);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumBoundaryNodes(), 8u);
        TS_ASSERT_EQUALS(delaunay_mesh.GetNumNodes() - delaunay_mesh.GetNumBoundaryNodes(), 5u);

        VertexMesh<3,3> voronoi_mesh(delaunay_mesh);

        // Check there are as many nodes in the Voronoi mesh as there are elements in the Delaunay mesh
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumNodes(), 32u);

        // Note that some Voronoi nodes are repeated - most of the points above are co-spherical
        // In particular, nodes 0, 11 and 14 are all at (0.5, 0.5, -2.5)
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(0u)->rGetLocation()[0],   0.5, 1e-8);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(0u)->rGetLocation()[1],   0.5, 1e-8);
        TS_ASSERT_DELTA(voronoi_mesh.GetNode(0u)->rGetLocation()[2],   -2.5, 1e-8);
        TS_ASSERT_DELTA(norm_1(voronoi_mesh.GetNode(0u)->rGetLocation() - voronoi_mesh.GetNode(11u)->rGetLocation()),   0.0, DBL_EPSILON);
        TS_ASSERT_DELTA(norm_1(voronoi_mesh.GetNode(0u)->rGetLocation() - voronoi_mesh.GetNode(14u)->rGetLocation()),   0.0, DBL_EPSILON);

        TS_ASSERT_EQUALS(voronoi_mesh.GetNumFaces(), 32u);
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(voronoi_mesh.GetFace(0u)), 2.7556, 1e-4); //Degenerate quad (is triangle)
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(voronoi_mesh.GetFace(1u)), 2.7556, 1e-4); //Five point, but is triangle
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(voronoi_mesh.GetFace(2u)), 2.3864, 1e-4); //Six point, but is triangle
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(voronoi_mesh.GetFace(4u)), 0.0000, 1e-12); //Degenerate triangle
        TS_ASSERT_DELTA(voronoi_mesh.CalculateAreaOfFace(voronoi_mesh.GetFace(5u)), 2.7556, 1e-4); //Triangle
        for (unsigned i=0; i<15; i++)
        {
            voronoi_mesh.CalculateAreaOfFace(voronoi_mesh.GetFace(i));
        }
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumElements(), 5u);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(0u), 7.5937, 1e-4);
        double volume=0.0;
        for (unsigned i=0; i< voronoi_mesh.GetNumElements(); i++)
        {
            volume += voronoi_mesh.GetVolumeOfElement(i);
        }
        TS_ASSERT_DELTA(volume, 31.5, 1e-4); // Agrees with Paraview
    }

    void TestGetMeshForVtk()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Test GetMeshForVtk() method
        VertexMesh<2, 2>* p_mesh_for_vtk = mesh.GetMeshForVtk();

        // The mesh for VTK should have the same number of elements and nodes
        TS_ASSERT_EQUALS(p_mesh_for_vtk->GetNumElements(), 9u);
        TS_ASSERT_EQUALS(p_mesh_for_vtk->GetNumNodes(), 30u);
    }
};

#endif /*TESTVERTEXMESH_HPP_*/
