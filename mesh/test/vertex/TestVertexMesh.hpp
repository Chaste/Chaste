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

public:

    void TestNodeIterator() throw (Exception)
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
    }

    void TestVertexElementIterator() throw (Exception)
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

    void TestBasic1dVertexMesh() throw(Exception)
    {
        // Make four nodes to assign to three elements
        std::vector<Node<1>*> nodes;
        for (unsigned i=0; i<4; i++)
        {
            nodes.push_back(new Node<1>(i, false, 0.5*(double)i));
        }

        // Make three elements out of these nodes
        std::vector<std::vector<Node<1>*> > nodes_elements(3);
        std::vector<VertexElement<1,1>*> elements;
        for (unsigned i=0; i<3; i++)
        {
            nodes_elements[i].push_back(nodes[i]);
            nodes_elements[i].push_back(nodes[i+1]);

            elements.push_back(new VertexElement<1,1>(i, nodes_elements[i]));
        }

        // Make a vertex mesh
        VertexMesh<1,1> mesh(nodes, elements);

        // Test the mesh has the correct number of nodes and elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);

        // Test the elements have the correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNodeLocation(0)[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNodeLocation(1)[0], 0.5, 1e-6);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNodeLocation(0)[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNodeLocation(1)[0], 1.0, 1e-6);

        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh.GetElement(2)->GetNodeLocation(0)[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetElement(2)->GetNodeLocation(1)[0], 1.5, 1e-6);

        // Test GetCentroidOfElement method
        TS_ASSERT_DELTA(mesh.GetCentroidOfElement(0)[0], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh.GetCentroidOfElement(1)[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh.GetCentroidOfElement(2)[0], 1.25, 1e-6);
    }

    void TestBasic2dVertexMesh() throw(Exception)
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
        std::vector<std::vector<Node<2>*> > nodes_elements(2);
        for (unsigned i=0; i<5; i++)
        {
            nodes_elements[0].push_back(basic_nodes[i]);
        }
        nodes_elements[1].push_back(basic_nodes[2]);
        nodes_elements[1].push_back(basic_nodes[5]);
        nodes_elements[1].push_back(basic_nodes[6]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elements[0]));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elements[1]));

        // Make a vertex mesh
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 7u);

        TS_ASSERT_DELTA(basic_vertex_mesh.GetNode(2)->rGetLocation()[0], 1.5, 1e-3);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(),6u);

        // Check that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);

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

        // Coverage
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveBoundaryElementMapping(0), 0u);
    }

    void TestBasic3dVertexMesh()
    {
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

        // Coverage
        TS_ASSERT_EQUALS(p_mesh->SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveBoundaryElementMapping(0), 0u);

        delete p_mesh;
    }

    void TestGetCentroidOfElement() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements);

        // Test GetCentroidOfElement() method
        c_vector<double, 2> centroid = mesh.GetCentroidOfElement(0);

        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6);
    }

    void TestGetAreaGradientOfElementAtNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements);

        // Test GetAreaGradientOfElementAtNode() method at each node
        VertexElement<2,2>* p_element = mesh.GetElement(0);

        c_vector<double, 2> element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 0);
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);

        element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 1);
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);

        element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 2);
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);

        element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 3);
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);
    }

    void TestVertexElementAreaAndPerimeter()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);

        // Check nodes have correct indices
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), i);
        }

        // Test area and perimeter calculations
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(0), 4.0, 1e-6);
    }

    void TestVertexElementAreaAndPerimeterOnCircle()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);

        //  Check nodes have correct indices
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), i);
        }

        // Test area and perimeter calculations
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), M_PI, 1e-4);
        TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(0), 2.0*M_PI, 1e-4);
    }

    void Test3dMethodsWithPrism()
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
        VertexMesh<3,3> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumFaces(), 5u);

        // Face 0 has four vertices, is perpendicular to the y axis, and has area 1*3 = 3
        VertexElement<2,3>* p_face_0 = mesh.GetFace(0);
        TS_ASSERT_EQUALS(p_face_0->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_0 = mesh.GetUnitNormalToFace(p_face_0);
        TS_ASSERT_DELTA(unit_normal_0[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_0[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_0[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_0), 3.0, 1e-6);

        // Face 1 has three vertices, is perpendicular to the x axis, and has area 0.5*2*3 = 3
        VertexElement<2,3>* p_face_1 = mesh.GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 3u);
        c_vector<double, 3> unit_normal_1 = mesh.GetUnitNormalToFace(p_face_1);
        TS_ASSERT_DELTA(unit_normal_1[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_1), 3.0, 1e-6);

        // Face 2 has four vertices, is at an angle theta to the y axis where tan(theta) = 2/3,
        // and has area 1*sqrt(2^2 + 3^2) = sqrt(13)
        VertexElement<2,3>* p_face_2 = mesh.GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_2 = mesh.GetUnitNormalToFace(p_face_2);
        TS_ASSERT_DELTA(unit_normal_2[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[1], -sin(atan2(3,2)), 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[2], -cos(atan2(3,2)), 1e-6);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_2), sqrt(13), 1e-6);

        // Face 1 has three vertices, is perpendicular to the x axis, and has area 0.5*2*3 = 3
        VertexElement<2,3>* p_face_3 = mesh.GetFace(3);
        TS_ASSERT_EQUALS(p_face_3->GetNumNodes(), 3u);
        c_vector<double, 3> unit_normal_3 = mesh.GetUnitNormalToFace(p_face_3);
        TS_ASSERT_DELTA(unit_normal_3[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_3[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_3[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_3), 3.0, 1e-6);

        // Face 4 has four vertices, is perpendicular to the z axis, and has area 1*2 = 2
        VertexElement<2,3>* p_face_4 = mesh.GetFace(4);
        TS_ASSERT_EQUALS(p_face_4->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_4 = mesh.GetUnitNormalToFace(p_face_4);
        TS_ASSERT_DELTA(unit_normal_4[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_4[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_4[2], -1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_4), 2.0, 1e-6);

        // The volume of the prism should be 0.5 * 3 * 2 * 1 = 3
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 3.0, 1e-6);

        // The volume of the prism should be the sum of the face areas
        TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(0), 11 + sqrt(13), 1e-6);

        // By symmetry, the centroid of the prism should lie in the plane x=0.5
        c_vector<double, 3> centroid = mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(centroid(0), 0.5, 1e-5);
    }

    void TestGetPerimeterGradientAtNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements);

        // Test gradient of area evaluated at each node
        VertexElement<2,2>* p_element = mesh.GetElement(0);

        c_vector<double, 2> element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 0);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 1);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 2);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);

        element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 3);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);
    }

    void TestMeshGetWidthAndBoundingBoxMethod()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Test CalculateBoundingBox() method
        ChasteCuboid<2> bounds = mesh.CalculateBoundingBox();
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], 3.50,   1e-4);
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], 2.8867, 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], 0.0,    1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], 0.0,    1e-4);

        // Test GetWidth() method
        double width = mesh.GetWidth(0);
        double height = mesh.GetWidth(1);

        TS_ASSERT_DELTA(height, 2.8867, 1e-4);
        TS_ASSERT_DELTA(width, 3.5000, 1e-4);
    }

    void TestMeshConstructionFromMeshReader()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

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

        // Create mesh in which elements have attributes
        VertexMeshReader<2,2> mesh_reader2("mesh/test/data/TestVertexMesh/vertex_mesh_with_attributes");
        VertexMesh<2,2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 2u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh2.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh2.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh2.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        // Check second element has the right nodes
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNode(1), mesh2.GetNode(5));

        // Check element attributes
        TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetRegion(), 76u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetRegion(), 89u);
    }

    void TestMeshConstructionFromMeshReaderIndexedFromOne()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh_elements_indexed_from_1");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        // Check first element has the right nodes
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
    }


    void TestMeshConstructionFromMeshReaderWithFaces()
    {
        // Create mesh
        VertexMeshReader<3,3> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_3d_with_faces");
        VertexMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes, elements and faces
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
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_0), 1.125, 1e-4);

        VertexElement<2,3>* p_face_1 = mesh.GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_1), 1.9485, 1e-4);

        VertexElement<2,3>* p_face_2 = mesh.GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_2), 1.125, 1e-4);

        VertexElement<2,3>* p_face_3 = mesh.GetFace(3);
        TS_ASSERT_EQUALS(p_face_3->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(mesh.GetAreaOfFace(p_face_3), 1.125, 1e-4);

        // Check Voronoi element is correct
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumFaces(), 4u);
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 0.6495, 1e-4);
        TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(0), 5.3235, 1e-4);

        // Create mesh in which elements have attributes
        VertexMeshReader<3,3> mesh_reader2("mesh/test/data/TestVertexMeshWriter/vertex_mesh_3d_with_attributes");
        VertexMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);

        // Check we have the right number of nodes, elements and faces
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumFaces(), 4u);

        // Check element attributes
        TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetRegion(), 49u);
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

            // Tidy up
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

            // Tidy up
            delete p_mesh_loaded;
        }

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

    void TestNeighbouringNodeAndElementMethods() throw(Exception)
    {
        // Create mesh
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

    void TestCalculateMomentOfElement() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.0, 1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> small_mesh(nodes, elements);

        // Test CalculateMomentOfElement() method
        c_vector<double, 3> moments = small_mesh.CalculateMomentsOfElement(0);

        TS_ASSERT_DELTA(moments(0), 5.0/90.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(moments(1), 2.0/9.0, 1e-6);   // Iyy
        TS_ASSERT_DELTA(moments(2), -5.0/90.0, 1e-6); // Ixy

        // Hexagonal mesh from mesh generator
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 48u);

        // Test area and perimeter calculations for all elements
        for (VertexMesh<2,2>::VertexElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();

            TS_ASSERT_DELTA(mesh.GetVolumeOfElement(elem_index), 0.8660, 1e-4);
            TS_ASSERT_DELTA(mesh.GetSurfaceAreaOfElement(elem_index), 3.4641, 1e-4);
        }

        // Test centroid calculations for random elements
        c_vector<double, 2> centroid = mesh.GetCentroidOfElement(5);
        TS_ASSERT_DELTA(centroid(0), 2.0, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 1.4433, 1e-4);

        centroid = mesh.GetCentroidOfElement(7);
        TS_ASSERT_DELTA(centroid(0), 4.0, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 1.4433, 1e-4);

        // Test CalculateMomentOfElement() for all elements
        // all elements are regular hexagons with edge 1/sqrt(3)
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            moments = mesh.CalculateMomentsOfElement(i);

            TS_ASSERT_DELTA(moments(0), 5*sqrt(3)/16/9, 1e-6); // Ixx
            TS_ASSERT_DELTA(moments(1), 5*sqrt(3)/16/9, 1e-6); // Iyy
            TS_ASSERT_DELTA(moments(2), 0.0, 1e-6); // Ixy
        }
    }

    void TestGetShortAxisOfElement() throw(Exception)
    {
        // First test

        // Create nodes: this is a rectangle, centre (0,0), width 4, height 2, parallel to x axis
        std::vector<Node<2>*> nodes1;
        nodes1.push_back(new Node<2>(0, false,  2.0,  1.0));
        nodes1.push_back(new Node<2>(1, false, -2.0,  1.0));
        nodes1.push_back(new Node<2>(2, false, -2.0, -1.0));
        nodes1.push_back(new Node<2>(3, false,  2.0, -1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements1;
        elements1.push_back(new VertexElement<2,2>(0, nodes1));

        // Create mesh
        VertexMesh<2,2> mesh1(nodes1, elements1);

        // Test GetShortAxisOfElement() method
        c_vector<double, 2> short_axis = mesh1.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), 1.0, 1e-6);

        // Second test

        // Create nodes: this is a rectangle, centre (0,0), width 2, height 4, parallel to x axis
        std::vector<Node<2>*> nodes2;
        nodes2.push_back(new Node<2>(0, false,  1.0,  2.0));
        nodes2.push_back(new Node<2>(1, false, -1.0,  2.0));
        nodes2.push_back(new Node<2>(2, false, -1.0, -2.0));
        nodes2.push_back(new Node<2>(3, false,  1.0, -2.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements2;
        elements2.push_back(new VertexElement<2,2>(0, nodes2));

        // Create mesh
        VertexMesh<2,2> mesh2(nodes2, elements2);

        // Test GetShortAxisOfElement() method
        short_axis = mesh2.GetShortAxisOfElement(0);

        TS_ASSERT_DELTA(short_axis(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), 0.0, 1e-6);

        // Third test

        /*
         * Create nodes: this is a trapezoid, width 1, top length 3*sqrt(3), bottom length sqrt(3),
         * rotated by 30 degrees anticlockwise
         */
        std::vector<Node<2>*> nodes3;
        nodes3.push_back(new Node<2>(0, false,  1.0, 0.0));
        nodes3.push_back(new Node<2>(1, false,  2.0, sqrt(3.0)));
        nodes3.push_back(new Node<2>(2, false, -2.5, -sqrt(3.0)/2.0));
        nodes3.push_back(new Node<2>(3, false, -0.5, -sqrt(3.0)/2.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements3;
        elements3.push_back(new VertexElement<2,2>(0, nodes3));

        // Create mesh
        VertexMesh<2,2> mesh3(nodes3, elements3);

        // Test GetShortAxisOfElement() method
        short_axis = mesh3.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), -sqrt(3.0)*0.5, 1e-6);

        // Fourth test

        // Test on a regular polygon (generates a random vector)
        std::vector<Node<2>*> nodes4;
        unsigned num_nodes = 6;   // vertices
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes4.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        std::vector<VertexElement<2,2>*> elements4;
        elements4.push_back(new VertexElement<2,2>(0, nodes4));

        // Create mesh
        VertexMesh<2,2> mesh4(nodes4, elements4);

        // Test GetShortAxisOfElement() method
        short_axis = mesh4.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(short_axis(0)*short_axis(0)+short_axis(1)*short_axis(1), 1.0, 1e-6);

        // This is the same as seeding the random axis
        TS_ASSERT_DELTA(short_axis(0), 0.8401, 1e-4);
        TS_ASSERT_DELTA(short_axis(1), 0.5422, 1e-4);
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

        // Test the translate method
        // Pick a certain node and store spatial position
        Node<3>* p_node = mesh3d.GetNode(7);
        ChastePoint<3> original_coordinate = p_node->GetPoint();

        const double x_movement = 1.0;
        const double y_movement = 2.5;
        const double z_movement = 2.5;

        mesh3d.Translate(x_movement, y_movement, z_movement);

        ChastePoint<3>  new_coordinate = p_node->GetPoint();

        TS_ASSERT_DELTA(original_coordinate[0], new_coordinate[0] - x_movement, 1e-6);
        TS_ASSERT_DELTA(original_coordinate[1], new_coordinate[1] - y_movement, 1e-6);
        TS_ASSERT_DELTA(original_coordinate[2], new_coordinate[2] - z_movement, 1e-6);
    }

    void TestBoundaryNodes()
    {
        // Create a mesh with just boundary nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        VertexMesh<2,2> mesh1(nodes, elements);

        // Test boundary property of nodes
        for (unsigned i=0; i<mesh1.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh1.GetNode(i)->IsBoundaryNode(), false);
        }

        // Create a mesh with some interior nodes
        VertexMeshReader<2,2> mesh_reader1("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_2_by_2");
        VertexMesh<2,2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader1);

        // Test boundary property of nodes
        for (unsigned i=0; i<mesh2.GetNumNodes(); i++)
        {
            bool expected_boundary_node = (i==6 || i==9) ? false : true;
            TS_ASSERT_EQUALS(mesh2.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

        // Create mesh
        VertexMeshReader<2,2> mesh_reader2("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> mesh3;
        mesh3.ConstructFromMeshReader(mesh_reader2);

        // Test boundary property of nodes
        for (unsigned i=0; i<mesh3.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==8 || i==9 || i==12 || i==13 || i==16 || i==17 || i==20 || i==21)
            {
                expected_boundary_node = false;
            }

            TS_ASSERT_EQUALS(mesh3.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
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

    void TestTranslation2DMethod() throw (Exception)
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

    void TestGenerateVerticesFromElementCircumcentres() throw (Exception)
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

    void TestTessellationConstructor2d() throw (Exception)
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

    void TestTessellationConstructor2dWithGhostNodes() throw (Exception)
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

        // Test the number of nodes owned by each Voronoi element
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(4)->GetNumNodes(), 4u);

        // Test element areas
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(4), 0.5, 1e-6);
    }

    void TestGetEdgeLengthWithSimpleMesh() throw (Exception)
    {
        // Create a simple 2D tetrahedral mesh, the Delaunay triangulation
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0.0, 0));
        nodes.push_back(new Node<2>(0, true,  0.0, 1));
        nodes.push_back(new Node<2>(0, true, -1.0, 0));
        nodes.push_back(new Node<2>(0, true,  1.0, 0));
        nodes.push_back(new Node<2>(0, true,  0.5, -sqrt(3.0)/2.0));
        nodes.push_back(new Node<2>(0, true, -0.5, -sqrt(3.0)/2.0));

        MutableMesh<2,2> delaunay_mesh(nodes);
        TS_ASSERT_EQUALS(delaunay_mesh.CheckIsVoronoi(), true);

        // Create Voronoi tessellation
        VertexMesh<2,2> voronoi_mesh(delaunay_mesh);

        TS_ASSERT_EQUALS(voronoi_mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetNumNodes(), 5u);

        //Measure the length of the edge separating the centre element and each of its neighbours
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,1), 1.0, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,2), 0.5 + 1./(sqrt(3.0)*2.0), 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,3), 0.5 + 1./(sqrt(3.0)*2.0), 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,4), 1./sqrt(3.0), 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetEdgeLength(0,5), 1./sqrt(3.0), 1e-6);

        //All other neighbouring elements share an infinite edge
        TS_ASSERT_THROWS_THIS(voronoi_mesh.GetEdgeLength(1,2), "Elements 1 and  2 share only one node.");

        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(0), sqrt(3.0)/4.0+0.5, 1e-6);
        TS_ASSERT_DELTA(voronoi_mesh.GetSurfaceAreaOfElement(0), 2.0 + sqrt(3.0), 1e-6);
    }

    void TestTessellationConstructor3dWithGhostNode() throw (Exception)
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

        // Create Voronoi tessellation
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

        // Check Voronoi faces are correct
        VertexElement<2,3>* p_face_0 = voronoi_mesh.GetFace(0);
        TS_ASSERT_EQUALS(p_face_0->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_0->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(voronoi_mesh.GetAreaOfFace(p_face_0), 1.125, 1e-4);

        VertexElement<2,3>* p_face_1 = voronoi_mesh.GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_1->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_DELTA(voronoi_mesh.GetAreaOfFace(p_face_1), 1.9485, 1e-4);

        VertexElement<2,3>* p_face_2 = voronoi_mesh.GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(p_face_2->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(voronoi_mesh.GetAreaOfFace(p_face_2), 1.125, 1e-4);

        VertexElement<2,3>* p_face_3 = voronoi_mesh.GetFace(3);
        TS_ASSERT_EQUALS(p_face_3->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_face_3->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_DELTA(voronoi_mesh.GetAreaOfFace(p_face_3), 1.125, 1e-4);

        // Check Voronoi element is correct
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(voronoi_mesh.GetElement(0)->GetNumFaces(), 4u);
        TS_ASSERT_DELTA(voronoi_mesh.GetVolumeOfElement(0), 0.6495, 1e-4);
        TS_ASSERT_DELTA(voronoi_mesh.GetSurfaceAreaOfElement(0), 5.3235, 1e-4);
    }
};

#endif /*TESTVERTEXMESH_HPP_*/
