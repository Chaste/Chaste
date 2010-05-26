/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef TESTMUTABLEVERTEXMESH_HPP_
#define TESTMUTABLEVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "HoneycombMutableVertexMeshGenerator.hpp"
#include "VertexMeshWriter.hpp"
#include "MutableVertexMesh.hpp"
#include "ArchiveOpener.hpp"

class TestMutableVertexMesh : public CxxTest::TestSuite
{
public:

    void TestMutableVertexElementIterator() throw (Exception)
    {
        // Create mesh
        HoneycombMutableVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMutableMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 9u);

        unsigned counter = 0;
        for (MutableVertexMesh<2,2>::VertexElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), counter);

        // For coverage, test with an empty mesh
        MutableVertexMesh<2,2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        MutableVertexMesh<2,2>::VertexElementIterator iter = empty_mesh.GetElementIteratorBegin();

        // Check that the iterator is now at the end (we need to check this as a double-negative,
        // as we only have a NOT-equals operator defined on the iterator).
        bool iter_is_not_at_end = (iter != empty_mesh.GetElementIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);

        // Delete an element from mesh and test the iterator
        p_mesh->DeleteElementPriorToReMesh(0);

        counter = 0;
        for (MutableVertexMesh<2,2>::VertexElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter+1, element_index); // assumes the iterator will give elements 1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), counter);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), counter+1);
        TS_ASSERT_EQUALS(p_mesh->IsMeshChanging(), true);
    }

    void TestBasic2dMutableVertexMesh() throw(Exception)
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

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;

        // Make two triangular elements out of these nodes
        nodes_elem_0.push_back(basic_nodes[0]);
        nodes_elem_0.push_back(basic_nodes[1]);
        nodes_elem_0.push_back(basic_nodes[2]);
        nodes_elem_0.push_back(basic_nodes[3]);
        nodes_elem_0.push_back(basic_nodes[4]);

        nodes_elem_1.push_back(basic_nodes[2]);
        nodes_elem_1.push_back(basic_nodes[5]);
        nodes_elem_1.push_back(basic_nodes[6]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

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

        // Test Set and Get methods
        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.01, 1e-4); // Default value
        TS_ASSERT_DELTA(basic_vertex_mesh.GetEdgeDivisionThreshold(), DBL_MAX, 1e-4); // Default value
        TS_ASSERT_DELTA(basic_vertex_mesh.GetT2Threshold(), 0.001, 1e-4); // Default value

        basic_vertex_mesh.SetCellRearrangementThreshold(0.03);
        basic_vertex_mesh.SetEdgeDivisionThreshold(3.0);
        basic_vertex_mesh.SetT2Threshold(0.003);

        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.03, 1e-4);
        TS_ASSERT_DELTA(basic_vertex_mesh.GetEdgeDivisionThreshold(), 3.0, 1e-4);
        TS_ASSERT_DELTA(basic_vertex_mesh.GetT2Threshold(), 0.003, 1e-4);

        // Coverage
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveBoundaryElementMapping(0), 0u);
    }

    void TestMeshConstructionFromMeshReader()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Test Get methods
        TS_ASSERT_DELTA(mesh.GetCellRearrangementThreshold(), 0.01, 1e-4); // Default value
        TS_ASSERT_DELTA(mesh.GetEdgeDivisionThreshold(), DBL_MAX, 1e-4); // Default value
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
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh");
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

    void TestAddNodeAndReMesh() throw (Exception)
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh");
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

    void TestAddElement() throw (Exception)
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

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;

        // Make two triangular elements out of these nodes
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[4]);

        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[6]);

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

    void TestDivideVertexElementGivenNodes() throw(Exception)
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
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 0u);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 1u);

        // For coverage, divide an element when mDeletedElementIndices is not empty
        basic_vertex_mesh.DeleteElementPriorToReMesh(0);
        new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(1), 2, 3, true);

        TS_ASSERT_EQUALS(new_element_index, 0u);
    }

    void TestDivideVertexElementGivenNodesForCoverage() throw(Exception)
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
        nodes_elem.push_back(basic_nodes[0]);
        nodes_elem.push_back(basic_nodes[1]);
        nodes_elem.push_back(basic_nodes[2]);
        nodes_elem.push_back(basic_nodes[3]);

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

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 3u);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 3u);

        // For coverage, divide an element when mDeletedElementIndices is not empty
        basic_vertex_mesh.DeleteElementPriorToReMesh(0);
        new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(1), 2, 3, true);

        TS_ASSERT_EQUALS(new_element_index, 0u);
    }


    void TestDivideVertexElementAbove() throw(Exception)
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
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 1u);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 0u);
    }

    // This also tests that boundary nodes are updated on element division
    void TestDivideVertexElementGivenAxisOfDivision() throw(Exception)
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
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);

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

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 3u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 2u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        //Test boundary nodes updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(0)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(1)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(2)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(3)->IsBoundaryNode(), false);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->IsBoundaryNode(), false);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->IsBoundaryNode(), false);

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

    void TestDivideVertexElementWithBoundaryNodes() throw(Exception)
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

        // Make a square element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);

        // Make a square element out of nodes 1,4,5,2
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[2]);

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

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 7u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(4), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(3), 7u);

        // Test boundary nodes updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(0)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(1)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(2)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(3)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->IsBoundaryNode(), false);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(7)->IsBoundaryNode(), true);

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

    void TestDeleteElementWithBoundaryNodes() throw(Exception)
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
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[7]);
        nodes_elem_0.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[6]);
        nodes_elem_1.push_back(nodes[7]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[5]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Create mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Test mesh has correct numbers of elements and nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Test correct nodes are boundary nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(0)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(1)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(2)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(3)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(7)->IsBoundaryNode(), false);

        // Delete smaller element
        vertex_mesh.DeleteElementPriorToReMesh(1);

        // Test mesh has correct numbers of elements and nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Test ownership
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 7u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(4), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(5), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(6), 6u);

        // Test correct nodes are boundary nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(0)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(1)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(2)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(3)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->IsBoundaryNode(), true);
    }

    /**
     * Test that in the case where the given axis of division does not
     * cross two edges of the element, an exception is thrown.
     */
    void TestDivideVertexElementGivenAxisOfDivisionFailsForBadElement() throw(Exception)
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

    void TestDivideVertexElementAlongShortAxis() throw(Exception)
    {
        // Make five nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, false, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 2.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);

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

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);

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

    void TestDivideVertexElementWithNonRegularElement() throw(Exception)
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
        nodes_elem.push_back(nodes[0]);
        nodes_elem.push_back(nodes[1]);
        nodes_elem.push_back(nodes[2]);
        nodes_elem.push_back(nodes[3]);
        nodes_elem.push_back(nodes[4]);

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);

        // Divide element using two given nodes
        unsigned new_element_index = mesh.DivideElementAlongShortAxis(mesh.GetElement(0), true);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(4)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(0)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(3)->GetIndex(), 6u);

        // Test locations of new nodes
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[0], 2.3735, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[1], 1.3735, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 1.6520, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 2.3260, 1e-4);
    }

    void TestDivideVertexElementWhereNewNodesAreCloseToOldNodes1() throw(Exception)
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
        nodes_elem.push_back(nodes[0]);
        nodes_elem.push_back(nodes[1]);
        nodes_elem.push_back(nodes[2]);
        nodes_elem.push_back(nodes[3]);
        nodes_elem.push_back(nodes[4]);
        nodes_elem.push_back(nodes[5]);

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
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(0)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(3)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(4)->GetIndex(), 7u);

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(4)->GetIndex(), 5u);

        // Test locations of new nodes
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], -0.009+1.5*mesh.GetCellRearrangementThreshold(), 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[1], 1.0, 1e-4);
    }

    void TestDivideVertexElementWhereNewNodesAreCloseToOldNodes2() throw(Exception)
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
        nodes_elem.push_back(nodes[0]);
        nodes_elem.push_back(nodes[1]);
        nodes_elem.push_back(nodes[2]);
        nodes_elem.push_back(nodes[3]);
        nodes_elem.push_back(nodes[4]);
        nodes_elem.push_back(nodes[5]);

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
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(0)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(3)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(4)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(5)->GetIndex(), 7u);

        // Test locations of new nodes
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 0.009-1.5*mesh.GetCellRearrangementThreshold(), 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(7)->rGetLocation()[1], 1.0, 1e-4);
    }

};

#endif /*TESTMUTABLEVERTEXMESH_HPP_*/
