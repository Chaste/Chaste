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

#ifndef TESTMUTABLEVERTEXMESHREMESH_HPP_
#define TESTMUTABLEVERTEXMESHREMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "VertexMeshWriter.hpp"
#include "MutableVertexMesh.hpp"

class TestMutableVertexMeshReMesh : public CxxTest::TestSuite
{
public:

    /*
     * This tests both PerformNodeMerge and IdentifySwapType.
     */
    void TestPerformNodeMerge() throw(Exception)
    {
        // Create seven nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, -0.1, -0.1));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, -1.0, 1.0));
        nodes.push_back(new Node<2>(5, false, -1.0, 0.0));
        nodes.push_back(new Node<2>(6, false, 0.1, -0.1));
        nodes.push_back(new Node<2>(7, false, 0.0, 0.1));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;

        // Create three elements containing these nodes
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[7]);

        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[6]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[7]);
        nodes_elem_1.push_back(nodes[3]);
        nodes_elem_1.push_back(nodes[4]);

        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Merge nodes 0 and 6 (node 0 is in elements 1 and 2, node 6 is in element 1)
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(0), vertex_mesh.GetNode(6), map);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Test nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(0)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(0)->rGetLocation()[1], -0.1, 1e-8);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(4)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 5u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.95, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 2.9+sqrt(1.01), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.65,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.9+sqrt(2.21)+2.0*sqrt(1.01), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.5,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.0+sqrt(2.21)+sqrt(1.01), 1e-6);
    }

    /*
     * This test provides coverage of the case in which, when the elements
     * previously containing the high-index node are updated to contain the
     * low-index node, at least one of these elements did not already contain
     * the low-index node.
     */
    void TestPerformNodeMergeWhenLowIndexNodeMustBeAddedToElement() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 1.01, 1.0));
        nodes.push_back(new Node<2>(5, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(6, false, 0.0, 2.0));

        // Create two elements containing nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[5]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Merge nodes 4 and 5
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), map);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 1.005, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 1.0, 1e-8);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 4u);
    }

    // This tests both PerformNodeMerge and IdentifySwapType
    void TestPerformNodeMergeOnEdge() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.4, 0.0));
        nodes.push_back(new Node<2>(5, false, 0.6, 0.0));
        nodes.push_back(new Node<2>(6, false, 0.4, 0.4));
        nodes.push_back(new Node<2>(7, false, 0.6, 0.6));

        // Create two elements containing nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[7]);
        nodes_elem_0.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[6]);
        nodes_elem_1.push_back(nodes[7]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Merge nodes 6 and 7
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(6), vertex_mesh.GetNode(7), map);

        // Merge nodes 4 and 5
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), map);

        // Test that the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test merged node is in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.0, 1e-3);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(4)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 3u);

        // Test Areas and Perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 2+sqrt(2), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.5,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 2.0+sqrt(2), 1e-6);
    }

    void TestAnotherPerformNodeMerge() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 3.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 3.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 3.0, 2.0));
        nodes.push_back(new Node<2>(5, false, 2.0, 2.0));
        nodes.push_back(new Node<2>(6, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(7, false, 0.99, 1.0));
        nodes.push_back(new Node<2>(8, false, 0.0, 1.0));

        // Create three elements containing nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[6]);
        nodes_elem_0.push_back(nodes[7]);
        nodes_elem_0.push_back(nodes[8]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        nodes_elem_1.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);
        nodes_elem_2.push_back(nodes[6]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));

        // Create a mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 9u);

        // Merge nodes 6 and 7
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(6), vertex_mesh.GetNode(7), map);

        // Test that the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Test merged node is in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 0.995, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 1.0, 1e-3);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 7u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 6u);
    }

    // This tests both PerformT1Swap and IdentifySwapType
    void TestPerformT1Swap() throw(Exception)
    {
        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 0.4));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.6));

        // Make two triangular and two rhomboid elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[5]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[1]);

        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[0]);

        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[4]);
        nodes_elem_3.push_back(nodes[5]);
        nodes_elem_3.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(3), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(3), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(3), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(3), 1.2+0.2*sqrt(41.0), 1e-6);

        // Perform a T1 swap on nodes 4 and 5
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), map);

        // Test moved nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 1u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 0u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 3u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.2,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(3), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(3), 1.0+0.2*sqrt(41.0), 1e-6);
    }

    // This tests both PerformT1Swap and IdentifySwapType
    void TestPerformT1SwapOnBoundary() throw(Exception)
    {
        /* Make 6 nodes to assign to 3 elements all boundary nodes
         *
         * Note: this tests ensures coverage
         *  _____
         * |\   /
         * | \ /
         * |  |
         * | / \
         * |/___\
         *
         */

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.5, 0.4));
        nodes.push_back(new Node<2>(5, true, 0.5, 0.6));

        // Make two triangular and one rhomboid elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[5]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[0]);

        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);
        nodes_elem_2.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);

        // Perform a T1 swap on nodes 5 and 4. Note: this way round to ensure coverage of boundary node tracking.
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4), map);

        // Test moved nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 0u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 3u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if ( i==5 )
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    // This tests both PerformT1Swap and IdentifySwapType
    void TestPerformT1SwapOnBoundary2() throw(Exception)
    {
        /* Make 6 nodes to assign to 3 elements with 3 boundary nodes
         *
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         *
         */

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 0.4));
        nodes.push_back(new Node<2>(5, true, 0.5, 0.6));

        // Make one triangular and two rhomboid elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[4]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[0]);

        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);
        nodes_elem_2.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);

        // Perform a T1 swap on nodes 5 and 4. Note: this way round to ensure coverage of boundary node tracking.
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4), map);

        // Test moved nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 4u);


        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 0u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 3u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if ( i==0 || i==1 )
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    // This tests both PerformT1Swap and IdentifySwapType
    void TestPerformT1SwapToSeparate() throw(Exception)
    {
        /* This tests the following setup
         *
         * |\   /|     |\      /|
         * | \ / |     | \    / |
         * |  |  |  => | /    \ |
         * | / \ |     |/      \|
         * |/   \|
         *
         * Make 6 nodes to assign to 2 elements all boundary nodes
         */

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.5, 0.4));
        nodes.push_back(new Node<2>(5, true, 0.5, 0.6));

        // Make two rhomboid elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[3]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

        // Perform a T1 swap on nodes 5 and 4.
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4), map);

        // Test moved nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 3u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 2u);

        // Test areas and perimeters of elements
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }


    /*
     * This tests that T1Swaps rearange to form a Triangular element for a T2 Swap
     */
    void TestPrepareForT2Swap() throw(Exception)
    {
        // Make 8 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, -1.0, -1.0));
        nodes.push_back(new Node<2>(1, false,  1.0, -1.0));
        nodes.push_back(new Node<2>(2, false,  1.0,  1.0));
        nodes.push_back(new Node<2>(3, false, -1.0,  1.0));
        nodes.push_back(new Node<2>(4, false, -0.1, -0.1));
        nodes.push_back(new Node<2>(5, false,  0.1, -0.1));
        nodes.push_back(new Node<2>(6, false,  0.1,  0.1));
        nodes.push_back(new Node<2>(7, false, -0.1,  0.1));

        /*
         *  Make Four trapezium elements with a central square element out of these nodes
         *  _______
         * |\  2  /|
         * | \___/ |
         * | |   | |
         * |3| 4 |1|
         * | |___| |
         * | / 0 \ |
         * |/_____\|
         *
         */

        // Trapezium element
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[4]);

        // Trapezium element
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[6]);
        nodes_elem_1.push_back(nodes[5]);

        // Trapezium element
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[2]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[7]);
        nodes_elem_2.push_back(nodes[6]);

        // Trapezium element
        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[4]);
        nodes_elem_3.push_back(nodes[7]);
        nodes_elem_3.push_back(nodes[3]);

        // Central square element
        std::vector<Node<2>*> nodes_elem_4;
        nodes_elem_4.push_back(nodes[4]);
        nodes_elem_4.push_back(nodes[5]);
        nodes_elem_4.push_back(nodes[6]);
        nodes_elem_4.push_back(nodes[7]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(4, nodes_elem_4));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.25);// T1 threshold distance is 0.25 so inner edges are too short
        vertex_mesh.SetT2Threshold(0.001); //T2 threshold small so doesnt occur

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        vertex_mesh.ReMesh();

        /*
         *  T1 swap occurs on nodes 4 and 5, mesh now looks like
         *
         *  ______
         * |\ 2  /|
         * | \__/ |
         * |  \/  |
         * | 3 | 1|
         * |  /\  |
         * | / 0\ |
         * |/____\|
         *
         */

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Test moved nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], -0.2875, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.0875, 1e-3);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(4)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(3)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(4)->GetIndex(), 3u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(0)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(1)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(2)->GetIndex(), 7u);

        vertex_mesh.SetT2Threshold(0.1); //T2 threshold larger so swap does occur
        vertex_mesh.ReMesh();

        /*
         *  T2 swap occurs on Element 4, mesh now looks like.
         *  ______
         * |\ 2  /|
         * | \  / |
         * |  \/  |
         * |3  | 1|
         * |  /\  |
         * | / 0\ |
         * |/____\|
         *
         */

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(),4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(),4u); // Elements are deleted not just marked for deletion.
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(),6u);

        // Test nodes are merged in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.2875/3.0, 1e-3);

        // Test elements are OK
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(3)->GetIndex(), 3u);
    }

    void TestPerformT2Swap() throw(Exception)
    {
        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.4, 0.25));
        nodes.push_back(new Node<2>(4, false, 0.6, 0.25));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.3));

        // Make one triangular and three trapezium elements out of these nodes

        // Triangle element
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[5]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[4]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[2]);
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[5]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[1]);
        nodes_elem_3.push_back(nodes[4]);
        nodes_elem_3.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1);// Threshold distance set to ease calculations.
        vertex_mesh.SetT2Threshold(0.01);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Perform a T2 swap on the middle triangle element
        VertexElement<2,2>* p_element_0 = vertex_mesh.GetElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 6u);

        for (unsigned j=1; j<4; j++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumNodes(), 3u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(0)->GetIndex(), j%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(1)->GetIndex(), (j+1)%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(2)->GetIndex(), 3u);
        }
    }

    void TestT2SwapsDontOccurWithTriangularNeighbours() throw(Exception)
    {
        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.4, 0.25));
        nodes.push_back(new Node<2>(4, false, 0.6, 0.25));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.3));

        // Make two triangles and two trapezium elements out of these nodes

        // Triangle element
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[5]);

        // Triangle element
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[4]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[2]);
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[5]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[1]);
        nodes_elem_3.push_back(nodes[4]);
        nodes_elem_3.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements, 0.1);
        vertex_mesh.SetCellRearrangementThreshold(0.1);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);

        // Attempt to perform a T2 swap on the middle triangle element
        VertexElement<2,2>* p_element_0 = vertex_mesh.GetElement(0);
        TS_ASSERT_THROWS_THIS( vertex_mesh.PerformT2Swap(*p_element_0),
                "One of the neighbours of a small triangular element is also a triangle - "
                "dealing with this has not been implemented yet" );
    }

    /**
     * This tests the ReMesh method for preforming T2Swaps (element removal).
     */
    void TestRemeshForT2Swap() throw(Exception)
    {
        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.1, 0.05));
        nodes.push_back(new Node<2>(4, false, 0.9, 0.05));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.475));

        // Make one triangular and three trapezium elements out of these nodes

        // Triangle element
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[5]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[4]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[2]);
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[5]);

        // Trapezium
        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[1]);
        nodes_elem_3.push_back(nodes[4]);
        nodes_elem_3.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        vertex_mesh.SetT2Threshold(0.01);
        vertex_mesh.SetCellRearrangementThreshold(0.00001); //So T1Swaps dont happen

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        vertex_mesh.ReMesh(); // Elements too big so nothing happens

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        c_vector<double, 2>& new_location_0 = vertex_elements[0]->GetNode(0)->rGetModifiableLocation();
        new_location_0(0) = 0.499;
        new_location_0(1) = 0.249;

        c_vector<double, 2>& new_location_1 = vertex_elements[0]->GetNode(1)->rGetModifiableLocation();
        new_location_1(0) = 0.501;
        new_location_1(1) = 0.249;

        c_vector<double, 2>& new_location_2 = vertex_elements[0]->GetNode(2)->rGetModifiableLocation();
        new_location_2(0) = 0.5;
        new_location_2(1) = 0.251;

        // T2 swaps should now happen
        vertex_mesh.ReMesh();
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        // Test merged node is in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[0], 0.4999, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[1], 0.2496, 1e-3);

        // Test elements have correct nodes
        // note nodes are renumbered as element 0 is deleted

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 3u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 3u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 3u);
    }


    /**
     * This tests the ReMesh method for preforming T1Swaps, both internaly and on the boundary.
     * In this test we read in a vertex mesh that contains several pairs of nodes that
     * are close enough for T1Swaps to be performed. The mesh consists of 6 elements and all
     * T1Swaps are performed on all horizontal edges.
     *
     *      /\    /\
     *     /  \__/  \
     *    /   /  \   \
     *    \__/\__/\__/
     *    /  \/  \/  \
     *    \   \__/   /
     *     \  /  \  /
     *      \/    \/
     *
     * Note: this also tests that boundary nodes are updated accordingly
     */

    void TestReMeshForT1Swaps() throw(Exception)
    {
        // This also tests IdentifySwapType

        // LoadMesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_remesh_mesh_all");
        MutableVertexMesh<2,2> vertex_mesh;

        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 22u);

        // assign boundary nodes \todo #1076 - once reading/writing of boundary elements is done
        // properly for vertex meshes this can be added to the .node file
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            if (i==4 || i==5 || i==6 || i==7 || i==8 || i==9 || i==10 || i==11 || i==12 || i==13 || i==15 || i==16 || i==17 || i==18)
            {
                vertex_mesh.GetNode(i)->SetAsBoundaryNode(true);
            }
        }

        // Calls ReMesh to identify all T1 swaps and perform them.
        vertex_mesh.ReMesh();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 22u);

        std::string dirname = "vertex_remeshing_mesh";
        std::string mesh_filename = "vertex_mesh_all";

        // Save the mesh data using mesh writers
        VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh(vertex_mesh);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 20u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(4)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(5)->GetIndex(), 20u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 21u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(3)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(4)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(5)->GetIndex(), 21u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(0)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(1)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(3)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(4)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(5)->GetIndex(), 1u);


        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(1)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(2)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(3)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(4)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(5)->GetIndex(), 16u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(0)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(1)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(2)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(3)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(4)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(5)->GetIndex(), 19u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(1)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(2)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(3)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(4)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(5)->GetIndex(), 13u);

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = false;
            if (i==4 || i==5 || i==6 || i==7 || i==8 || i==9 || i==10 || i==11 || i==12 || i==14 || i==15 || i==17 || i==18 || i==19)
            {
                expected_boundary_node = true;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

    }

    /**
     * This tests the ReMesh method for preforming node merges, both internaly and on the boundary.
     *
     * In this test we read in a vertex mesh that contains several pairs of nodes that
     * are close enough to be merged.
     */
    void TestRemeshForMerge() throw(Exception)
    {
        // This also tests IdentifySwapType

        // Load in mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_merge_mesh_all");
        MutableVertexMesh<2,2> vertex_mesh;
        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 14u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 10u);

        // Identify and perform any cell rearrangments
        vertex_mesh.ReMesh();

        /* We should have performed five node merges and
         * 1 and 12 merge to 1
         * 5 and 11 merge to 5
         * 9 and 10 merge to 9 becomes 7 on renumbering
         * 4 and 8  merge to 4
         * 6 and 7  merge to 6
         * 13 becomes 8 on renumbering
         *
         * Nodes 9, 10, 11, 12 and 13 should have been removed
         */

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 9u); // this should be 8

        std::string dirname = "vertex_remeshing_mesh";
        std::string mesh_filename = "vertex_merge_mesh_all";

        // Save the mesh data using mesh writers
        VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh(vertex_mesh);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(4)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(4)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(5)->GetIndex(), 5u);
    }

    void TestReMeshExceptions() throw(Exception)
    {
        // This also tests IdentifySwapType

        // Create some nodes
        Node<2>* p_node0 = new Node<2>(0, false, 0.0, 0.0);
        Node<2>* p_node1 = new Node<2>(1, false, 1.0, 0.0);
        Node<2>* p_node2 = new Node<2>(2, false, 1.0, 1.0);
        Node<2>* p_node3 = new Node<2>(3, false, 0.0, 1.0);
        Node<2>* p_node4 = new Node<2>(4, false, 0.5, 0.5);
        Node<2>* p_node5 = new Node<2>(5, false, 0.49, 0.49);
        Node<2>* p_node6 = new Node<2>(6, false, 0.75, 0.75); // so that all elements have at least 4 nodes

        std::vector<Node<2>*> nodes_in_element0;
        nodes_in_element0.push_back(p_node0);
        nodes_in_element0.push_back(p_node1);
        nodes_in_element0.push_back(p_node4);
        nodes_in_element0.push_back(p_node5);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(p_node1);
        nodes_in_element1.push_back(p_node2);
        nodes_in_element1.push_back(p_node6);
        nodes_in_element1.push_back(p_node4);

        std::vector<Node<2>*> nodes_in_element2;
        nodes_in_element2.push_back(p_node2);
        nodes_in_element2.push_back(p_node3);
        nodes_in_element2.push_back(p_node4);
        nodes_in_element2.push_back(p_node6);

        std::vector<Node<2>*> nodes_in_element3;
        nodes_in_element3.push_back(p_node0);
        nodes_in_element3.push_back(p_node5);
        nodes_in_element3.push_back(p_node4);
        nodes_in_element3.push_back(p_node3);

        std::vector<Node<2>*> nodes;
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);
        nodes.push_back(p_node2);
        nodes.push_back(p_node3);
        nodes.push_back(p_node4);
        nodes.push_back(p_node5);
        nodes.push_back(p_node6);

        /* Create 4 joined triangular elements with an extra nodes at 'o'.
         *  ______
         * |\    /|
         * | \  o |
         * |  \/  |
         * |  o\  |
         * | /  \ |
         * |/____\|
         *
         */
        VertexElement<2,2>* p_element0 = new VertexElement<2,2>(0, nodes_in_element0);
        VertexElement<2,2>* p_element1 = new VertexElement<2,2>(1, nodes_in_element1);
        VertexElement<2,2>* p_element2 = new VertexElement<2,2>(2, nodes_in_element2);
        VertexElement<2,2>* p_element3 = new VertexElement<2,2>(3, nodes_in_element3);

        std::vector<VertexElement<2,2>* > elements;
        elements.push_back(p_element0);
        elements.push_back(p_element1);
        elements.push_back(p_element2);
        elements.push_back(p_element3);

        // Create mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);

        // Call remesh
        TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "A node is contained in more than three elements");
    }

    void TestReMeshDivideEdgeIfTooBig() throw(Exception)
    {
        // Create some nodes
        Node<2>* p_node0 = new Node<2>(0, false, 0.0, 0.0);
        Node<2>* p_node1 = new Node<2>(1, false, 0.5, -1.0);
        Node<2>* p_node2 = new Node<2>(2, false, 1.0, 0.0);
        Node<2>* p_node3 = new Node<2>(3, false, 0.5, 1.0);

        std::vector<Node<2>*> nodes_in_element0;
        nodes_in_element0.push_back(p_node0);
        nodes_in_element0.push_back(p_node1);
        nodes_in_element0.push_back(p_node3);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(p_node1);
        nodes_in_element1.push_back(p_node2);
        nodes_in_element1.push_back(p_node3);

        std::vector<Node<2>*> nodes;
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);
        nodes.push_back(p_node2);
        nodes.push_back(p_node3);

        // Create 2 joined triangular elements
        VertexElement<2,2>* p_element0 = new VertexElement<2,2>(0, nodes_in_element0);
        VertexElement<2,2>* p_element1 = new VertexElement<2,2>(1, nodes_in_element1);
        std::vector<VertexElement<2,2>* > elements;
        elements.push_back(p_element0);
        elements.push_back(p_element1);

        // Create mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);
        mesh.SetEdgeDivisionThreshold(1.5); // This needs to be set to allow edge division.

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Call remesh
        mesh.ReMesh();

        // Check that the edge between nodes 1 and 2 has divided
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[0], 0.5, 1e-8);
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[1], 0.0, 1e-8);
    }

    //\TODO include boundary nodes in the tests
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
        Node<2>* p_node0 = new Node<2>(0, false, 1.0, 1.0);
        Node<2>* p_node1 = new Node<2>(1, false, 2.0, 1.0);
        Node<2>* p_node2 = new Node<2>(2, false, 2.0, 2.0);
        Node<2>* p_node3 = new Node<2>(3, false, 1.0, 2.0);
        Node<2>* p_node4 = new Node<2>(4, false, 0.5, 1.5);
        Node<2>* p_node5 = new Node<2>(5, false, 2.5, 1.5);

        std::vector<Node<2>*> nodes_in_element0;
        nodes_in_element0.push_back(p_node0);
        nodes_in_element0.push_back(p_node3);
        nodes_in_element0.push_back(p_node4);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(p_node1);
        nodes_in_element1.push_back(p_node5);
        nodes_in_element1.push_back(p_node2);

        std::vector<Node<2>*> nodes_in_element2;
        nodes_in_element2.push_back(p_node0);
        nodes_in_element2.push_back(p_node1);
        nodes_in_element2.push_back(p_node2);
        nodes_in_element2.push_back(p_node3);

        std::vector<Node<2>*> nodes;
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);
        nodes.push_back(p_node2);
        nodes.push_back(p_node3);
        nodes.push_back(p_node4);
        nodes.push_back(p_node5);

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

        TS_ASSERT_DELTA(mesh.GetAreaOfElement(2), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(2), 4.0, 1e-6);

        // Test other nodes are updated

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(), 2u);

        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(1)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(2)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(3)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(4)->GetIndex(), 3u);

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
        mesh.DivideEdge(mesh.GetNode(3), mesh.GetNode(0));

        // Divide the edge joining nodes 2 and 1
        mesh.DivideEdge(mesh.GetNode(2), mesh.GetNode(1));

        // Test edges are divided
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
        TS_ASSERT_DELTA(mesh.GetAreaOfElement(2), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(2), 4.0, 1e-6);

        // Test other nodes are updated
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(3)->GetIndex(), 8u);

        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(1)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(2)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(3)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(4)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(5)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNode(6)->GetIndex(), 7u);

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
    }

    void TestElementIncludesPointAndGetLocalIndexForElementEdgeClosestToPoint()
    {
        // Make four nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Make element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Make mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);

        // Make some test points and test ElementIncludesPoint()

        // A point far outside the element
        c_vector<double, 2> test_point1;
        test_point1[0] = -1.0;
        test_point1[1] = -1.0;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point1, 0), false);

        // A point far inside the element
        c_vector<double, 2> test_point2;
        test_point2[0] = 0.5;
        test_point2[1] = 0.5;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point2, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point2, 0), 0u);

        // A point on a non-horizontal edge
        c_vector<double, 2> test_point3;
        test_point3[0] = 0.0;
        test_point3[1] = 0.5;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point3, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point3, 0), 3u);

        // A point on a horizontal edge
        c_vector<double, 2> test_point4;
        test_point4[0] = 0.5;
        test_point4[1] = 0.0;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point4, 0), false);

        // A point just inside the element
        c_vector<double, 2> test_point5;
        test_point5[0] = 0.999;
        test_point5[1] = 0.998;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point5, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point5, 0), 1u);

        // A point just outside the element
        c_vector<double, 2> test_point6;
        test_point6[0] = 1.001;
        test_point6[1] = 0.5;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point6, 0), false);

        // A point coinciding with a vertex
        c_vector<double, 2> test_point7;
        test_point7[0] = 1.0;
        test_point7[1] = 1.0;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point7, 0), false);
    }

    void TestT3Swap()
    {
        /*
         * Make a small mesh consisting of five elements:
         * a square and three triangles sat on top of a rectangle.
         *         _____
         *    |\  |     |  /|
         *    | \ |     | /_|
         *    | / |     | \ |
         *    |/__|_____|__\|
         *    |             |
         *    |_____________|
         */

        // Make all nodes boundary nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(5, true, 2.0, 1.0));
        nodes.push_back(new Node<2>(6, true, 1.1, 0.5));
        nodes.push_back(new Node<2>(7, true, -1.0, 0.0));
        nodes.push_back(new Node<2>(8, true, -0.1, 0.5));
        nodes.push_back(new Node<2>(9, true, -1.0, 1.0));
        nodes.push_back(new Node<2>(10, true, -1.0, -1.0));
        nodes.push_back(new Node<2>(11, true, 2.0, -1.0));
        nodes.push_back(new Node<2>(12, true, 2.0, 0.5));


        std::vector<Node<2>*> nodes_in_element0;
        nodes_in_element0.push_back(nodes[0]);
        nodes_in_element0.push_back(nodes[1]);
        nodes_in_element0.push_back(nodes[2]);
        nodes_in_element0.push_back(nodes[3]);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(nodes[4]);
        nodes_in_element1.push_back(nodes[12]);
        nodes_in_element1.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_in_element2;
        nodes_in_element2.push_back(nodes[12]);
        nodes_in_element2.push_back(nodes[5]);
        nodes_in_element2.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_in_element3;
        nodes_in_element3.push_back(nodes[7]);
        nodes_in_element3.push_back(nodes[8]);
        nodes_in_element3.push_back(nodes[9]);

        std::vector<Node<2>*> nodes_in_element4;
        nodes_in_element4.push_back(nodes[10]);
        nodes_in_element4.push_back(nodes[11]);
        nodes_in_element4.push_back(nodes[4]);
        nodes_in_element4.push_back(nodes[1]);
        nodes_in_element4.push_back(nodes[0]);
        nodes_in_element4.push_back(nodes[7]);

        // Make elements
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
        elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
        elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));
        elements.push_back(new VertexElement<2,2>(3, nodes_in_element3));
        elements.push_back(new VertexElement<2,2>(4, nodes_in_element4));

        // Make mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);
        mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);// Threshold distance set to ease calculations.

        // Node 6 is close to, but not overlapping, an edge of element 0
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(6)->rGetLocation(), 0), false);

        // Move node 6 to the left so that it overlaps element 1
        ChastePoint<2> point = mesh.GetNode(6)->GetPoint();
        point.SetCoordinate(0u, 0.9);
        mesh.SetNode(6, point);

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(6)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(mesh.GetNode(6)->rGetLocation(), 0), 1u);

        // Node 8 is close to, but not overlapping, an edge of element 0
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(8)->rGetLocation(), 0), false);

        // Move node 8 to the left so that it overlaps element 1
        point.SetCoordinate(0u, 0.1);
        mesh.SetNode(8, point);

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(8)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(mesh.GetNode(8)->rGetLocation(), 0), 3u);


        // Call method to update mesh in this situation
        mesh.ReMesh();

        // Check that node 6 has been moved onto the edge a new node has been created and both added to elements 0 amd 1
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 16u);

        // Test locations of moved and new nodes due to node 6
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.5, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(13)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(13)->rGetLocation()[1], 0.4, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(14)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(14)->rGetLocation()[1], 0.6, 1e-4);

         // Test locations of moved and new nodes due to node 8
        TS_ASSERT_DELTA(mesh.GetNode(8)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8)->rGetLocation()[1], 0.45, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(15)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(15)->rGetLocation()[1], 0.55, 1e-4);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 13u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(3), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(4), 14u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(5), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(6), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(7), 15u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(8), 8u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 12u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(3), 13u);

        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(0), 12u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(2), 14u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNodeGlobalIndex(0), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNodeGlobalIndex(1), 8u);
        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNodeGlobalIndex(2), 15u);
        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNodeGlobalIndex(3), 9u);


        //Other elements remain the same so get a void
        TS_ASSERT_EQUALS(mesh.GetElement(4)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(4)->GetNodeGlobalIndex(0), 10u);
        TS_ASSERT_EQUALS(mesh.GetElement(4)->GetNodeGlobalIndex(1), 11u);
        TS_ASSERT_EQUALS(mesh.GetElement(4)->GetNodeGlobalIndex(2), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(4)->GetNodeGlobalIndex(3), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(4)->GetNodeGlobalIndex(4), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(4)->GetNodeGlobalIndex(5), 7u);

        // Test boundary property of nodes. All are boundary nodes except node 6.
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==6)
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestPerformT3SwapExceptions() throw(Exception)
    {
        /* Create 3 joined triangular elements intesecting at a node inside a square element
         *  ______
         * |      |   /|
         * |      |  /_|
         * |      | // |
         * |      | \\_|
         * |      |  \ |
         * |______|   \|
         *
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.9, 0.5));
        nodes.push_back(new Node<2>(5, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, true, 2.0, 0.3));
        nodes.push_back(new Node<2>(7, true, 2.0, 0.7));
        nodes.push_back(new Node<2>(8, true, 2.0, 1.0));

        std::vector<Node<2>*> nodes_in_element0;
        nodes_in_element0.push_back(nodes[0]);
        nodes_in_element0.push_back(nodes[1]);
        nodes_in_element0.push_back(nodes[2]);
        nodes_in_element0.push_back(nodes[3]);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(nodes[5]);
        nodes_in_element1.push_back(nodes[6]);
        nodes_in_element1.push_back(nodes[4]);

        std::vector<Node<2>*> nodes_in_element2;
        nodes_in_element2.push_back(nodes[6]);
        nodes_in_element2.push_back(nodes[7]);
        nodes_in_element2.push_back(nodes[4]);

        std::vector<Node<2>*> nodes_in_element3;
        nodes_in_element3.push_back(nodes[7]);
        nodes_in_element3.push_back(nodes[8]);
        nodes_in_element3.push_back(nodes[4]);

        // Make elements
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
        elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
        elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));
        elements.push_back(new VertexElement<2,2>(3, nodes_in_element3));

        // Make mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);

        // Call remesh which in turn calls PerformT3Swap
        TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "Trying to merge a node, contained in more than 2 elements, into another element, this is not possible with the vertex mesh.");


        /* Create a rectangualar and a triangular node to test interscting on an
         * edge that is too small
         *
         *
         *  ______  /|
         * |      |/ | <---
         * |______|\ |
         *          \|
         *
         *
         */
        nodes.clear();
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 0.1));
        nodes.push_back(new Node<2>(3, true, 0.0, 0.1));
        nodes.push_back(new Node<2>(4, true, 0.99, 0.05));
        nodes.push_back(new Node<2>(5, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, true, 2.0, 0.1));

        nodes_in_element0.clear();
        nodes_in_element0.push_back(nodes[0]);
        nodes_in_element0.push_back(nodes[1]);
        nodes_in_element0.push_back(nodes[2]);
        nodes_in_element0.push_back(nodes[3]);

        nodes_in_element1.clear();
        nodes_in_element1.push_back(nodes[5]);
        nodes_in_element1.push_back(nodes[6]);
        nodes_in_element1.push_back(nodes[4]);

        // Make elements
        elements.clear();
        elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
        elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));

        // Make mesh
        MutableVertexMesh<2,2> vertex_mesh_2(nodes, elements);
        vertex_mesh_2.SetCellRearrangementThreshold(1.0);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh_2.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh_2.GetNumElements(), 2u);

        // Call remesh which in turn calls PerformT3Swap
        TS_ASSERT_THROWS_THIS(vertex_mesh_2.ReMesh(), "Trying to merge a node onto an edge which is too small.");



    }

    void TestT3SwapForNeighboringElements()
    {
        /*
         * Make a small mesh consisting of 4 elements:
         * a square and a three triangles.
         *         _____
         *        |     |
         *     /\ |     | /|\
         *    /__\|_____|/_|_\
         *
         */

        // Make all nodes boundary nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 1.5, 0.0));
        nodes.push_back(new Node<2>(5, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, true, 1.5, 0.5));
        nodes.push_back(new Node<2>(7, true, -0.1, 0.0));
        nodes.push_back(new Node<2>(8, true, -0.5, 0.5));

        std::vector<Node<2>*> nodes_in_element0;
        nodes_in_element0.push_back(nodes[0]);
        nodes_in_element0.push_back(nodes[1]);
        nodes_in_element0.push_back(nodes[2]);
        nodes_in_element0.push_back(nodes[3]);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(nodes[1]);
        nodes_in_element1.push_back(nodes[4]);
        nodes_in_element1.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_in_element2;
        nodes_in_element2.push_back(nodes[4]);
        nodes_in_element2.push_back(nodes[5]);
        nodes_in_element2.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_in_element3;
        nodes_in_element3.push_back(nodes[7]);
        nodes_in_element3.push_back(nodes[0]);
        nodes_in_element3.push_back(nodes[8]);

        // Make elements
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
        elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
        elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));
        elements.push_back(new VertexElement<2,2>(3, nodes_in_element3));

        // Make mesh
        MutableVertexMesh<2,2> mesh(nodes, elements);
        mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);// Threshold distance set to ease calculations.

        // Node 6 and 8 are close to, but not overlapping, an edge of element 0
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(6)->rGetLocation(), 0), false);
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(8)->rGetLocation(), 0), false);

        // Move node 6 to the left so that it overlaps element 0
        ChastePoint<2> point = mesh.GetNode(6)->GetPoint();
        point.SetCoordinate(0u, 0.9);
        mesh.SetNode(6, point);

        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 0.9, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.5, 1e-4);

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(6)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(mesh.GetNode(6)->rGetLocation(), 0), 1u);

        // Move node 8 to the right so that it overlaps element 0
        point.SetCoordinate(0u, 0.1);
        mesh.SetNode(8, point);

        TS_ASSERT_DELTA(mesh.GetNode(8)->rGetLocation()[0], 0.1, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8)->rGetLocation()[1], 0.5, 1e-4);

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(8)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(mesh.GetNode(8)->rGetLocation(), 0), 3u);

        // Call method to update mesh in this situation
        mesh.ReMesh();//MoveOverlappingNodeOntoEdgeOfElement(mesh.GetNode(6), 0);

        // Check that node 6 has been moved onto the edge a new node has been created and both added to elements 0 amd 1
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10u);

        // Test locations of moved and new nodes
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.45, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(9)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(9)->rGetLocation()[1], 0.55, 1e-4);

        TS_ASSERT_DELTA(mesh.GetNode(8)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8)->rGetLocation()[1], 0.5, 1e-4);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(3), 9u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(4), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(5), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(6), 8u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);

        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(0), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(2), 9u);
        TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNodeGlobalIndex(0), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(3)->GetNodeGlobalIndex(2), 8u);

        // Test boundary property of nodes. All are boundary nodes except node 6.
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==6)
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

   /**
     * This tests the ReMesh method for preforming T3Swaps, In this test we read in a vertex mesh
     * that contains several nodes that are inside other elements.
     *
     *     ____
     *    _\ | /_
     * |\|  \|/  |/|
     * | \       /_|
     * | /       \ |
     * |/|__/\___|\|
     *     /__\
     *
     *
     *      |\                /|
     *      |_\ |          | /_|
     *      | / v          v \ |
     *  ____|/________________\|
     *  \  /|                  |
     *   \/ |                  |
     *   -->|                  |<--
     *      |                  | /\
     *      |__________________|/__\
     *      |\                /|
     *      |_\ ^          ^ /_|
     *      | / |          | \ |
     *      |/                \|
     *
     * Note: this also tests that boundary nodes are updated accordingly
     */

    void TestReMeshForT3Swaps() throw(Exception)
    {
        // This also tests IdentifySwapType

        // LoadMesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_remesh_T3");
        MutableVertexMesh<2,2> vertex_mesh;

        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        vertex_mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 38u);

        // assign boundary nodes \todo #1076 - once reading/writing of boundary elements is done
        // properly for vertex meshes this can be added to the .node file
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            vertex_mesh.GetNode(i)->SetAsBoundaryNode(true);
        }

        // Calls ReMesh to identify all T3 swaps (element overlaps) and perform them.
        vertex_mesh.ReMesh();


        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 48u);

        std::string dirname = "vertex_remeshing_mesh";
        std::string mesh_filename = "vertex_mesh_T3";

        // Save the mesh data using mesh writers
        VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh(vertex_mesh);


        //Test Moved Nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(14)->rGetLocation()[0], 0.55, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(14)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(15)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(15)->rGetLocation()[1], 0.5, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(16)->rGetLocation()[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(16)->rGetLocation()[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(17)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(17)->rGetLocation()[1], 0.45, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(23)->rGetLocation()[0], 4.95, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(23)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(27)->rGetLocation()[0], 5.25, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(27)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(29)->rGetLocation()[0], 6.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(29)->rGetLocation()[1], 0.85, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(32)->rGetLocation()[0], 5.05, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(32)->rGetLocation()[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(33)->rGetLocation()[0], 4.75, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(33)->rGetLocation()[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(37)->rGetLocation()[0], 4.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(37)->rGetLocation()[1], 0.15, 1e-4);

        //Test Added Nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(38)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(38)->rGetLocation()[1], 0.55, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(39)->rGetLocation()[0], 0.45, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(39)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(40)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(40)->rGetLocation()[1], 0.4, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(41)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(41)->rGetLocation()[1], 0.6, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(42)->rGetLocation()[0], 0.6, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(42)->rGetLocation()[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(43)->rGetLocation()[0], 0.4, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(43)->rGetLocation()[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(44)->rGetLocation()[0], 5.05, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(44)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(45)->rGetLocation()[0], 5.15, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(45)->rGetLocation()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(46)->rGetLocation()[0], 4.95, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(46)->rGetLocation()[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(47)->rGetLocation()[0], 4.85, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(47)->rGetLocation()[1], 1.0, 1e-4);

        // Test elements have correct nodes (1st Block)
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 39u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(4)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(5)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(6)->GetIndex(), 41u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(7)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(8)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(9)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(10)->GetIndex(), 43u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(11)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(12)->GetIndex(), 38u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(13)->GetIndex(), 17u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 38u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 10u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 39u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(3)->GetIndex(), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(0)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(1)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(2)->GetIndex(), 41u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(3)->GetIndex(), 15u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(0)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(1)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(2)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(3)->GetIndex(), 43u);


        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(0)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(1)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(2)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(3)->GetIndex(), 13u);


         // Test elements have correct nodes (2nd Block)
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNumNodes(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(0)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(1)->GetIndex(), 23u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(2)->GetIndex(), 44u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(3)->GetIndex(), 45u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(4)->GetIndex(), 27u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(5)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(6)->GetIndex(), 29u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(7)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(8)->GetIndex(), 32u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(9)->GetIndex(), 46u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(10)->GetIndex(), 47u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(11)->GetIndex(), 33u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(12)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(7)->GetNode(13)->GetIndex(), 37u);


        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 38u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 10u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 39u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(3)->GetIndex(), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(0)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(1)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(2)->GetIndex(), 41u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(3)->GetIndex(), 15u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(0)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(1)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(2)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(5)->GetNode(3)->GetIndex(), 43u);


        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(0)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(1)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(2)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(6)->GetNode(3)->GetIndex(), 13u);

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==15 || i==16 || i==23 || i==27 || i==32 || i==33)
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

    }

    void TestReMeshForRemovingVoids() throw(Exception)
    {

        /*
         * 3 elementswith a central void
         *  ______       _______
         * |     /|     |      /|
         * |___/| |     |_____/ |
         * |   \| | --> |     \ |
         * |_____\|     |______\|
         */

        // Make 7 nodes to assign to three elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.4, 0.5));
        nodes.push_back(new Node<2>(5, true, 0.55, 0.4));
        nodes.push_back(new Node<2>(6, true, 0.55, 0.6));
        nodes.push_back(new Node<2>(7, true, 0.0, 0.5));

        // Make three elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[7]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[6]);
        nodes_elem_1.push_back(nodes[5]);

        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[7]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[6]);
        nodes_elem_2.push_back(nodes[2]);
        nodes_elem_2.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        vertex_mesh.ReMesh(); // Edges too long so nothing happens

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        c_vector<double, 2>& new_location_0 = vertex_mesh.GetNode(5)->rGetModifiableLocation();
        new_location_0(1) = 0.51;

        c_vector<double, 2>& new_location_1 = vertex_mesh.GetNode(6)->rGetModifiableLocation();
        new_location_1(1) = 0.49;

        // T1 swap should now happen, removing the void
        vertex_mesh.ReMesh();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test merged node is in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-4);

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_4;
        expected_elements_containing_node_4.insert(0);
        expected_elements_containing_node_4.insert(1);
        expected_elements_containing_node_4.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->rGetContainingElementIndices(), expected_elements_containing_node_4);

        // Test elements have correct nodes
        // Note: nodes are renumbered void is removed and nodes reordered.
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 3u);
    }

    void TestReMeshForRemovingVoidsForCoverage() throw(Exception)
    {

        /*
         * 3 elementswith a central void
         *           _________
         *          |        /|
         *         |/|\___/| |
         *         |\|/   \| |
         *          |________\|
         *
         * The 2 central trinagles are voids
         */

        // Make 11 nodes to assign to three elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.45, 0.49));
        nodes.push_back(new Node<2>(5, true, 0.6, 0.5));
        nodes.push_back(new Node<2>(6, true, 0.45, 0.51));
        nodes.push_back(new Node<2>(7, true, 0.8, 0.5));
        nodes.push_back(new Node<2>(8, true, 0.9, 0.51));
        nodes.push_back(new Node<2>(9, true, 0.9, 0.49));
        nodes.push_back(new Node<2>(10, true, 1.0, 0.5));

        // Make 4 elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[10]);
        nodes_elem_0.push_back(nodes[9]);
        nodes_elem_0.push_back(nodes[7]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[4]);

        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[6]);
        nodes_elem_1.push_back(nodes[3]);

        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[7]);
        nodes_elem_2.push_back(nodes[8]);
        nodes_elem_2.push_back(nodes[10]);
        nodes_elem_2.push_back(nodes[2]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[6]);
        nodes_elem_2.push_back(nodes[5]);

        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[10]);
        nodes_elem_3.push_back(nodes[8]);
        nodes_elem_3.push_back(nodes[9]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 11u);

        // Call IdentifySwapType on nodes 6 and 4  (ordering for coverage)
        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(6),vertex_mesh.GetNode(4),map);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 9u);

        // Test merged node is in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-4);

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_4;
        expected_elements_containing_node_4.insert(0);
        expected_elements_containing_node_4.insert(1);
        expected_elements_containing_node_4.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->rGetContainingElementIndices(), expected_elements_containing_node_4);

        // Call IdentifySwapType on nodes 6 and 7  (originaly nodes 8 and 9)
        VertexElementMap map_2(vertex_mesh.GetNumElements());
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(6), vertex_mesh.GetNode(7), map_2),
                              "Triangular element next to triangular void, not implemented yet.");
    }
};

#endif /*TESTMUTABLEVERTEXMESHREMESH_HPP_*/
