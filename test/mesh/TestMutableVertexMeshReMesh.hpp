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

#include "Warnings.hpp"

class TestMutableVertexMeshReMesh : public CxxTest::TestSuite
{
public:

	/*
	 * This tests both PerformNodeMerge and IdentifySwapType
	 *
	 *
	 *       /|
	 *      / |
	 *     /  |
	 *    /   |
	 *   /    |
	 *  /     |
	 *  --xx--
	 *
	 *  The nodes marked with an x are merged
	 *
	 */
	void TestPerformNodeMerge() throw(Exception)
	{
		// Create nodes
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
		nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
		nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
		nodes.push_back(new Node<2>(3, true, 0.4, 0.0));
		nodes.push_back(new Node<2>(4, true, 0.6, 0.0));

		// Create two elements containing nodes
		std::vector<Node<2>*> nodes_elem_0;
		nodes_elem_0.push_back(nodes[0]);
		nodes_elem_0.push_back(nodes[3]);
		nodes_elem_0.push_back(nodes[4]);
		nodes_elem_0.push_back(nodes[1]);
		nodes_elem_0.push_back(nodes[2]);

		std::vector<VertexElement<2,2>*> vertex_elements;
		vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));

		// Make a vertex mesh
		MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

		// Test mesh has the correct numbers of elements and nodes
		TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
		TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);

		// Test the correct nodes are boundary nodes
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(0)->IsBoundaryNode(), true);
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(1)->IsBoundaryNode(), true);
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(2)->IsBoundaryNode(), true);
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(3)->IsBoundaryNode(), true);
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(4)->IsBoundaryNode(), true);


		// Merge nodes 3 and 4
		VertexElementMap map(vertex_mesh.GetNumElements());
		vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(3), vertex_mesh.GetNode(4), map);

		// Test that the mesh is correctly updated
		TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
		TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

		// Test the correct nodes are boundary nodes
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(0)->IsBoundaryNode(), true);
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(1)->IsBoundaryNode(), true);
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(2)->IsBoundaryNode(), true);
		TS_ASSERT_EQUALS(vertex_mesh.GetNode(3)->IsBoundaryNode(), true);

		// Test merged node is in the correct place
		TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[0], 0.5, 1e-3);
		TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[1], 0.0, 1e-3);

		// Test elements have correct nodes
		TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
		TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
		TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
		TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 1u);
		TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 2u);

		// Test Area and Perimeter of element
		TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.5, 1e-6);
		TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 2+sqrt(2), 1e-6);
	}

    /*
     * This test provides coverage of the case in which, when the elements
     * previously containing the high-index node are updated to contain the
     * low-index node, at least one of these elements did not already contain
     * the low-index node.
     *
     *   -------x-x-----
     *  |       |       |
     *  |       |       |
     *  |       |       |
     *   ------- -------
     *
     *  The nodes marked with an x are merged
     *
     * \TODO this could be merged with the earlier test to shorten the test file.
     * \TODO I think this should be a T1Swap see #1263
     *
     */
    void TestPerformNodeMergeWhenLowIndexNodeMustBeAddedToElement() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 2.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 1.01, 1.0));
        nodes.push_back(new Node<2>(5, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(6, true, 0.0, 2.0));

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

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

    }

    // This tests both PerformT1Swap and IdentifySwapType
    void TestPerformT1Swap() throw(Exception)
    {

        /* Make 6 nodes to assign to 4 elements.
         *
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         *
         */

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 1.2+0.2*sqrt(41.0), 1e-6);

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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test T1Swap Location tracking
        std::vector< c_vector<double, 2> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);

        // Test T1Swap Location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);

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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);

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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

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
        /* Make 6 nodes to assign to 3 elements with one non boundary nodes
         *
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         *
         */

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);

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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);

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
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test boundary property of nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestPerformT1SwapExceptions() throw(Exception)
    {
        // Make 6 nodes to assign to four elements where nodes 4 and 5 are the same position
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.5));

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

        // Perform a T1 swap on nodes 4 and 5 by using identify swap type
        VertexElementMap map(vertex_mesh.GetNumElements());

        // Call remesh
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), map), "Nodes are too close together, this shouldn't happen");
    }

    /*
     * This tests that T1Swaps rearrange to form a Triangular element for a T2 Swap
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
        //vertexElement<2,2>* p_element_4 = vertex_mesh.GetElement(4);
        //vertex_mesh.PerformT2Swap(*p_element_4);
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
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(),4u); // Elements are deleted not just marked for deletion, as calling ReMesh().
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(),6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(),6u); // Nodes are deleted not just marked for deletion, as calling ReMesh().

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
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.4, 0.25));
        nodes.push_back(new Node<2>(4, false, 0.6, 0.25));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.3));

        /*
        *  Make three trapezium elements with a central triangular element out of these nodes
        *
        *      /|\
        *     / | \
        *    /  |  \
        *   /2 /_\ 1\   Triangular element is element zero
        *  / _/   \_ \
        * /_/___3___\_\
        *
        */


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

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Perform a T2 swap on the middle triangle element
        VertexElement<2,2>* p_element_0 = vertex_mesh.GetElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 7u);

        for (unsigned j=1; j<4; j++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumNodes(), 3u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(0)->GetIndex(), j%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(1)->GetIndex(), (j+1)%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(2)->GetIndex(), 6u);
        }

        // Test boundary property of nodes. All are boundary nodes except node 3.
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==3)
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestPerformT2SwapWithBoundaryNodes() throw(Exception)
    {
    	{
			// Make 6 nodes to assign to three elements
			std::vector<Node<2>*> nodes;
			nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
			nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
			nodes.push_back(new Node<2>(2, true, 0.5, 0.5));
			nodes.push_back(new Node<2>(3, true, 0.4, 0.25));
			nodes.push_back(new Node<2>(4, true, 0.6, 0.25));
			nodes.push_back(new Node<2>(5, false, 0.5, 0.3));

			/*
			*  Make two trapezium elements with a central triangular element out of these nodes
			*
			*      /|\
			*     / | \
			*    /  |  \
			*   /2 /_\ 1\   Triangular element is element zero
			*  / _/   \_ \
			* /_/       \_\
			*
			*/


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


			std::vector<VertexElement<2,2>*> vertex_elements;
			vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
			vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
			vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));

			// Make a vertex mesh
			MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

			TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
			TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

			// Perform a T2 swap on the middle triangle element
			VertexElement<2,2>* p_element_0 = vertex_mesh.GetElement(0);
			vertex_mesh.PerformT2Swap(*p_element_0);

			TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
			TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

			TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 3u);
			TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 7u);

			for (unsigned j=1; j<3; j++)
			{
				TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumNodes(), 3u);
				TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(0)->GetIndex(), j%3);
				TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(1)->GetIndex(), (j+1)%3);
				TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(2)->GetIndex(), 6u);
			}

			// Test boundary property of nodes. All are boundary nodes.
			for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
			{
				TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
			}
    	}

    	{
			// Make 5 nodes to assign to two elements
			std::vector<Node<2>*> nodes;
			nodes.push_back(new Node<2>(0, true, 1.0, 0.0));
			nodes.push_back(new Node<2>(1, true, 0.5, 0.5));
			nodes.push_back(new Node<2>(2, true, 0.4, 0.25));
			nodes.push_back(new Node<2>(3, true, 0.6, 0.25));
			nodes.push_back(new Node<2>(4, true, 0.5, 0.3));

			/*
			*  Make one trapezium element with a central triangular element out of these nodes
			*
			*       |\
			*       | \
			*       |  \
			*      /_\ 1\   Triangular element is element zero
			*         \_ \
			*           \_\
			*
			*/


			// Triangle element
			std::vector<Node<2>*> nodes_elem_0;
			nodes_elem_0.push_back(nodes[2]);
			nodes_elem_0.push_back(nodes[3]);
			nodes_elem_0.push_back(nodes[4]);

			// Trapezium
			std::vector<Node<2>*> nodes_elem_1;
			nodes_elem_1.push_back(nodes[0]);
			nodes_elem_1.push_back(nodes[1]);
			nodes_elem_1.push_back(nodes[4]);
			nodes_elem_1.push_back(nodes[3]);

			std::vector<VertexElement<2,2>*> vertex_elements;
			vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
			vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

			// Make a vertex mesh
			MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

			TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
			TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);

			// Perform a T2 swap on the middle triangle element
			VertexElement<2,2>* p_element_0 = vertex_mesh.GetElement(0);
			vertex_mesh.PerformT2Swap(*p_element_0);

			TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
			TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 3u);

			TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 2u);
			TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 6u);


		    TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
			TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 0u);
			TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 1u);
			TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 5u);

			// Test boundary property of nodes. All are boundary nodes.
			for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
			{
				TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
			}
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
     * This tests the ReMesh method for performing T2Swaps (element removal).
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
     * This tests the ReMesh method for performing T1Swaps, both internally and on the boundary.
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

        // Load mesh
        VertexMeshReader<2,2> mesh_reader("notforrelease_cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1");
        MutableVertexMesh<2,2> vertex_mesh;

        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 22u);

        // Calls ReMesh to identify all T1 swaps and perform them.
        vertex_mesh.ReMesh();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 22u);

        std::string dirname = "TestVertexMeshReMesh";
        std::string mesh_filename = "vertex_remesh_T1";

        // Save the mesh data using mesh writers
        VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh(vertex_mesh);


        //Check the positions are updated correctly
        OutputFileHandler handler("TestVertexMeshReMesh", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T1.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T1.cell";

        // To ignore the provenance data we only go as far as
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file1 + " notforrelease_cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1_after_remesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file2 + " notforrelease_cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1_after_remesh.cell").c_str()), 0);
    }

    void TestReMeshExceptions() throw(Exception)
    {
        // This also tests IdentifySwapType
    	{
    		/*
			 *   ______
			 *  |     /|
			 *  |    / |
			 *  |   x  |
			 *  |  x   |
			 *  | /    |
			 *  |/_____|
			 *
			 *  The nodes marked with an x are merged
			 *
			 */

    		// Create nodes
    		std::vector<Node<2>*> nodes;
    		nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
    		nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
    		nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
    		nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
    		nodes.push_back(new Node<2>(4, false, 0.49, 0.49));
    		nodes.push_back(new Node<2>(5, false, 0.51, 0.51));

    		// Create two elements containing nodes
    		std::vector<Node<2>*> nodes_elem_0;
    		nodes_elem_0.push_back(nodes[0]);
    		nodes_elem_0.push_back(nodes[1]);
    		nodes_elem_0.push_back(nodes[2]);
    		nodes_elem_0.push_back(nodes[5]);
    		nodes_elem_0.push_back(nodes[4]);

    		std::vector<Node<2>*> nodes_elem_1;
    		nodes_elem_1.push_back(nodes[0]);
    		nodes_elem_1.push_back(nodes[4]);
    		nodes_elem_1.push_back(nodes[5]);
    		nodes_elem_1.push_back(nodes[2]);
    		nodes_elem_1.push_back(nodes[3]);

    		std::vector<VertexElement<2,2>*> vertex_elements;
    		vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
    		vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

    		// Make a vertex mesh
    		MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
    		vertex_mesh.SetCellRearrangementThreshold(0.1);

    		// Try to Merge nodes 4 and 5
    		TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "There are non boundary nodes contained in only in 2 elements something has gone wrong.");
    	}

    	{

			/* Create 4 joined elements.
			 *  _________
			 * |    |   /
			 * |    |  /
			 * |____|_/
			 * |    | \
			 * |    |  \
			 * |____|___\
			 *
			 */

			// Create some nodes all boundary nodes except the central one.
			std::vector<Node<2>*> nodes;
			nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
			nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
			nodes.push_back(new Node<2>(2, true, 2.0, 0.0));
			nodes.push_back(new Node<2>(3, true, 1.1, 1.0));
			nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
			nodes.push_back(new Node<2>(5, true, 1.0, 2.0));
			nodes.push_back(new Node<2>(6, true, 0.0, 2.0));
			nodes.push_back(new Node<2>(7, true, 0.0, 1.0));
			nodes.push_back(new Node<2>(8, false, 1.0, 1.0));

			std::vector<Node<2>*> nodes_in_element0;
			nodes_in_element0.push_back(nodes[0]);
			nodes_in_element0.push_back(nodes[1]);
			nodes_in_element0.push_back(nodes[8]);
			nodes_in_element0.push_back(nodes[7]);

			std::vector<Node<2>*> nodes_in_element1;
			nodes_in_element1.push_back(nodes[1]);
			nodes_in_element1.push_back(nodes[2]);
			nodes_in_element1.push_back(nodes[3]);
			nodes_in_element1.push_back(nodes[8]);

			std::vector<Node<2>*> nodes_in_element2;
			nodes_in_element2.push_back(nodes[8]);
			nodes_in_element2.push_back(nodes[3]);
			nodes_in_element2.push_back(nodes[4]);
			nodes_in_element2.push_back(nodes[5]);

			std::vector<Node<2>*> nodes_in_element3;
			nodes_in_element3.push_back(nodes[7]);
			nodes_in_element3.push_back(nodes[8]);
			nodes_in_element3.push_back(nodes[5]);
			nodes_in_element3.push_back(nodes[6]);

			std::vector<VertexElement<2,2>* > elements;
			elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
			elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
			elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));
			elements.push_back(new VertexElement<2,2>(3, nodes_in_element3));

			// Create mesh
			MutableVertexMesh<2,2> vertex_mesh(nodes, elements);
			vertex_mesh.SetCellRearrangementThreshold(0.2); // so will call IdentifySwapType()

			TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 9u);
			TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);

			// Call remesh

	        // Attempt to merge nodes 7 and 8.
	        TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "A node is contained in more than three elements");
    	}

    	{
			/*
    	     *             ______
    	     *            /      |
    	     *           /       |
    	     *          /        |
    	     *   -----xx---------|
    	     *  |       \        |
    	     *  |        \       |
    	     *  |_________\______|
    	     *
    	     *  The nodes marked with an x are to be merged.
    	     *  This will throw an exception in IdentifySwapType.
    	     */

			// Create nodes, all boundary nodes
			std::vector<Node<2>*> nodes;
			nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
			nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
			nodes.push_back(new Node<2>(2, true, 3.0, 0.0));
			nodes.push_back(new Node<2>(3, true, 3.0, 1.0));
			nodes.push_back(new Node<2>(4, true, 3.0, 2.0));
			nodes.push_back(new Node<2>(5, true, 2.0, 2.0));
			nodes.push_back(new Node<2>(6, true, 1.0, 1.0));
			nodes.push_back(new Node<2>(7, true, 0.99, 1.0));
			nodes.push_back(new Node<2>(8, true, 0.0, 1.0));

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

			// Attempt to Merge nodes 6 and 7
			VertexElementMap map(vertex_mesh.GetNumElements());
			TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(6), vertex_mesh.GetNode(7), map), "There is a boundary node contained in three elements something has gone wrong.");
		}

    	{
			/* Make 7 nodes to assign to 2 elements with 6 boundary nodes
			 *
			 * |\   /|
			 * | \ / |
			 * |  x  |
			 * |  x  |
			 * |__|__|
			 *
			 */

			std::vector<Node<2>*> nodes;
			nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
			nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
			nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
			nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
			nodes.push_back(new Node<2>(4, false, 0.5, 0.49));
			nodes.push_back(new Node<2>(5, true, 0.5, 0.51));
			nodes.push_back(new Node<2>(6, true, 0.5, 0.0));

			// Make two elements out of these nodes
			std::vector<Node<2>*> nodes_elem_0;
			nodes_elem_0.push_back(nodes[1]);
			nodes_elem_0.push_back(nodes[2]);
			nodes_elem_0.push_back(nodes[5]);
			nodes_elem_0.push_back(nodes[4]);
			nodes_elem_0.push_back(nodes[6]);

			std::vector<Node<2>*> nodes_elem_1;
			nodes_elem_1.push_back(nodes[0]);
			nodes_elem_1.push_back(nodes[6]);
			nodes_elem_1.push_back(nodes[4]);
			nodes_elem_1.push_back(nodes[5]);
			nodes_elem_1.push_back(nodes[3]);

			std::vector<VertexElement<2,2>*> vertex_elements;
			vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
			vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

			// Make a vertex mesh
			MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
    		vertex_mesh.SetCellRearrangementThreshold(0.1);

			// Attempt to Merge nodes 4 and 5
			TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "There is a non boundary node contained only in 2 elements something has gone wrong.");
		}

    	{
			/* Make 8 nodes to assign to 3 elements with 5 boundary nodes
			 *  __x__
			 * |\   /|
			 * | \ / |
			 * |  x  |
			 * |  x  |
			 * |__|__|
			 *
			 * Not the extra node on the top is to stop the element only having 3 nodes,
			 * otherwise wont go into IdentifySwapType().
			 *
			 */

			std::vector<Node<2>*> nodes;
			nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
			nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
			nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
			nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
			nodes.push_back(new Node<2>(4, false, 0.5, 0.49));
			nodes.push_back(new Node<2>(5, false, 0.5, 0.51));
			nodes.push_back(new Node<2>(6, true, 0.5, 0.0));
			nodes.push_back(new Node<2>(7, true, 0.5, 1.0));

			// Make two elements out of these nodes
			std::vector<Node<2>*> nodes_elem_0;
			nodes_elem_0.push_back(nodes[1]);
			nodes_elem_0.push_back(nodes[2]);
			nodes_elem_0.push_back(nodes[5]);
			nodes_elem_0.push_back(nodes[4]);
			nodes_elem_0.push_back(nodes[6]);

			std::vector<Node<2>*> nodes_elem_1;
			nodes_elem_1.push_back(nodes[0]);
			nodes_elem_1.push_back(nodes[6]);
			nodes_elem_1.push_back(nodes[4]);
			nodes_elem_1.push_back(nodes[5]);
			nodes_elem_1.push_back(nodes[3]);

			std::vector<Node<2>*> nodes_elem_2;
			nodes_elem_2.push_back(nodes[2]);
			nodes_elem_2.push_back(nodes[3]);
			nodes_elem_2.push_back(nodes[5]);
			nodes_elem_2.push_back(nodes[7]);

			std::vector<VertexElement<2,2>*> vertex_elements;
			vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
			vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
			vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));

			// Make a vertex mesh
			MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
			vertex_mesh.SetCellRearrangementThreshold(0.1);

			// Attempt to Merge nodes 4 and 5
			TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "There are non boundary nodes contained only in 2 elements something has gone wrong.");
		}
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

    void TestPerformT3Swap()
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

        // Save the mesh data using mesh writers
				std::string dirname = "TempyTempy";
                std::string mesh_filename = "vertex_remesh_T3";
                VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
                mesh_writer.WriteFilesUsingMesh(mesh);

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

        // Test T3Swap Location tracking
        std::vector< c_vector<double, 2> > t3_locations = mesh.GetLocationsOfT3Swaps();
        TS_ASSERT_EQUALS(t3_locations.size(), 2u);
        TS_ASSERT_DELTA(t3_locations[0][0], 1.0, 1e-6);
        TS_ASSERT_DELTA(t3_locations[0][1], 0.5, 1e-6);
        TS_ASSERT_DELTA(t3_locations[1][0], 0.0, 1e-6);
        TS_ASSERT_DELTA(t3_locations[1][1], 0.5, 1e-6);

        // Test T1Swap Location clearing
        mesh.ClearLocationsOfT3Swaps();
        t3_locations = mesh.GetLocationsOfT3Swaps();
        TS_ASSERT_EQUALS(t3_locations.size(), 0u);
    }

    void TestPerformT3SwapExceptions() throw(Exception)
    {
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
        }

        {
            /*
             * Make a small mesh consisting of 3 elements: and move the top left node to intersect the central horizontal edge.
             *   __
             *  |  |\
             *  |__| \
             *   \ | /
             *    \|/
             */

            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
            nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
            nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
            nodes.push_back(new Node<2>(4, true, 1.0, -1.0));
            nodes.push_back(new Node<2>(5, true, 2.0, 0.0));

            std::vector<Node<2>*> nodes_in_element0;
            nodes_in_element0.push_back(nodes[0]);
            nodes_in_element0.push_back(nodes[1]);
            nodes_in_element0.push_back(nodes[2]);
            nodes_in_element0.push_back(nodes[3]);

            std::vector<Node<2>*> nodes_in_element1;
            nodes_in_element1.push_back(nodes[4]);
            nodes_in_element1.push_back(nodes[5]);
            nodes_in_element1.push_back(nodes[2]);

            std::vector<Node<2>*> nodes_in_element2;
            nodes_in_element2.push_back(nodes[0]);
            nodes_in_element2.push_back(nodes[4]);
            nodes_in_element2.push_back(nodes[1]);

            // Make elements
            std::vector<VertexElement<2,2>*> elements;
            elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
            elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
            elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));

            // Make mesh
            MutableVertexMesh<2,2> mesh(nodes, elements);

            // Move node 3  so that it overlaps element 2 across an internal edge
            ChastePoint<2> point = mesh.GetNode(3)->GetPoint();
            point.SetCoordinate(0u, 0.5);
            point.SetCoordinate(1u, -0.1);
            mesh.SetNode(3, point);

            TS_ASSERT_DELTA(mesh.GetNode(3)->rGetLocation()[0], 0.5, 1e-4);
            TS_ASSERT_DELTA(mesh.GetNode(3)->rGetLocation()[1], -0.1, 1e-4);

            TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(3)->rGetLocation(), 2), true);
            TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(mesh.GetNode(3)->rGetLocation(), 0), 2u);

            // Call ReMesh which in turn calls PerformT3Swap
            TS_ASSERT_THROWS_THIS(mesh.ReMesh(), "A boundary node has intersected a non boundary edge, this is because the boundary element has become concave you need to rerun the simulation with a smaller time step to prevent this.");
        }
    }


    void TestT3SwapOnSmallEdge()
    {
        /* Create a rectangular and a two triangular elements to test intersecting on an
         * edge that is too small
         *
         *
         *  ______  /|
         * |      |/_| <---
         * |______|\ |
         *          \|
         *
         *
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 0.1));
        nodes.push_back(new Node<2>(3, true, 0.0, 0.1));
        nodes.push_back(new Node<2>(4, true, 0.99, 0.05));
        nodes.push_back(new Node<2>(5, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, true, 2.0, 0.1));
        nodes.push_back(new Node<2>(7, true, 2.0, 0.0));

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
        nodes_in_element2.push_back(nodes[7]);
        nodes_in_element2.push_back(nodes[4]);
        nodes_in_element2.push_back(nodes[6]);

        // Make elements
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
        elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
        elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));


        // Make mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1/1.5);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);

        // Call PerformT3Swap
        // Note we don't call ReMesh as this would also perform T1Swaps.
        vertex_mesh.PerformT3Swap(vertex_mesh.GetNode(4), 0u);

        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Trying to merge a node onto an edge which is too small.");
        Warnings::QuietDestroy();

        // Check that node 4 has been moved onto the edge and a new node has been created and both added to elements 0 and 1
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);

        // Test locations of moved and new nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.05, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[1], -0.05, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(9)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(9)->rGetLocation()[1], 0.15, 1e-8);


        // Test locations of edges of intersected edge
        TS_ASSERT_DELTA(vertex_mesh.GetNode(1)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(1)->rGetLocation()[1], -0.15, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(2)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(2)->rGetLocation()[1], 0.25, 1e-8);


        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(4), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(5), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(6), 3u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), 8u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        // Test boundary property of nodes. All are boundary nodes except 4.
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==4)
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }


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
        nodes.push_back(new Node<2>(7, true, -1.0, 0.0));
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
        MutableVertexMesh<2,2> vertex_mesh(nodes, elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);// Threshold distance set to ease calculations.

        // Node 6 and 8 are close to, but not overlapping, an edge of element 0
        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(6)->rGetLocation(), 0), false);
        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(8)->rGetLocation(), 0), false);

        // Move node 6 to the left so that it overlaps element 0
        ChastePoint<2> point = vertex_mesh.GetNode(6)->GetPoint();
        point.SetCoordinate(0u, 0.9);
        vertex_mesh.SetNode(6, point);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 0.9, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.5, 1e-4);

        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(6)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetLocalIndexForElementEdgeClosestToPoint(vertex_mesh.GetNode(6)->rGetLocation(), 0), 1u);

        // Move node 8 to the right so that it overlaps element 0
        point.SetCoordinate(0u, 0.1);
        vertex_mesh.SetNode(8, point);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[0], 0.1, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[1], 0.5, 1e-4);

        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(8)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetLocalIndexForElementEdgeClosestToPoint(vertex_mesh.GetNode(8)->rGetLocation(), 0), 3u);

        // Call method to update vertex_mesh in this situation
        vertex_mesh.ReMesh();//MoveOverlappingNodeOntoEdgeOfElement(vertex_mesh.GetNode(6), 0);

        // Check that node 6 has been moved onto the edge a new node has been created and both added to elements 0 amd 1
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);

        // Test locations of moved and new nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.45, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(9)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(9)->rGetLocation()[1], 0.55, 1e-4);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[1], 0.5, 1e-4);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(4), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(5), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(6), 8u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(0), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(2), 8u);

        // Test boundary property of nodes. All are boundary nodes except node 6.
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==6)
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }


    void TestT3SwapForNeighboringElementsWithTwoCommonNodes()
    {
        /*
         * Make a small mesh consisting of 4 elements:
         * a square and 2 quadrilaterals and a triangles.
         *        _____
         *       |     |
         *     /\|     |/|\
         *    /__|_____|_|_\
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
        nodes.push_back(new Node<2>(7, true, -1.0, 0.0));
        nodes.push_back(new Node<2>(8, true, -0.5, 0.5));
        nodes.push_back(new Node<2>(9, true, 1.0, 0.25));
        nodes.push_back(new Node<2>(10, true, 0.0, 0.25));

        std::vector<Node<2>*> nodes_in_element0;
        nodes_in_element0.push_back(nodes[0]);
        nodes_in_element0.push_back(nodes[1]);
        nodes_in_element0.push_back(nodes[9]);
        nodes_in_element0.push_back(nodes[2]);
        nodes_in_element0.push_back(nodes[3]);
        nodes_in_element0.push_back(nodes[10]);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(nodes[1]);
        nodes_in_element1.push_back(nodes[4]);
        nodes_in_element1.push_back(nodes[6]);
        nodes_in_element1.push_back(nodes[9]);

        std::vector<Node<2>*> nodes_in_element2;
        nodes_in_element2.push_back(nodes[4]);
        nodes_in_element2.push_back(nodes[5]);
        nodes_in_element2.push_back(nodes[6]);

        std::vector<Node<2>*> nodes_in_element3;
        nodes_in_element3.push_back(nodes[7]);
        nodes_in_element3.push_back(nodes[0]);
        nodes_in_element3.push_back(nodes[10]);
        nodes_in_element3.push_back(nodes[8]);

        // Make elements
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
        elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
        elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));
        elements.push_back(new VertexElement<2,2>(3, nodes_in_element3));

        // Make mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);// Threshold distance set to ease calculations.

        // Node 6 and 8 are close to, but not overlapping, an edge of element 0
        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(6)->rGetLocation(), 0), false);
        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(8)->rGetLocation(), 0), false);

        // Move node 6 to the left so that it overlaps element 0
        ChastePoint<2> point = vertex_mesh.GetNode(6)->GetPoint();
        point.SetCoordinate(0u, 0.9);
        vertex_mesh.SetNode(6, point);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 0.9, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.5, 1e-4);

        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(6)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetLocalIndexForElementEdgeClosestToPoint(vertex_mesh.GetNode(6)->rGetLocation(), 0), 2u);

        // Move node 8 to the right so that it overlaps element 0
        point.SetCoordinate(0u, 0.1);
        vertex_mesh.SetNode(8, point);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[0], 0.1, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[1], 0.5, 1e-4);

        TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(8)->rGetLocation(), 0), true);
        TS_ASSERT_EQUALS(vertex_mesh.GetLocalIndexForElementEdgeClosestToPoint(vertex_mesh.GetNode(8)->rGetLocation(), 0), 4u);

        // Call method to update vertex_mesh in this situation
        vertex_mesh.ReMesh();//MoveOverlappingNodeOntoEdgeOfElement(vertex_mesh.GetNode(6), 0);

        // Check that node 6 has been moved onto the edge a new node has been created and both added to elements 0 and 1
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);


        // Test locations of moved and new nodes (9 is the next free node when 6 is merged)
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.45, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(9)->rGetLocation()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(9)->rGetLocation()[1], 0.55, 1e-4);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[1], 0.5, 1e-4);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(4), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(5), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(6), 8u);


        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(0), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(1), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(2), 8u);

        // Test boundary property of nodes. All are boundary nodes except node 6.
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = true;
            if (i==6)
            {
                expected_boundary_node = false;
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

   /**
     * This tests the ReMesh method for performing T3Swaps, In this test we read in a vertex mesh
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
     *
     *
     *      |\  |          |  /|
     *      |_\ v          v /_|
     *   ___|_/______________\_|
     *   \  |                  |
     *    \/|                  |
     *   -->|                  |<--
     *      |                  |/\
     *      |__________________|__\
     *      |_\              /_|
     *      | / ^         ^  \ |
     *      |/  |         |   \|
     *
     * Note: this also tests that boundary nodes are updated accordingly
     */

    void TestReMeshForT3Swaps() throw(Exception)
    {
        // This also tests IdentifySwapType

        // Load mesh
        VertexMeshReader<2,2> mesh_reader("notforrelease_cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T3");
        MutableVertexMesh<2,2> vertex_mesh;

        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        vertex_mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);// Threshold distance set to ease calculations.

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 29u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 64u);

        // Calls ReMesh to identify all T3 swaps (element overlaps) and perform them.
        vertex_mesh.ReMesh();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 29u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 72u);

        // Save the mesh data using mesh writers
        std::string dirname = "TestVertexMeshReMesh";
        std::string mesh_filename = "vertex_remesh_T3";
        VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh(vertex_mesh);

        //Check the positions are updated correctly
        OutputFileHandler handler("TestVertexMeshReMesh", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T3.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T3.cell";

        // To ignore the provenance data we only go as far as
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file1 + " notforrelease_cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T3_after_remesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file2 + " notforrelease_cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T3_after_remesh.cell").c_str()), 0);
    }

    void TestReMeshForRemovingVoids() throw(Exception)
    {

        /*
         * 3 elements with a central void
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
         * 3 elements with a central void
         *           _________
         *          |        /|
         *          |/|\___/| |
         *          |\|/   \| |
         *          |________\|
         *
         * The 2 central triangles are voids
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


    void TestT3SwapForRemovingVoids() throw(Exception)
        {

            /*
             * 3 elements with a central void central node moved to overlap riight element
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
            nodes.push_back(new Node<2>(5, true, 0.5, 0.4));
            nodes.push_back(new Node<2>(6, true, 0.5, 0.6));
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

            //vertex_mesh.SetCellRearrangementThreshold(0.1);

            TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
            TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

            vertex_mesh.ReMesh(); // Edges too long so nothing happens

            TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
            TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);


            // Move node 4 to inside of element 1
            c_vector<double, 2>& new_location = vertex_mesh.GetNode(4)->rGetModifiableLocation();
            new_location(0) = 0.6;

            TS_ASSERT_EQUALS(vertex_mesh.ElementIncludesPoint(vertex_mesh.GetNode(4)->rGetLocation(), 1), true);
            TS_ASSERT_EQUALS(vertex_mesh.GetLocalIndexForElementEdgeClosestToPoint(vertex_mesh.GetNode(4)->rGetLocation(), 1), 2u);


            // T3 swap should now happen, removing the void
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


    void TestT3SwapWithConcaveElements()
       {
           /*
            * Make a small mesh consisting of 2 elements:
            * a rectangle sat upon an l shape (so concave)
            *     ______       _______
            *    |____|/      |_ _ | /
            *    | |____  --> | |_\|/_
            *    |______|     |_______|
            *
            */

           // Make all nodes boundary nodes
           std::vector<Node<2>*> nodes;
           nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
           nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
           nodes.push_back(new Node<2>(2, true, 2.0, 1.0));
           nodes.push_back(new Node<2>(3, true, 1.0, 1.0));
           nodes.push_back(new Node<2>(4, true, 1.0, 2.0));
           nodes.push_back(new Node<2>(5, true, 0.0, 2.0));
           nodes.push_back(new Node<2>(6, true, 2.0, 2.0));
           nodes.push_back(new Node<2>(7, true, 2.0, 3.0));
           nodes.push_back(new Node<2>(8, true, 0.0, 3.0));
           nodes.push_back(new Node<2>(9, true, 3.0, 3.0));

           std::vector<Node<2>*> nodes_in_element0;
           nodes_in_element0.push_back(nodes[0]);
           nodes_in_element0.push_back(nodes[1]);
           nodes_in_element0.push_back(nodes[2]);
           nodes_in_element0.push_back(nodes[3]);
           nodes_in_element0.push_back(nodes[4]);
           nodes_in_element0.push_back(nodes[5]);

           std::vector<Node<2>*> nodes_in_element1;
           nodes_in_element1.push_back(nodes[5]);
           nodes_in_element1.push_back(nodes[4]);
           nodes_in_element1.push_back(nodes[6]);
           nodes_in_element1.push_back(nodes[7]);
           nodes_in_element1.push_back(nodes[8]);

           std::vector<Node<2>*> nodes_in_element2;
           nodes_in_element2.push_back(nodes[7]);
           nodes_in_element2.push_back(nodes[6]);
           nodes_in_element2.push_back(nodes[9]);

           // Make elements
           std::vector<VertexElement<2,2>*> elements;
           elements.push_back(new VertexElement<2,2>(0, nodes_in_element0));
           elements.push_back(new VertexElement<2,2>(1, nodes_in_element1));
           elements.push_back(new VertexElement<2,2>(2, nodes_in_element2));

           // Make mesh
           MutableVertexMesh<2,2> mesh(nodes, elements);
           mesh.SetCellRearrangementThreshold(0.1*1.0/1.5);// Threshold distance set to ease calculations.

           // Move node 6  so that it overlaps element 0
           ChastePoint<2> point = mesh.GetNode(6)->GetPoint();
           point.SetCoordinate(0u, 1.5);
           point.SetCoordinate(1u, 0.9);
           mesh.SetNode(6, point);

           TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 1.5, 1e-4);
           TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 0.9, 1e-4);

           TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(mesh.GetNode(6)->rGetLocation(), 0), true);
           TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(mesh.GetNode(6)->rGetLocation(), 0), 2u);

           // Call method to update mesh in this situation
           mesh.ReMesh();//MoveOverlappingNodeOntoEdgeOfElement(mesh.GetNode(6), 0);

           // Check that node 6 has been moved onto the edge a new node has been created and both added to elements 0 and 1
           TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);
           TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);


           // Test locations of moved and new nodes (10 and 11 are the next free node when 6 is merged)
           TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 1.5, 1e-4);
           TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 1.0, 1e-4);
           TS_ASSERT_DELTA(mesh.GetNode(10)->rGetLocation()[0], 1.6, 1e-4);
           TS_ASSERT_DELTA(mesh.GetNode(10)->rGetLocation()[1], 1.0, 1e-4);
           TS_ASSERT_DELTA(mesh.GetNode(11)->rGetLocation()[0], 1.4, 1e-4);
           TS_ASSERT_DELTA(mesh.GetNode(11)->rGetLocation()[1], 1.0, 1e-4);


           // Test elements have correct nodes
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 9u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(2), 2u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(3), 10u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(4), 6u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(5), 11u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(6), 3u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(7), 4u);
           TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(8), 5u);


           TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 6u);
           TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 5u);
           TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
           TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 11u);
           TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(3), 6u);
           TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(4), 7u);
           TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(5), 8u);

           TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNumNodes(), 4u);
           TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(0), 7u);
           TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(1), 6u);
           TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(2), 10u);
           TS_ASSERT_EQUALS(mesh.GetElement(2)->GetNodeGlobalIndex(3), 9u);



           // Test boundary property of nodes. All are boundary nodes.
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

};

#endif /*TESTMUTABLEVERTEXMESHREMESH_HPP_*/
