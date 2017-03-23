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

#ifndef TESTMUTABLEVERTEXMESH33REMESH_HPP_
#define TESTMUTABLEVERTEXMESH33REMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <string>
#include "FileComparison.hpp"
#include "HexagonalPrism3dVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexMeshWriter.hpp"
#include "FakePetscSetup.hpp"

#define OUTPUT_NAME std::string("TestMutableVertexMesh33ReMesh")

class TestMutableVertexMesh33ReMesh : public CxxTest::TestSuite
{
public:
    void TestPerformNodeMerge() throw(Exception)
    {
        /*
         * Create a mesh comprising a single triangular element, as shown below.
         * We will test that the nodes marked with an x are merged correctly.
         *
         *      /|
         *     / |
         *    /  |
         *   /   |
         *  /    |
         *  --xx-
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 0.4, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.6, 0.0));

        const unsigned node_indices_elem_0[5] = { 0, 3, 4, 1, 2 };
        MonolayerVertexMeshGenerator builder(nodes, "NodeMerge");
        builder.BuildElementWith(5, node_indices_elem_0);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 7u);

        // Merge nodes 3 and 4
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(3), vertex_mesh.GetNode(4));
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");

        // Test the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFaces(), 6u);

        // Test the correct nodes are boundary nodes
        for (unsigned i = 0; i < vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        // Test the merged node is in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[0], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[1], 0.0, 1e-3);

        // Test the elements own the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 2u);

        // Test the element's area and perimeter are computed correctly
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 2 + sqrt(2.0) + 1.0, 1e-6);
    }

    void TestPerformNodeMergeWhenLowIndexNodeMustBeAddedToElement() throw(Exception)
    {
        /**
         * Create a mesh comprising two square elements, as shown below. We will test that the
         * nodes marked with an x are merged correctly. We will test node merging in the case
         * where, when the elements previously containing the high-index node are updated to
         * contain the low-index node, at least one of these elements did not already contain
         * the low-index node.
         *
         *   -----x-x---
         *  |     |     |
         *  |     |     |
         *   ----- -----
         *
         * \todo I think this should be a T1 swap (see #1263)
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.00, 0.00));
        nodes.push_back(new Node<3>(1, true, 1.00, 0.00));
        nodes.push_back(new Node<3>(2, true, 2.00, 0.00));
        nodes.push_back(new Node<3>(3, true, 2.00, 1.00));
        nodes.push_back(new Node<3>(4, true, 1.01, 1.00));
        nodes.push_back(new Node<3>(5, true, 1.00, 1.00));
        nodes.push_back(new Node<3>(6, true, 0.00, 2.00));

        unsigned node_indices_elem_0[4] = { 0, 1, 5, 6 };
        unsigned node_indices_elem_1[5] = { 1, 2, 3, 4, 5 };
        MonolayerVertexMeshGenerator builder(nodes, "NodeMergeWhenLowIndexNodeMustBeAddedToElement");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(5, node_indices_elem_1);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 12u);

        // Merge nodes 4 and 5
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");

        // Test the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 11u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFaces(), 11u);

        // Test the correct nodes are boundary nodes
        for (unsigned i = 0; i < vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        // Test that the moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 1.005, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 1.0, 1e-8);

        // Test the elements own the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 6u);
        unsigned node_indices_element_0[4] = { 0, 1, 4, 5 };
        unsigned node_indices_element_1[4] = { 1, 2, 3, 4 };
        for (unsigned i = 0; i < 4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
        }
    }

    void TestCheckForSwapsAndIdentifySwapType() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[3] = { 2, 3, 5 };
        unsigned node_indices_elem_1[4] = { 4, 1, 2, 5 };
        unsigned node_indices_elem_2[3] = { 0, 1, 4 };
        unsigned node_indices_elem_3[4] = { 4, 5, 3, 0 };

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapWith4Elements");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(3, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        // Set the threshold distance between vertices for a T1 swap as follows
        // so that it will trigger CheckForSwapsFromShortEdges
        vertex_mesh.SetCellRearrangementThreshold(0.3);
        vertex_mesh.CheckForSwapsFromShortEdges();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "AfterOnce");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.275, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.725, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.275, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.725, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1, 1e-8);

        // Test that each element contains the correct nodes following the rearrangement
        const unsigned node_indices_element_0[4] = { 2, 5, 4, 3 };
        const unsigned node_indices_element_1[3] = { 1, 2, 5 };
        const unsigned node_indices_element_2[4] = { 0, 4, 5, 1 };
        const unsigned node_indices_element_3[3] = { 4, 3, 0 };
        for (unsigned i = 0; i < 4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);
        // Perform a T1 swap on nodes 4 and 5
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5));
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "AfterTwice");

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0 + 0.2 * sqrt(41.0) + 2 * 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2 + 0.2 * sqrt(41.0) + 2 * 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0 + 0.2 * sqrt(41.0) + 2 * 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 1.2 + 0.2 * sqrt(41.0) + 2 * 0.3, 1e-6);

        // Test T1 swap location tracking
        std::vector<c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 4u);
        for (unsigned i = 0; i < 4; ++i)
        {
            TS_ASSERT_DELTA(t1_locations[i][0], 0.5, 1e-6);
            TS_ASSERT_DELTA(t1_locations[i][1], 0.5, 1e-6);
            TS_ASSERT_DELTA(t1_locations[i][2], i % 2, 1e-6);
        }

        // Keep testing...
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 17u);
    }

    void TestT1SwapAsynchronous() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap doesn't occur when max(l1,l2)>mCellRearrangementThreshold.
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.48, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.52, 0.0));

        unsigned node_indices_elem_0[3] = { 2, 3, 5 };
        unsigned node_indices_elem_1[4] = { 4, 1, 2, 5 };
        unsigned node_indices_elem_2[3] = { 0, 1, 4 };
        unsigned node_indices_elem_3[4] = { 4, 5, 3, 0 };

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapAsynchronous");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(3, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        // Set the threshold distance between vertices for a T1 swap as follows
        // so that it will not trigger CheckForSwapsFromShortEdges
        vertex_mesh.GetNode(11)->rGetModifiableLocation()[1] = 0.65;
        vertex_mesh.GetNode(10)->rGetModifiableLocation()[1] = 0.35;
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Initial_1");
        vertex_mesh.SetCellRearrangementThreshold(0.05);
        vertex_mesh.SetCellRearrangementRatio(2);
        vertex_mesh.CheckForSwapsFromShortEdges();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After_1");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.55, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.45, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);

        vertex_mesh.GetNode(4)->rGetModifiableLocation()[1] = 0.36;
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Initial_2");
        vertex_mesh.SetCellRearrangementThreshold(0.23);
    }

    void TestT1SwapNonEvenFace() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap doesn't occur when max(l1,l2)>mCellRearrangementThreshold.
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[3] = { 2, 3, 5 };
        unsigned node_indices_elem_1[4] = { 4, 1, 2, 5 };
        unsigned node_indices_elem_2[3] = { 0, 1, 4 };
        unsigned node_indices_elem_3[4] = { 4, 5, 3, 0 };

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapNonEvenFace");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(3, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        vertex_mesh.GetNode(4)->rGetModifiableLocation()[0] = 0.45;
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        // Set the threshold distance between vertices for a T1 swap as follows
        // so that it will trigger CheckForSwapsFromShortEdges
        vertex_mesh.SetCellRearrangementThreshold(0.3);
        vertex_mesh.CheckForSwapsFromShortEdges();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.2518, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5278, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0055, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.6981, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.4721, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], -0.0055, 1e-4);
    }

    void TestPerformT1SwapOnBoundary() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in three elements such that all nodes are
         * boundary nodes, as shown below. We will test that that a T1 swap is correctly implemented.
         *  _____
         * |\   /
         * | \ /
         * |  |
         * | / \
         * |/___\
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.4));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6));

        unsigned node_indices_elem_0[3] = { 2, 3, 5 };
        unsigned node_indices_elem_1[3] = { 1, 4, 0 };
        unsigned node_indices_elem_2[4] = { 0, 4, 5, 3 };

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapOnBoundary");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(3, node_indices_elem_1);
        builder.BuildElementWith(4, node_indices_elem_2);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);
        // Perform a T1 swap on nodes 5 and 4 (this way round to ensure coverage of boundary node tracking)
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-8);

        // Test that each element contains the correct number nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 14u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[4] = { 2, 5, 4, 3 };
        unsigned node_indices_element_1[4] = { 0, 4, 5, 1 };
        unsigned node_indices_element_2[4] = { 0, 4, 3 };
        for (unsigned i = 0; i < 4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2 + 0.2 * sqrt(41.0) + 2 * 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2 + 0.2 * sqrt(41.0) + 2 * 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0 + 0.2 * sqrt(41.0) + 2 * 0.2, 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i = 0; i < vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = !(i == 4 || i == 10);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestPerformT1SwapOnBoundary2() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in three elements such that all but one node
         * are boundary nodes, as shown below. We will test that that a T1 swap is correctly implemented.
         *
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[4] = { 1, 2, 5, 4 };
        unsigned node_indices_elem_1[3] = { 1, 4, 0 };
        unsigned node_indices_elem_2[4] = { 0, 4, 5, 3 };

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapOnBoundary2");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(3, node_indices_elem_1);
        builder.BuildElementWith(4, node_indices_elem_2);

        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 14u);

        // Perform a T1 swap on nodes 5 and 4 (this way round to ensure coverage of boundary node tracking)
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-8);

        // Test that each element contains the correct number nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumFaces(), 5u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[3] = { 1, 2, 5 };
        unsigned node_indices_element_1[4] = { 0, 4, 5, 1 };
        unsigned node_indices_element_2[3] = { 0, 4, 3 };
        for (unsigned i = 0; i < 4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0 + 0.2 * sqrt(41.0) + 2 * 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2 + 0.2 * sqrt(41.0) + 2 * 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0 + 0.2 * sqrt(41.0) + 2 * 0.2, 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i = 0; i < vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    void TestPerformT1SwapWhenVoidForms() throw(Exception)
    {
        /*
         * Create a mesh containing six nodes containing in two elements. We will test that
         * a T1 swap is correctly performed in the case where a void forms as a result of
         * the rearrangement, as shown below.
         *
         * |\   /|     |\      /|
         * | \ / |     | \    / |
         * |  |  |  => | /    \ |
         * | / \ |     |/      \|
         * |/   \|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[4] = { 0, 4, 5, 3 };
        unsigned node_indices_elem_1[4] = { 4, 1, 2, 5 };

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapWhenVoidForms");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);
        // Perform a T1 swap on nodes 5 and 4.
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4));
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-8);

        // Test that each element contains the correct number of nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 10u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[3] = { 0, 4, 3 };
        unsigned node_indices_element_1[3] = { 1, 2, 5 };
        for (unsigned i = 0; i < 3; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0 + 0.2 * sqrt(41.0) + 2 * 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0 + 0.2 * sqrt(41.0) + 2 * 0.2, 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i = 0; i < vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    void TestPerformT1SwapExceptions() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes containing in two triangle and two rhomboid elements,
         * where two nodes (those with indices 4 and 5) have the same location. We will test that
         * trying to perform a T1 swap on these nodes throws the correct exception.
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.5, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.5, 0.0));

        unsigned node_indices_elem_0[3] = { 2, 3, 5 };
        unsigned node_indices_elem_1[4] = { 2, 5, 4, 1 };
        unsigned node_indices_elem_2[3] = { 1, 4, 0 };
        unsigned node_indices_elem_3[4] = { 0, 4, 5, 3 };

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapExceptions");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(3, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

        // Test that trying to perform a T1 swap on nodes 4 and 5 throws the correct exception
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5)), "Nodes(5&4) are too close together, this shouldn't happen");
    }

    void TestDoNotPerforT1SwapWithRemovingEdgeFromTriangularElement() throw(Exception)
    {
        /**
         * In this test we check that a T1 swap does not occur if one of the elements is triangular
         * and would loose an edge by swapping nodes. The mesh looks like this
         *
         *       ______________
         *      |\             |
         *      | \ _________  |
         *      |  |          \| ...where the funny shaped element in the middle is supposed to be
         *      |  |_________ /|    a very long triangle that has the third vertex on the right hand boundary.
         *      | /            |
         *      |/_____________|
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 2.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 2.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.3, 0.95, 0.0));
        nodes.push_back(new Node<3>(5, true, 2.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(6, false, 0.3, 1.05, 0.0));

        unsigned node_indices_elem_0[4] = { 0, 1, 5, 4 };
        unsigned node_indices_elem_1[4] = { 5, 2, 3, 6 };
        unsigned node_indices_elem_2[4] = { 0, 4, 6, 3 };
        unsigned node_indices_elem_3[3] = { 4, 5, 6 };

        MonolayerVertexMeshGenerator builder(nodes, "T1NoSwapWithTriangularPrism");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(4, node_indices_elem_2);
        builder.BuildElementWith(3, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        // Ensure that the inner edge will be considered for a swap
        vertex_mesh.SetCellRearrangementThreshold(0.11);

        // Check for T1 swaps and carry them out if allowed - the short edge should not swap!
        vertex_mesh.CheckForSwapsFromShortEdges();

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 6u);

        // Test that each element still contains the correct nodes following the rearrangement
        for (unsigned i = 0; i < 4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_elem_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_elem_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_elem_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_elem_3[i]);
            }
        }
    }

    void TestExceptionForVoidRemovalWithRemovingEdgeFromTriangularElement() throw(Exception)
    {
        /**
         * In this test we check that void removal does not occur if one of the adjacent elements is triangular
         * and would loose an edge by swapping nodes. The code should throw and exception in this case.
         * The mesh looks like this
         *
         *       ______________./This corner is not a node.
         *      |\      1      |
         *      | \ _________  |
         *      |  |   void   \| ...where elements 1, and 2 are triangles that share the right hand vertex
         *      |  |_________ /|    with the triangular void in the middle.
         *      | /     2      |
         *      |/_____________|.This corner is not a node either.
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.3, 0.95, 0.0));
        nodes.push_back(new Node<3>(3, true, 2.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.3, 1.05, 0.0));

        unsigned node_indices_elem_0[4] = { 0, 2, 4, 1 };
        unsigned node_indices_elem_1[3] = { 1, 4, 3 };
        unsigned node_indices_elem_2[3] = { 0, 3, 2 };

        MonolayerVertexMeshGenerator builder(nodes, "NoT1SwapWithTriangularVoid");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(3, node_indices_elem_1);
        builder.BuildElementWith(3, node_indices_elem_2);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        // Ensure that the inner edge will be considered for a swap
        vertex_mesh.SetCellRearrangementThreshold(0.11);

        // Check for possible swaps and carry them out if allowed - the short edge should not swap and
        // the void should not be removed!
        TS_ASSERT_THROWS_THIS(vertex_mesh.CheckForSwapsFromShortEdges(),
                              "Triangular element next to triangular void, not implemented yet.");
    }

    void TestReMeshForT1Swaps() throw(Exception)
    {
        /*
         * Read in a vertex mesh that contains several pairs of nodes that are close enough for
         * T1 swaps to be performed, as shown below. The mesh consists of six elements and all
         * T1 swaps are performed on all horizontal edges. We will test that the ReMesh() method
         * correctly performs T1 swaps for internal and boundary elements, and correctly updates
         * which nodes are labelled as boundary nodes.
         *
         *      /\    /\
         *     /  \__/  \
         *    /   /  \   \
         *    \__/\__/\__/
         *    /  \/  \/  \
         *    \   \__/   /
         *     \  /  \  /
         *      \/    \/
         */
        VertexMeshReader<2, 2> mesh_reader("cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1");
        MutableVertexMesh<2, 2> vertex_2mesh;
        vertex_2mesh.ConstructFromMeshReader(mesh_reader);

        MonolayerVertexMeshGenerator builder("T1SwapReMesh");
        MutableVertexMesh<3, 3>& vertex_mesh = *(builder.MakeMeshUsing2dMesh(vertex_2mesh));
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 44u);

        // Calls ReMesh() to identify and perform any T1 swaps
        vertex_mesh.SetCellRearrangementThreshold(0.1);
        vertex_mesh.ReMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 44u);

        std::string dirname = OUTPUT_NAME + std::string("/T1SwapReMesh");
        std::string mesh_filename = "vertex33_remesh_T1_after_remesh";

        // Save the mesh data using mesh writers
        VertexMeshWriter<3, 3> mesh_writer(dirname, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh(vertex_mesh);

        // Check the positions are updated correctly
        OutputFileHandler handler(dirname, false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex33_remesh_T1_after_remesh.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex33_remesh_T1_after_remesh.cell";

        FileComparison comparer1(results_file1, "cell_based/test/data/TestMutableVertexMesh/vertex33_remesh_T1_after_remesh.node");
        TS_ASSERT(comparer1.CompareFiles());
        FileComparison comparer2(results_file2, "cell_based/test/data/TestMutableVertexMesh/vertex33_remesh_T1_after_remesh.cell");
        TS_ASSERT(comparer2.CompareFiles());
    }

    void TestPerformT2Swap() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in three trapezium element and
         * a central triangle element, as shown below. We will test that a T2 swap
         * correctly removes the triangle element from the mesh.
         *
         *      /|\
         *     / | \
         *    /  |  \    (the triangular element has index zero)
         *   /2 /_\ 1\
         *  /  /   \  \
         * /__/__3__\__\
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.5, 0.5, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.4, 0.2, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.6, 0.2, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.3, 0.0));

        unsigned node_indices_elem_0[3] = { 3, 4, 5 };
        unsigned node_indices_elem_1[4] = { 1, 2, 5, 4 };
        unsigned node_indices_elem_2[4] = { 2, 0, 3, 5 };
        unsigned node_indices_elem_3[4] = { 0, 1, 4, 3 };

        MonolayerVertexMeshGenerator builder(nodes, "T2Swap");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(4, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        // Perform a T2 swap on the central triangle element
        VertexElement<3, 3>* p_element_0 = vertex_mesh.GetElement(0);
        c_vector<double, 3> centroid_of_element_0_before_swap = vertex_mesh.GetCentroidOfElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 12u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFaces(), 17u);

        for (unsigned j = 1; j < 4; j++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumNodes(), 6u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumFaces(), 5u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNodeGlobalIndex(0), j % 3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNodeGlobalIndex(1), (j + 1) % 3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNodeGlobalIndex(2), 12u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumFaces(), 5u);
        }

        // Test boundary property of nodes. All are boundary nodes except node 3.
        for (unsigned i = 0; i < vertex_mesh.GetNumAllNodes(); i++)
        {
            if (vertex_mesh.GetNode(i)->IsDeleted())
            {
                continue;
            }
            bool expected_boundary_node = !(i == 12 || i == 13);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

        // Test the location of the new node:
        TS_ASSERT_DELTA(vertex_mesh.GetNode(12)->rGetLocation()[0], centroid_of_element_0_before_swap[0], 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(12)->rGetLocation()[1], centroid_of_element_0_before_swap[1], 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(12)->rGetLocation()[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(13)->rGetLocation()[0], centroid_of_element_0_before_swap[0], 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(13)->rGetLocation()[1], centroid_of_element_0_before_swap[1], 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(13)->rGetLocation()[2], 1.0, 1e-10);
        // Test the tracking of the T2 swap location:
        TS_ASSERT_DELTA(vertex_mesh.GetLastT2SwapLocation()[0], centroid_of_element_0_before_swap[0], 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetLastT2SwapLocation()[1], centroid_of_element_0_before_swap[1], 1e-10);

        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.RemoveDeletedNodesAndElements(map);
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "AfterRemove");
    }

    void TestPerformT2SwapOnBoundary() throw(Exception)
    {
        /*
         * Create a mesh comprising six nodes contained in two trapezium elements
         * and one triangle element, as shown below. We will test that a T2 swap
         * is performed correctly when boundary nodes are involved.
         *
         *       /|\
         *      / | \
         *     /  |  \    (the triangular element has index zero)
         *    /2 /_\ 1\
         *   /  /   \  \
         *  /__/     \__\
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.5, 0.5));
        nodes.push_back(new Node<3>(3, true, 0.4, 0.25));
        nodes.push_back(new Node<3>(4, true, 0.6, 0.25));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.3));

        unsigned node_indices_elem_0[3] = { 3, 4, 5 };
        unsigned node_indices_elem_1[4] = { 1, 2, 5, 4 };
        unsigned node_indices_elem_2[4] = { 2, 0, 3, 5 };

        MonolayerVertexMeshGenerator builder(nodes, "T2SwapOnBoundary");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(4, node_indices_elem_2);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 14u);

        // Perform a T2 swap on the central triangle element
        VertexElement<3, 3>* p_element_0 = vertex_mesh.GetElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 9u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFaces(), 14u);

        for (unsigned j = 1; j < 3; j++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumNodes(), 6u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumFaces(), 5u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNodeGlobalIndex(0), j % 3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNodeGlobalIndex(1), (j + 1) % 3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNodeGlobalIndex(2), 12u);
        }

        // Test boundary property of nodes. All are boundary nodes.
        for (unsigned i = 0; i < vertex_mesh.GetNumAllNodes(); i++)
        {
            if (vertex_mesh.GetNode(i)->IsDeleted())
                continue;
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.RemoveDeletedNodesAndElements(map);
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "AfterRemove");
    }

    void TestPerformT2OnBoundary2() throw(Exception)
    {
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
        // Make five nodes to assign to two elements
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 0.5, 0.5));
        nodes.push_back(new Node<3>(2, true, 0.4, 0.25));
        nodes.push_back(new Node<3>(3, true, 0.6, 0.25));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.3));

        const unsigned node_indices_elem_0[3] = { 2, 3, 4 };
        const unsigned node_indices_elem_1[4] = { 0, 1, 4, 3 };

        MonolayerVertexMeshGenerator builder(nodes, "T2SwapOnBoundary2");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 10u);

        // Perform a T2 swap on the central triangle element
        VertexElement<3, 3>* p_element_0 = vertex_mesh.GetElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFaces(), 10u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(4), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(5), 11u);

        // Test boundary property of nodes. All are boundary nodes.
        for (unsigned i = 0; i < vertex_mesh.GetNumAllNodes(); i++)
        {
            if (vertex_mesh.GetNode(i)->IsDeleted())
                continue;
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.RemoveDeletedNodesAndElements(map);
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "AfterRemove");
    }

    void TestT2SwapsDontOccurWithTriangularNeighbours() throw(Exception)
    {
        // Make 6 nodes to assign to four elements
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.5, 0.5));
        nodes.push_back(new Node<3>(3, false, 0.4, 0.25));
        nodes.push_back(new Node<3>(4, false, 0.6, 0.25));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.3));

        // Make two triangles and two trapezium elements out of these nodes
        const unsigned node_indices_elem_0[3] = { 3, 4, 5 };
        const unsigned node_indices_elem_1[3] = { 2, 5, 4 };
        const unsigned node_indices_elem_2[4] = { 2, 0, 3, 5 };
        const unsigned node_indices_elem_3[4] = { 0, 1, 4, 3 };

        MonolayerVertexMeshGenerator builder(nodes, "T2NoSwapWithTriangularNeighbours");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(3, node_indices_elem_1);
        builder.BuildElementWith(4, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);

        // Attempt to perform a T2 swap on the middle triangle element
        VertexElement<3, 3>* p_element_0 = vertex_mesh.GetElement(0);
        TS_ASSERT_THROWS_THIS(vertex_mesh.PerformT2Swap(*p_element_0),
                              "One of the neighbours of a small triangular element is also a triangle - "
                              "dealing with this has not been implemented yet");
    }

    void TestPerformT2SwapWithRosettes() throw(Exception)
    {
        /* Create a mesh containing a smaller triangular element, each of whose nodes are
         * 'rosette' nodes. Test that a T2 swap correctly removes the triangular element
         * from the mesh.
         *  _________
         *  |\      |
         *  | \     |
         *  | |_\___|
         *  | / |   |
         *  |/__|___|
         */

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.99, 1.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 2.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.99, 1.01));
        nodes.push_back(new Node<3>(7, false, 0.0, 2.0));
        nodes.push_back(new Node<3>(8, false, 2.0, 2.0));

        // Make two triangles and two trapezium elements out of these nodes
        const unsigned node_indices_elem_0[4] = { 0, 1, 4, 3 };
        const unsigned node_indices_elem_1[4] = { 1, 2, 5, 4 };
        const unsigned node_indices_elem_2[4] = { 0, 3, 6, 7 };
        const unsigned node_indices_elem_3[3] = { 3, 4, 6 };
        const unsigned node_indices_elem_4[5] = { 4, 5, 8, 7, 6 };

        MonolayerVertexMeshGenerator builder(nodes, "T2SwapWithRosette");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(4, node_indices_elem_2);
        builder.BuildElementWith(3, node_indices_elem_3);
        builder.BuildElementWith(5, node_indices_elem_4);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 23u);

        // Perform a T2 swap on the central triangle element
        VertexElement<3, 3>* p_element_3 = vertex_mesh.GetElement(3);
        c_vector<double, 3> centroid_of_element_3_before_swap = vertex_mesh.GetCentroidOfElement(3);
        vertex_mesh.PerformT2Swap(*p_element_3);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 18u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 20u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 18u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), 18u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 7u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumFaces(), 0u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(0), 18u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(2), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(3), 7u);

        VertexElementMap map(vertex_mesh.GetNumElements());
        vertex_mesh.RemoveDeletedNodesAndElements(map);
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "AfterRemove");
    }

    void TestDivideElement() throw(Exception)
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, true, 1.0, 2.0));
        nodes.push_back(new Node<3>(5, true, 0.0, 2.0));

        const unsigned node_indices_elem_0[4] = { 0, 1, 2, 3 };
        const unsigned node_indices_elem_1[4] = { 2, 4, 5, 3 };

        MonolayerVertexMeshGenerator builder(nodes, "DivideElement");
        builder.BuildElementWith(4, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>* p_mesh = builder.GenerateMesh();
        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

        c_vector<double, 3> axis_of_division = zero_vector<double>(3);
        axis_of_division[0] = 1;
        p_mesh->DivideElementAlongGivenAxis(p_mesh->GetElement(0), axis_of_division);
        p_mesh->ReMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumNodes(), 8u);

        builder.WriteVtkWithSubfolder(OUTPUT_NAME, "After");
    }

    void TestDivideElement2() throw(Exception)
    {
        HexagonalPrism3dVertexMeshGenerator generator(1, 1);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        {
            VertexMeshWriter<3, 3> writer(OUTPUT_NAME + "/DivideElement2", "DivideElement2");
            writer.WriteVtkUsingMesh(*p_mesh, "Initial");
        }

        c_vector<double, 3> axis_of_division = zero_vector<double>(3);
        axis_of_division[0] = 1;
        axis_of_division[1] = 0;
        axis_of_division[2] = 0.2;
        axis_of_division /= norm_2(axis_of_division);
        p_mesh->DivideElementAlongGivenAxis(p_mesh->GetElement(0), axis_of_division);
        p_mesh->ReMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 13u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNumNodes(), 8u);

        {
            VertexMeshWriter<3, 3> writer(OUTPUT_NAME + "/DivideElement2", "DivideElement2", false);
            writer.WriteVtkUsingMesh(*p_mesh, "After");
        }
    }

    // // Commented this test as T2 Swap should only happen to triangular prism
    // // element in vertex model
    // void TestPerformT2SwapWithRosettes2() throw(Exception)
    // {
    //     /* Create a mesh containing a smaller square element.
    //         * Test that a T2 swap correctly removes the triangular element
    //         * from the mesh and create rosette out of it.
    //         *  ___________
    //         *  |    |    |
    //         *  |    |    |
    //         *  |___/ \___|
    //         *  |   \ /   |
    //         *  |    |    |
    //         *  |____|____|
    //         */

    //     std::vector<Node<3>*> nodes;
    //     nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
    //     nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
    //     nodes.push_back(new Node<3>(2, true, 2.0, 0.0));
    //     nodes.push_back(new Node<3>(3, true, 2.0, 1.0));
    //     nodes.push_back(new Node<3>(4, true, 2.0, 2.0));
    //     nodes.push_back(new Node<3>(5, true, 1.0, 2.0));
    //     nodes.push_back(new Node<3>(6, true, 0.0, 2.0));
    //     nodes.push_back(new Node<3>(7, true, 0.0, 1.0));
    //     nodes.push_back(new Node<3>(8, false, 1.0, 0.8));
    //     nodes.push_back(new Node<3>(9, false, 1.2, 1.0));
    //     nodes.push_back(new Node<3>(10, false, 1.0, 1.2));
    //     nodes.push_back(new Node<3>(11, false, 0.8, 1.0));

    //     const unsigned node_indices_elem_0[5] = { 0, 1, 8, 11, 7 };
    //     const unsigned node_indices_elem_1[5] = { 2, 3, 9, 8, 1 };
    //     const unsigned node_indices_elem_2[5] = { 4, 5, 10, 9, 3 };
    //     const unsigned node_indices_elem_3[5] = { 6, 7, 11, 10, 5 };
    //     const unsigned node_indices_elem_4[4] = { 8, 9, 10, 11 };

    //     Helper3dVertexMeshBuilder builder(nodes, "T2SwapWithRosette2");
    //     builder.BuildElementWith(5, node_indices_elem_0);
    //     builder.BuildElementWith(5, node_indices_elem_1);
    //     builder.BuildElementWith(5, node_indices_elem_2);
    //     builder.BuildElementWith(5, node_indices_elem_3);
    //     builder.BuildElementWith(4, node_indices_elem_4);

    //     // A reference variable as mesh is noncopyable
    //     MutableVertexMesh<3, 3>& vertex_mesh = *builder.GenerateMesh();
    //     builder.WriteVtkWithSubfolder(OUTPUT_NAME, "Before");

    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 5u);
    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 24u);
    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 26u);

    //     // Perform a T2 swap on the central triangle element
    //     VertexElement<3, 3>* p_element_4 = vertex_mesh.GetElement(4);
    //     c_vector<double, 3> centroid_of_element_4_before_swap = vertex_mesh.GetCentroidOfElement(3);
    //     vertex_mesh.PerformT2Swap(*p_element_4);

    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 18u);
    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 20u);

    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 5u);
    //     TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 26u);

    //     for (unsigned i = 0; i < 4; ++i)
    //     {
    //         TS_ASSERT_EQUALS(vertex_mesh.GetElement(i)->GetNumNodes(), 8u);
    //         TS_ASSERT_EQUALS(vertex_mesh.GetElement(i)->GetNumFaces(), 6u);
    //         TS_ASSERT_EQUALS(vertex_mesh.GetElement(i)->GetNodeGlobalIndex(0), i * 2);
    //         TS_ASSERT_EQUALS(vertex_mesh.GetElement(i)->GetNodeGlobalIndex(1), i * 2 + 1);
    //         TS_ASSERT_EQUALS(vertex_mesh.GetElement(i)->GetNodeGlobalIndex(2), 24u);
    //         TS_ASSERT_EQUALS(vertex_mesh.GetElement(i)->GetNodeGlobalIndex(3), (i * 2 - 1 + 8) % 8);
    //     }

    //     TS_ASSERT_EQUALS(vertex_mesh.GetNode(24)->IsBoundaryNode(), false);
    //     TS_ASSERT_EQUALS(vertex_mesh.GetNode(25)->IsBoundaryNode(), false);

    //     VertexElementMap map(vertex_mesh.GetNumElements());
    //     vertex_mesh.RemoveDeletedNodesAndElements(map);
    //     builder.WriteVtkWithSubfolder(OUTPUT_NAME, "AfterRemove");
    // }
};

#endif /*TESTMUTABLEVERTEXMESH33REMESH_HPP_*/
