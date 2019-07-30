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

#ifndef TESTCYLINDRICAL2DMESH_HPP_
#define TESTCYLINDRICAL2DMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

#include "UblasCustomFunctions.hpp"
#include "VertexMesh.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "ArchiveOpener.hpp"
#include "TrianglesMeshWriter.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCylindrical2dMesh : public CxxTest::TestSuite
{
public:

    void TestCreateMirrorCellsAndAlignmentTester()
    {
        // Note that elements are not created (and boundary elements are not changed),
        // this just creates a set of new nodes

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Reset the mesh
        p_mesh = generator.GetCylindricalMesh();
        p_mesh->CreateMirrorNodes();

        // Check the vectors are the right size...
        TS_ASSERT_EQUALS(p_mesh->mLeftOriginals.size(), 36u);
        TS_ASSERT_EQUALS(p_mesh->mRightOriginals.size(), 36u);
        TS_ASSERT_EQUALS(p_mesh->mLeftOriginals.size(), p_mesh->mLeftImages.size());
        TS_ASSERT_EQUALS(p_mesh->mRightOriginals.size(), p_mesh->mRightImages.size());

        // Check that the image nodes are where they should be.
        for (unsigned i=0; i<p_mesh->mLeftOriginals.size(); i++)
        {
            c_vector<double,2> original_location;
            original_location = p_mesh->GetNode(p_mesh->mLeftOriginals[i])->rGetLocation();
            c_vector<double,2> image_location;
            image_location = p_mesh->GetNode(p_mesh->mLeftImages[i])->rGetLocation();

            TS_ASSERT_DELTA(original_location[0] + crypt_width, image_location[0], 1e-7);
            TS_ASSERT_DELTA(original_location[1], image_location[1], 1e-7);
        }

        for (unsigned i=0; i<p_mesh->mRightOriginals.size(); i++)
        {
            c_vector<double,2> original_location;
            original_location = p_mesh->GetNode(p_mesh->mRightOriginals[i])->rGetLocation();
            c_vector<double,2> image_location;
            image_location = p_mesh->GetNode(p_mesh->mRightImages[i])->rGetLocation();

            TS_ASSERT_DELTA(original_location[0] - crypt_width, image_location[0], 1e-7);
            TS_ASSERT_DELTA(original_location[1], image_location[1], 1e-7);
        }

        // Check that we've got the correct number
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up*2);

        // Cheat to put something into mTopHaloNodes to make exception throw
        p_mesh->mTopHaloNodes.push_back(234u);
        unsigned corresponding_node_index = p_mesh->GetCorrespondingNodeIndex(0u);

        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[0]+crypt_width, p_mesh->GetNode(corresponding_node_index)->rGetLocation()[0], 1e-9);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[1], p_mesh->GetNode(corresponding_node_index)->rGetLocation()[1], 1e-9);
    }

    void TestReconstructCylindricalMesh()
    {
        // This test takes in a new mesh created using the mirror function above
        // and a ReMesh call, then removes nodes, elements and boundary elements

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;

        // Set up a mesh which can be mirrored (no ghost nodes in this case)
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Create a mirrored load of nodes for the normal remesher to work with
        p_mesh->CreateMirrorNodes();

        // Call the normal re-mesh
        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->MutableMesh<2,2>::ReMesh(map);

        // Re-Index the vectors regarding left/right nodes with the node map
        for (unsigned i = 0; i<p_mesh->mLeftOriginals.size(); i++)
        {
             p_mesh->mLeftOriginals[i] = map.GetNewIndex(p_mesh->mLeftOriginals[i]);
             p_mesh->mLeftImages[i] = map.GetNewIndex(p_mesh->mLeftImages[i]);
        }
        for (unsigned i = 0; i<p_mesh->mRightOriginals.size(); i++)
        {
             p_mesh->mRightOriginals[i] = map.GetNewIndex(p_mesh->mRightOriginals[i]);
             p_mesh->mRightImages[i] = map.GetNewIndex(p_mesh->mRightImages[i]);
        }

        p_mesh->ReconstructCylindricalMesh();

        unsigned elements_for_node_0 = 0;
        unsigned elements_for_node_11 = 0;
        unsigned elements_for_node_12 = 0;
        unsigned elements_for_node_18 = 0;

        unsigned checksum_for_node_0 = 0;
        unsigned checksum_for_node_11 = 0;
        unsigned checksum_for_node_12 = 0;
        unsigned checksum_for_node_18 = 0;

        for (Cylindrical2dMesh::ElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
             elem_iter != p_mesh->GetElementIteratorEnd();
             ++elem_iter)
        {
            for (unsigned i=0; i<3; i++)
            {
                unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);

                if (this_node_index==0)
                {
                    elements_for_node_0++;
                    for (unsigned j=0; j<3; j++)
                    {
                        checksum_for_node_0 += elem_iter->GetNodeGlobalIndex(j);
                    }
                }
                if (this_node_index==11)
                {
                    elements_for_node_11++;
                    for (unsigned j=0; j<3; j++)
                    {
                        checksum_for_node_11 += elem_iter->GetNodeGlobalIndex(j);
                    }
                }
                if (this_node_index==12)
                {
                    elements_for_node_12++;
                    for (unsigned j=0; j<3; j++)
                    {
                        checksum_for_node_12 += elem_iter->GetNodeGlobalIndex(j);
                    }
                }
                if (this_node_index==18)
                {
                    elements_for_node_18++;
                    for (unsigned j=0; j<3; j++)
                    {
                        checksum_for_node_18 += elem_iter->GetNodeGlobalIndex(j);
                    }
                }
            }
        }

        TS_ASSERT_EQUALS(elements_for_node_0, 3u);
        TS_ASSERT_EQUALS(elements_for_node_11, 6u);
        TS_ASSERT_EQUALS(elements_for_node_12, 6u);
        TS_ASSERT_EQUALS(elements_for_node_18, 6u);

        // These are nodes on the edge of the cylindrical region.
        // If the mesh is periodic they will be joined to the following other nodes.
        unsigned checksum_target_for_node_0 = (0 + 5 + 11) + (0 + 6 + 11) + (0 + 1 + 6);
        TS_ASSERT_EQUALS(checksum_for_node_0, checksum_target_for_node_0);

        unsigned checksum_target_for_node_11 = (0 + 5 + 11) + (0 + 6 + 11) + (11 + 12 + 6)
                                              +(11 + 17 + 12) + (10 + 11 + 17) + (5 + 10 + 11);
        TS_ASSERT_EQUALS(checksum_for_node_11, checksum_target_for_node_11);

        unsigned checksum_target_for_node_12 = (12 + 13 + 18) + (12 + 18 + 23) + (12 + 17 + 23)
                                              +(12 + 11 + 17) + (12 + 6 + 11) + (12 + 6 + 13);
        TS_ASSERT_EQUALS(checksum_for_node_12, checksum_target_for_node_12);

        unsigned checksum_target_for_node_18 = (18 + 19 + 25) + (18 + 24 + 25) + (18 + 23 + 24)
                                                +(18 + 12 + 23) + (18 + 12 + 13) + (18 + 13 + 19);
        TS_ASSERT_EQUALS(checksum_for_node_18, checksum_target_for_node_18);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 12u);
    }

    void TestCylindricalReMesh()
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;

        // Set up a mesh which can be mirrored (no ghosts in this case)
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u);  // No boundary elements now the halo nodes are removed
    }

    /*
     * Failing test for ReMesh (see #1275)
     */
    void noTestCylindricalReMeshFailingTest()
    {
        // Load a problematic mesh
        TrianglesMeshReader<2,2> mesh_reader("cell_based/test/data/TestCylindricalMeshBug/mesh");
        Cylindrical2dMesh mesh(20);
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA(mesh.GetWidth(0), 20, 1e-3);

        TrianglesMeshWriter<2,2> mesh_writer("TestCylindricalMeshBug", "mesh", false);
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Use showme to view this mesh

        NodeMap map(mesh.GetNumNodes());
        mesh.ReMesh(map);
    }

    void TestCylindricalReMeshAfterDelete()
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;

        // Set up a mesh which can be mirrored (no ghosts in this case)
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        unsigned num_old_nodes = p_mesh->GetNumNodes();

        p_mesh->DeleteNode(15);

        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up-1);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u);  // No boundary elements now the halo nodes are removed

        TS_ASSERT_EQUALS(map.GetSize(), num_old_nodes);
        TS_ASSERT_EQUALS(map.IsDeleted(15), true);

        for (unsigned i=0; i<num_old_nodes; i++)
        {
            if (i<15)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), i);
            }
            if (i>15)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-1));
            }
        }
   }

    void TestCylindricalReMeshOnSmallMesh()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;

        // Set up a mesh which can be mirrored (no ghosts in this case)
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u); // boundary elements removed now halo nodes are used
    }

    void TestGetVectorBetweenCyclindricalPoints()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;

        // Set up a mesh which can be mirrored (no ghosts in this case)
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        c_vector<double, 2> location1 = p_mesh->GetNode(1)->rGetLocation();
        c_vector<double, 2> location2 = p_mesh->GetNode(4)->rGetLocation();

        // Test a normal distance calculation...
        c_vector<double, 2> vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], sqrt(3.0)/2.0, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetDistanceBetweenNodes(1, 4), 1.0, 1e-7);

        // Cylindrical at x=3, so closest node should be node 0 at (0,0).
        c_vector<double, 2> test;
        test[0] = 2.9;
        test[1] = 0.0;
        TS_ASSERT_EQUALS(p_mesh->GetNearestNodeIndex(test), 0u);

        // ...and the opposite vector
        vector = p_mesh->GetVectorFromAtoB(location2, location1);
        TS_ASSERT_DELTA(vector[0], -0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], -sqrt(3.0)/2.0, 1e-4);

        // Test a periodic calculation
        location1[0] = 0.5;
        location1[1] = 3.0;
        location2[0] = 2.5;
        location2[1] = 4.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], -1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], 1.0, 1e-7);

        // Test a periodic calculation where points need to be swapped
        location1[0] = 2.5;
        location1[1] = 4.0;
        location2[0] = 0.5;
        location2[1] = 3.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], -1.0, 1e-7);

        // We also want GetVectorFromAtoB to work when the x coord of A and B is not
        // between 0 and crypt width, by first normalizing the x coord (by taking modulus by crypt width)
        location1[0] = -0.5;
        location1[1] = 0;
        location2[0] = 2.5;
        location2[1] = 1;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], 1.0, 1e-7);

        // Some tests where the location[0] is not between -0.25*crypt_width and 1.25*crypt_width
        location1[0] = -2.5;
        location1[1] = 0;
        location2[0] = 4.5;
        location2[1] = 1;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], 1.0, 1e-7);
    }

    void TestSetNodeLocationForCylindricalMesh()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        double crypt_width = 3.0;
        unsigned thickness_of_ghost_layer = 2;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        /*
         * New test to check that the top and bottom rows move together
         * bottom row = 0,1,2
         * top row = 18,19,20
         */
//        c_vector<double, 2> new_location = p_mesh->GetNode(0)->rGetLocation();
//        new_location[1] = -1.760;
//        ChastePoint<2> boundary_point(new_location);
//        // We just move one of the bottom boundary nodes and then...
//        p_mesh->SetNode(0, boundary_point, false);
//        // check that all the nodes on this boundary have moved down
//        for (unsigned i=0; i<3; i++)
//        {
//            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],-1.76000,1e-6);
//        }
//
//        // Same for one of the top boundary nodes
//        new_location = p_mesh->GetNode(19)->rGetLocation();
//        new_location[1] = 4.0;
//        ChastePoint<2> boundary_point2(new_location);
//        p_mesh->SetNode(19, boundary_point2, false);
//        // check that all the nodes on this boundary have moved up
//        for (unsigned i=18; i<21; i++)
//        {
//            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[1],4.0,1e-6);
//        }

        // Move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = 2.999999999;
        new_point_location[1] = -0.866025;
        ChastePoint<2> new_point(new_point_location);
        p_mesh->GetNode(5)->SetPoint(new_point);

        new_point.SetCoordinate(0, -0.0001);

        // This node was on left and is now near the right
        p_mesh->SetNode(0, new_point, false);
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 2.9999, 1e-4);

        new_point.SetCoordinate(0, 1.0000);
        p_mesh->SetNode(0, new_point, false);

        // This node has stayed close to where it was
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 1.0000, 1e-4);

        new_point.SetCoordinate(0, 3.0001);
        p_mesh->SetNode(0, new_point, false);

        // This node was on right and is now on the left
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 0.0001, 1e-4);
    }

    void TestAddNodeAndReMesh()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;

        // Set up a mesh which can be mirrored (no ghosts in this case)
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*(cells_up-1));
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u); // boundary elements removed now halo nodes are used

        c_vector<double,2> point;
        point[0] = -0.05;
        point[1] = 1.0;
        Node<2>* p_node = new Node<2>(p_mesh->GetNumNodes(), point);

        unsigned new_index = p_mesh->AddNode(p_node);
        NodeMap map(p_mesh->GetNumNodes());

        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up+1);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*(cells_up-1)+2);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u); // boundary elements removed now halo nodes are used

        // Check that we have moved the new node across
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 3.0+point[0], 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], point[1], 1e-7);

        // Test GetWidth
        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), 3.0, 1e-9);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), sqrt(3.0), 1e-6);
    }

    void TestHaloNodeInsertionAndRemoval()
    {
        unsigned cells_across = 5;
        unsigned cells_up = 3;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        c_vector<double,2>& rLocation1 = p_mesh->GetNode(1)->rGetModifiableLocation();
        rLocation1[1] -= 0.5;

        c_vector<double,2>& rLocation2 = p_mesh->GetNode(3)->rGetModifiableLocation();
        rLocation2[1] -= 0.4;

        c_vector<double,2>& rLocation3 = p_mesh->GetNode(12)->rGetModifiableLocation();
        rLocation3[1] += 0.8;

        double original_mesh_height = p_mesh->GetWidth(1);
        unsigned original_num_nodes = p_mesh->GetNumNodes();

        p_mesh->CreateHaloNodes();

        unsigned num_original_halo_nodes = p_mesh->GetNumNodes() - original_num_nodes;

        p_mesh->CreateMirrorNodes();

        double new_mesh_height = p_mesh->GetWidth(1);
        unsigned new_num_nodes = p_mesh->GetNumNodes();

        // Halo of nodes is added 0.5 above and below the original mesh.
        TS_ASSERT_DELTA(original_mesh_height, new_mesh_height, 1.0 + 1e-5);
        TS_ASSERT_EQUALS(new_num_nodes, original_num_nodes*2+2*num_original_halo_nodes);

        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->MutableMesh<2,2>::ReMesh(map);   // recreates the boundary elements

        // Test halo node removal
        p_mesh->GenerateVectorsOfElementsStraddlingPeriodicBoundaries();
        p_mesh->CorrectNonPeriodicMesh();
        p_mesh->ReconstructCylindricalMesh();
        p_mesh->DeleteHaloNodes();

        TS_ASSERT_DELTA(original_mesh_height, p_mesh->GetWidth(1), 1.1e-6);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), original_num_nodes);
    }

    void TestHaloNodeReMesh()
    {
        // This test checks that a Halo node remesh can handle a mesh of uneven height

        unsigned cells_across = 5;
        unsigned cells_up = 3;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        c_vector<double,2>& rLocation1 = p_mesh->GetNode(1)->rGetModifiableLocation();
        rLocation1[1] -= 0.5;

        c_vector<double,2>& rLocation2 = p_mesh->GetNode(3)->rGetModifiableLocation();
        rLocation2[1] -= 0.4;

        c_vector<double,2>& rLocation3 = p_mesh->GetNode(12)->rGetModifiableLocation();
        rLocation3[1] += 0.8;

        unsigned total_elements = p_mesh->GetNumElements();
        unsigned total_nodes = p_mesh->GetNumNodes();

        // Check that the ReIndex is working still
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh->GetNumElements());

        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        // Check that we haven't added any nodes or elements by doing this Halo Node ReMesh.
        TS_ASSERT_EQUALS(total_elements, p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(total_nodes, p_mesh->GetNumNodes());

        // Check the ReIndex is working
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), p_mesh->GetNumNodes());
    }

    // NB This checks that periodicity is maintained through archiving...
    void TestArchiving()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "cylindrical_mesh_base.arch";
        ArchiveLocationInfo::SetMeshFilename("cylindrical_mesh");

        // Set up a mesh
        unsigned cells_across = 5;
        unsigned cells_up = 3;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        AbstractTetrahedralMesh<2,2>* const p_mesh = generator.GetCylindricalMesh();

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects
         * changing during the save, and so object tracking leading to wrong results.
         * (e.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.)
         */

        {
            // Serialize the mesh
            double width = p_mesh->GetWidth(0);
            TS_ASSERT_DELTA(width, crypt_width, 1e-7);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost.
            (*p_arch) << p_mesh;
        }

        {
            // De-serialize and compare
            AbstractTetrahedralMesh<2,2>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            // Cylindrical2dMesh now remeshes itself on load (convert from TrianglesMeshReader to normal format)

            TS_ASSERT_DELTA(p_mesh2->GetWidth(0), crypt_width, 1e-7);

            // Compare the loaded mesh against the original
            TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), p_mesh2->GetNumAllNodes());
            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), p_mesh2->GetNumNodes());
            TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryNodes(), p_mesh2->GetNumBoundaryNodes());

            for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
            {
                Node<2>* p_node = p_mesh->GetNode(i);
                Node<2>* p_node2 = p_mesh2->GetNode(i);
                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-16);
                }
            }

            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), p_mesh2->GetNumElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh2->GetNumAllElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), p_mesh2->GetNumBoundaryElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumAllBoundaryElements(), p_mesh2->GetNumAllBoundaryElements());

            AbstractTetrahedralMesh<2,2>::ElementIterator iter2 = p_mesh2->GetElementIteratorBegin();

            for (AbstractTetrahedralMesh<2,2>::ElementIterator iter = p_mesh->GetElementIteratorBegin();
                 iter != p_mesh->GetElementIteratorEnd();
                 ++iter, ++iter2)
            {
                TS_ASSERT_EQUALS(iter->GetNumNodes(), iter2->GetNumNodes());
                for (unsigned i=0; i<iter->GetNumNodes(); i++)
                {
                    TS_ASSERT_EQUALS(iter->GetNodeGlobalIndex(i), iter2->GetNodeGlobalIndex(i));
                }
            }

            // We now need to free the mesh, since there is no honeycomb generator to do so.
            delete p_mesh2;
        }
    }

    void TestConstructFromNodeList()
    {
        std::vector<Node<2>*> nodes;

        nodes.push_back(new Node<2>(0, true, 0.1, -0.01));
        nodes.push_back(new Node<2>(1, true, 0.5, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.9, 0.01));
        nodes.push_back(new Node<2>(3, true, 0.1, 0.99));
        nodes.push_back(new Node<2>(4, true, 0.5, 1.0));
        nodes.push_back(new Node<2>(5, true, 0.9, 1.01));

        const double width = 1.0;
        Cylindrical2dMesh mesh(width, nodes);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6u);

        // Find the element with node indices 2,3,5 (which stradles the periodic boundary)
        unsigned element_index;
        std::set<unsigned> target_element_node_indices;
        for (element_index=0; element_index<mesh.GetNumElements(); element_index++)
        {
            target_element_node_indices.clear();
            target_element_node_indices.insert(2);
            target_element_node_indices.insert(3);
            target_element_node_indices.insert(5);

            for (unsigned node_local_index=0; node_local_index<=2; node_local_index++)
            {
                target_element_node_indices.erase(mesh.GetElement(element_index)->GetNodeGlobalIndex(node_local_index));
            }
            if (target_element_node_indices.empty())
            {
                break;
            }
        }
        TS_ASSERT_EQUALS(target_element_node_indices.empty(), true);

        // Calculate the circumsphere of the element
//        c_vector<double, 3> circumsphere = mesh.GetElement(element_index)->CalculateCircumsphere();
//
//        TS_ASSERT_DELTA(circumsphere[0], 0.9509, 1e-3);
//        TS_ASSERT_DELTA(circumsphere[1], 0.5100, 1e-3);
//        TS_ASSERT_DELTA(circumsphere[2], 0.2526, 1e-3);

        /* The reason that the circumsphere is calculated correctly for a periodic boundary
         * stradling element is somewhat obscure.
         * The Jacobian of the element is calculated when the element has a mirror node
         * The mirror node is then replaced with the node within the periodic mesh
         * The circumsphere is calculated based on the Jacobian and the replaced node within the periodic mesh
         *
         * uncommenting the following line of code causes an error:
         * Jacobian determinant is non-positive
         *
         * mesh.GetElement(element_index)->RefreshJacobianDeterminant();
         */
    }

    void TestGenerateVectorsOfElementsStraddlingPeriodicBoundaries()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/bad_cylindrical_9_1");
        Cylindrical2dMesh mesh(9.1);
        mesh.ConstructFromMeshReader(mesh_reader);

        // We now emulate the commands of the ReMesh function as far as it goes before generating the lists
        {
            mesh.CreateHaloNodes();
            mesh.CreateMirrorNodes();

            NodeMap big_map(mesh.GetNumAllNodes());
            mesh.MutableMesh<2,2>::ReMesh(big_map);
        }

        mesh.GenerateVectorsOfElementsStraddlingPeriodicBoundaries();

        TS_ASSERT_EQUALS(mesh.mLeftPeriodicBoundaryElementIndices.size(), 43u);

        // The commented test below fails as there are no nodes waiting to be deleted
        // and the current mesh is Voronoi, hence no call is made to triangle...
//        TS_ASSERT_EQUALS(mesh.mRightPeriodicBoundaryElementIndices.size(), 42u);

        // Test the GetCorrespondingNodeIndex() method

        // ... instead, there should still be the same number
        TS_ASSERT_EQUALS(mesh.mRightPeriodicBoundaryElementIndices.size(), 43u);

        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(393), 84u);
        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(84), 393u);
        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(188), 329u);
        TS_ASSERT_EQUALS(mesh.GetCorrespondingNodeIndex(329), 188u);

        mesh.CorrectNonPeriodicMesh();

        mesh.DeleteHaloNodes();
    }

    void TestCorrectNonPeriodicMeshMapLeftToRight()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/bad_cylindrical_9_1");
        Cylindrical2dMesh mesh(9.1);
        mesh.ConstructFromMeshReader(mesh_reader);

        NodeMap map(0);
        mesh.ReMesh(map);
        assert(map.IsIdentityMap());

        for (unsigned node_index=0; node_index<mesh.GetNumAllNodes(); node_index++)
        {
            std::vector<unsigned> indices;

            // Get the forward star from each node that isn't at the top or bottom boundary
            Node<2>* p_node = mesh.GetNode(node_index);
            if (p_node->rGetLocation()[1] < -2.5)
            {
                continue;
            }
            if (p_node->rGetLocation()[1] > 13.8)
            {
                continue;
            }

            // Iterate over countaining elements to get the elements of the forward star
            for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
                it != p_node->ContainingElementsEnd();
                ++it)
            {
                Element <2,2>* p_element = mesh.GetElement(*it);
                for (unsigned j=0; j<3; j++)
                {
                    unsigned index=p_element->GetNodeGlobalIndex(j);
                    if (index != node_index)
                    {
                        indices.push_back(index);
                    }
                }
            }

            // Each node in the forward star should appear exactly twice.  Sort and test.
            sort(indices.begin(), indices.end());
            for (unsigned i=0; i<indices.size(); i++)
            {
                if (i%2 == 0)
                {
                   TS_ASSERT_EQUALS(indices[i], indices[i+1]);
                }
            }
        }
    }

    void TestCorrectNonPeriodicMeshes()
    {
        std::vector<Node<2>*> nodes;
        // Generates a mesh which could be meshed in different ways.
        nodes.push_back(new Node<2>(0, true, 1.1, 0.0));
        nodes.push_back(new Node<2>(1, true, 3.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.1, 2.0));    // Stabilise mesh and prevent extra edge elements
        nodes.push_back(new Node<2>(3, true, 2.9, 2.0));
        nodes.push_back(new Node<2>(4, true, 1.0, 4.0));
        nodes.push_back(new Node<2>(5, true, 3.0, 4.0));

        const double width = 4.0;
        Cylindrical2dMesh mesh(width, nodes);
        // Create the mirrored nodes - double the size of the mesh
        mesh.CreateMirrorNodes();

        // Create elements for the new larger mesh
        NodeMap big_map(mesh.GetNumAllNodes());
        mesh.MutableMesh<2,2>::ReMesh(big_map);

        // We need the mesh in a certain configuration for this test
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 12u);

        mesh.GenerateVectorsOfElementsStraddlingPeriodicBoundaries();
        std::set<unsigned> left_elements = mesh.mLeftPeriodicBoundaryElementIndices;
        std::set<unsigned> right_elements = mesh.mRightPeriodicBoundaryElementIndices;
        // There should be four elements on each side which cross the periodic boundary
        TS_ASSERT_EQUALS(left_elements.size(), 4u);
        TS_ASSERT_EQUALS(right_elements.size(), 4u);

        /*
         * Swap one of the pairs of elements on the left around so that
         * CorrectNonPeriodicMesh() has some work to do.
         * Note this test could possibly make the mesh break the Voronoi condition,
         * but this is OK as the CorrectNonPeriodicMesh()
         * deals with cases where the Voronoi definition is ambiguous.
         */

        // A pair of elements on the left are flipped around.
        mesh.GetElement(0)->UpdateNode(0, mesh.GetNode(0));
        mesh.GetElement(1)->UpdateNode(1, mesh.GetNode(10));
        mesh.CorrectNonPeriodicMesh();

        // Check that these elements have been swapped around.
        // Note that the element indices are changed by each call and you
        // have to reexamine the mesh if changing this test
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(1)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(2)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(mesh.GetElement(11)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(11)->GetNode(1)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(mesh.GetElement(11)->GetNode(2)->GetIndex(), 0u);

        // A pair of elements on the right are flipped around.
        mesh.GetElement(4)->UpdateNode(0, mesh.GetNode(6));
        mesh.GetElement(6)->UpdateNode(1, mesh.GetNode(3));
        mesh.CorrectNonPeriodicMesh();

        // Check that these elements have been swapped around.
        // Note that the element indices are changed by each call and you
        // have to reexamine the mesh if changing this test
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(1)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(2)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(mesh.GetElement(11)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(11)->GetNode(1)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(mesh.GetElement(11)->GetNode(2)->GetIndex(), 0u);

        // We can now reconstruct the cylindrical mesh without any problems
        TS_ASSERT_THROWS_NOTHING(mesh.ReconstructCylindricalMesh());
    }

    void TestVoronoiTessellationUsesOverriddenMetric()
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        TS_ASSERT_EQUALS(p_mesh->CheckIsVoronoi(), true);

        // Create Voronoi tessellation
        VertexMesh<2, 2> tessellation(*p_mesh);

        //  Get two neighbouring nodes on boundary 48 and 53.
        //  Check that they have a common edge
        //  check it is a reasonable length (O(1)?)

        c_vector<double, 2> location_48 = p_mesh->GetNode(48)->rGetLocation();
        double common_edge_between_48_and_53 = tessellation.GetEdgeLength(48, 53);

        TS_ASSERT_DELTA(tessellation.GetEdgeLength(48, 49), pow(3.0, -0.5), 1e-4);

        TS_ASSERT_DELTA(common_edge_between_48_and_53,  pow(3.0, -0.5), 1e-4);

        //  Check that both cells have a reasonable sized area
        TS_ASSERT_DELTA(tessellation.GetVolumeOfElement(44),  0.5 * pow(3.0, 0.5), 1e-4);
        TS_ASSERT_DELTA(tessellation.GetSurfaceAreaOfElement(44), 2 * pow(3.0, 0.5), 1e-4);

        TS_ASSERT_DELTA(tessellation.GetVolumeOfElement(48),  0.5 * pow(3.0, 0.5), 1e-4);
        TS_ASSERT_DELTA(tessellation.GetSurfaceAreaOfElement(48), 2 * pow(3.0, 0.5), 1e-4);
    }
};

#endif /*TESTCYLINDRICAL2DMESH_HPP_*/
