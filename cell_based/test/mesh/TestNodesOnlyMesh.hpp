/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef TESTNODESONLYMESH_HPP_
#define TESTNODESONLYMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

#include "UblasCustomFunctions.hpp"
#include "NodesOnlyMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "ArchiveOpener.hpp"

class TestNodesOnlyMesh : public CxxTest::TestSuite
{
public:

    void TestConstructNodesWithoutMesh()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetGlobalNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), 8u);

        // Check that the nodes mapping is set up correctly
        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT(!(p_mesh->mNodesMapping.find(i) == p_mesh->mNodesMapping.end()));
            TS_ASSERT_EQUALS(p_mesh->SolveNodeMapping(i), p_mesh->mNodesMapping[i]);
        }

        TS_ASSERT_THROWS_CONTAINS(p_mesh->SolveNodeMapping(8), " does not belong to process ")

        // Avoid memory leak
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCalculateBoundingBox() throw (Exception)
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, -1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, -1.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));

        NodesOnlyMesh<3> mesh;
        ChasteCuboid<3> bounding_cuboid = mesh.CalculateBoundingBox(nodes);

        TS_ASSERT_DELTA(bounding_cuboid.rGetUpperCorner()[0], 1.0, 1e-10);
        TS_ASSERT_DELTA(bounding_cuboid.rGetUpperCorner()[1], 1.0, 1e-10);
        TS_ASSERT_DELTA(bounding_cuboid.rGetUpperCorner()[2], 0.0, 1e-10);

        TS_ASSERT_DELTA(bounding_cuboid.rGetLowerCorner()[0], -1.0, 1e-10);
        TS_ASSERT_DELTA(bounding_cuboid.rGetLowerCorner()[1], -1.0, 1e-10);
        TS_ASSERT_DELTA(bounding_cuboid.rGetLowerCorner()[2], 0.0, 1e-10);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestConstuctingAndSwellingInitialBoxCollection() throw (Exception)
    {
        double cut_off = 0.5;

        /*
         * Nodes chosen so to test the cases that the domain width in x is
         * "divisible" by the cut_off, the y-dimension is not "divisible" and
         * the z-direction has zero extent.
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, -1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, -1.1, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.1, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, cut_off);

        BoxCollection<3>* p_box_collection = mesh.mpBoxCollection;

        TS_ASSERT(p_box_collection != NULL);

        // 5x5x1 box collection
        TS_ASSERT_EQUALS(p_box_collection->GetNumBoxes(), 25u);

        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(mesh.GetNode(0)), 10u);
        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(mesh.GetNode(1)), 4u);
        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(mesh.GetNode(2)), 22u);
        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(mesh.GetNode(3)), 24u);

        mesh.EnlargeBoxCollection();

        BoxCollection<3>* p_new_collection = mesh.mpBoxCollection;

        // 7x7x3 box collection
        TS_ASSERT_EQUALS(p_new_collection->GetNumBoxes(), 147u);

        // The "old" boxes should be in the same spatial position.
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(mesh.GetNode(0)), 71u);
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(mesh.GetNode(1)), 61u);
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(mesh.GetNode(2)), 87u);
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(mesh.GetNode(3)), 89u);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCheckingBoxCollectionSizeAndSwelling() throw (Exception)
    {
        double cut_off = 0.5;

        /*
         * Nodes chosen so to test the cases that the domain width in x is
         * "divisible" by the cut_off, the y-dimension is not "divisible" and
         * the z-direction has zero extent.
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, -0.9, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.9, -0.9, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 0.9, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.9, 0.9, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, cut_off);

        TS_ASSERT_EQUALS(mesh.mpBoxCollection->GetNumBoxes(), 16u);

        bool enlarge = mesh.IsANodeCloseToDomainBoundary();

        TS_ASSERT(enlarge);

        // Enlarge twice to make domain large enough in z-direction for the enlarge test to fail.
        mesh.EnlargeBoxCollection();
        mesh.EnlargeBoxCollection();

        enlarge = mesh.IsANodeCloseToDomainBoundary();

        TS_ASSERT(!enlarge);

        mesh.SetMinimumNodeDomainBoundarySeparation(2.0);

        enlarge = mesh.IsANodeCloseToDomainBoundary();

        TS_ASSERT(enlarge);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestClearingNodesOnlyMesh()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 8u);

        p_mesh->Clear();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllBoundaryElements(), 0u);
        TS_ASSERT(!p_mesh->mpBoxCollection);

        // Avoid memory leak
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestConstructNodesWithoutMeshUsingMesh()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(generating_mesh, 1.5);

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllNodes(), 5u);

        // Avoid memory leak
        delete p_mesh;
    }

    void TestGetNextAvailableIndex()
    {
        // Construct a simple 2D mesh.
        ChastePoint<2> point1(0.0, 0.0);
        ChastePoint<2> point2(2.9, 2.9);
        ChasteCuboid<2> cuboid(point1, point2);

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(1, false, 1.5, 1.5));
        nodes.push_back(new Node<2>(2, false, 2.5, 2.5));
        nodes.push_back(new Node<2>(3, false, 1.5, 2.5));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // These resuls will change when the nodes are distributed.
        if(PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 4u);
        }
        else if(PetscTools::GetNumProcs() == 2)
        {
            if(PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 4u);
            }
            else if(PetscTools::GetMyRank() == 1)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 5u);
            }
        }
        else if(PetscTools::GetNumProcs() == 3)
        {
            if(PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 12u);
            }
            else if(PetscTools::GetMyRank() == 1)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 13u);
            }
            else if(PetscTools::GetMyRank() == 2)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 14u);
            }
        }

//        // Delete some nodes and make sure that their global indices become available.
//        if(mesh.GetNodeOwnership(0))
//        {
//            mesh.DeleteNode(0);
//            TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 0u);
//        }
//        else // Covers an exception
//        {
//            TS_ASSERT_THROWS_CONTAINS(mesh.GetNode(0), "Requested node 0 does not belong to processor ");
//        }
//
//        if(mesh.GetNodeOwnership(1))
//        {
//            mesh.DeleteNode(1);
//            unsigned next_index = mesh.GetNextAvailableIndex();
//            TS_ASSERT((next_index == 0u) || ( next_index == 1u) );
//        }
//
//        // Delete a node as if it moves off this process, and make sure its global index is not available
//        if(mesh.GetNodeOwnership(2))
//        {
//            mesh.DeleteMovedNode(2);
//            TS_ASSERT(mesh.GetNextAvailableIndex() != 1u);
//        }

        //Clean up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestWriteNodesWithoutMeshUsingVtk()
    {
 #ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        // Note in my version of Paraview, you need data on points before you can view with Glyphs
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "just_nodes", false);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            distance.push_back(norm_2(p_mesh->GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add boundary node "point" data
        std::vector<double> boundary;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            boundary.push_back(p_mesh->GetNode(i)->IsBoundaryNode());
        }
        writer.AddPointData("Boundary", boundary);

        // Add fibre type to "point" data
        std::vector< c_vector<double, 3> > location;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            location.push_back(p_mesh->GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        writer.WriteFilesUsingMesh(*p_mesh);

        // Avoid memory leak
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
    }

    void TestWriteNodesWithoutMesh()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        TrianglesMeshWriter<3,3> writer("TestMeshWriter", "3dNodesOnlyMesh");
        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(*p_mesh));

        // Avoid memory leak
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestGetSetMethods()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));

        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        p_mesh->GetNode(0)->SetRadius(1.0);
        p_mesh->GetNode(1)->SetRadius(2.0);

        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetRadius(), 1.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(1)->GetRadius(), 2.0, 1e-6);

        // Avoid memory leak
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestAddNode() throw (Exception)
    {
        std::vector<Node<2>*> nodes;
        Node<2> node0(0, true, 0.0, 0.0);
        nodes.push_back(&node0);
        Node<2> node1(1, true, 0.0, 0.5);
        nodes.push_back(&node1);

        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        // Test add node
        p_mesh->AddNode(new Node<2>(2, true, 0.0, 1.0));//This node pointer is added to the mesh and deleted by the destructor

        unsigned num_nodes = 3;
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), num_nodes);
        TS_ASSERT_DELTA(p_mesh->GetNode(2)->GetRadius(), 0.5, 1e-4);

        // Cover an exception.
        TS_ASSERT_THROWS_CONTAINS(p_mesh->GetNode(3)->SetRadius(1.0), " does not belong to process ");

        // Avoid memory leak
        delete p_mesh;
    }

    void TestDeleteNodesAndRemesh() throw (Exception)
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.0, 0.5));
        nodes.push_back(new Node<2>(2, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(5, false, 0.4, 0.5));
        nodes.push_back(new Node<2>(6, false, 0.6, 0.5));
        nodes.push_back(new Node<2>(7, false, 0.0, 0.5));

        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(nodes, 1.5);

        // Test that there are never any boundary nodes
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryNodes(), 0u);
        NodeMap node_map(8);
        p_mesh->ReMesh(node_map);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryNodes(), 0u);

        // Free memory - the constructor does a deep copy of its input
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

        // Set radii of cells from 1 to 8
        for (unsigned i=0; i<nodes.size(); i++)
        {
            p_mesh->GetNode(i)->SetRadius(i+1);
        }

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 8u);

        // Delete from interior
        TS_ASSERT_DELTA(p_mesh->GetNode(6)->GetRadius(), 7.0, 1e-4);
        p_mesh->DeleteNode(6);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 7u);

        // Delete from edge
        p_mesh->DeleteNode(1);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 6u);

        // Delete from corner
        TS_ASSERT_DELTA(p_mesh->GetNode(3)->GetRadius(), 4.0, 1e-4);
        p_mesh->DeleteNode(3);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 5u);

        // Deleting a deleted node should throw an exception
        TS_ASSERT_THROWS_THIS(p_mesh->DeleteNode(3),"Trying to delete a deleted node");

        /*
         * Check that the cell radii are updated correctly when a new cell
         * is added using the most recently deleted index.
         * (Index 3 is at the back of the deleted nodes list and is thus the one to be reused.)
         */
        p_mesh->AddNode(new Node<2>(0, true, 6.0, 6.0)); //This node pointer is added to the mesh and deleted by the destructor

        // Check the most recently deleted node now has the correct cell radius
        TS_ASSERT_DELTA(p_mesh->GetNode(3)->GetRadius(), 0.5, 1e-4);

        // Now we have deleted/reused 3, and deleted 1 and 6.
        // The new nodes are:
        // New:     0     1     2     3     4     5
        // Old:     0     2 (new3)    4     5     7
        // Radius:  1     3     0.5     5     6     8

        NodeMap map(8);
        p_mesh->ReMesh(map);
        TS_ASSERT_EQUALS(map.GetNewIndex(0), 0u);
        TS_ASSERT(map.IsDeleted(1));
        TS_ASSERT_EQUALS(map.GetNewIndex(2), 1u);
        TS_ASSERT_EQUALS(map.GetNewIndex(3), 2u);
        TS_ASSERT_EQUALS(map.GetNewIndex(4), 3u);
        TS_ASSERT_EQUALS(map.GetNewIndex(5), 4u);
        TS_ASSERT(map.IsDeleted(6));
        TS_ASSERT_EQUALS(map.GetNewIndex(7), 5u);

        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetRadius(), 1.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(1)->GetRadius(), 3.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(2)->GetRadius(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(3)->GetRadius(), 5.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(4)->GetRadius(), 6.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->GetRadius(), 8.0, 1e-4);

        // Avoid memory leak
        delete p_mesh;
    }

    void TestArchiving() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "nodes_only_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("nodes_only_mesh");

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
            TetrahedralMesh<2,2> generating_mesh;
            generating_mesh.ConstructFromMeshReader(mesh_reader);

            // Convert this to a NodesOnlyMesh
            NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
            p_mesh->ConstructNodesWithoutMesh(generating_mesh, 1.5);

            p_mesh->GetNode(0)->SetRadius(1.12);
            p_mesh->GetNode(1)->SetRadius(2.34);

            TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetRadius(), 1.12, 1e-6);
            TS_ASSERT_DELTA(p_mesh->GetNode(1)->GetRadius(), 2.34, 1e-6);

            AbstractTetrahedralMesh<2,2>* const p_const_mesh = p_mesh;

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) << p_const_mesh;

            // Avoid memory leak
            delete p_mesh;
        }

        {
            /*
             * Should archive the most abstract class possible to check that
             * boost knows what individual classes are (but here AbstractMesh
             * doesn't have the methods below).
             */
            AbstractTetrahedralMesh<2,2>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            NodesOnlyMesh<2>* p_nodes_only_mesh = dynamic_cast<NodesOnlyMesh<2>*>(p_mesh2);

            // Check we have the right number of nodes & elements
            TS_ASSERT_EQUALS(p_nodes_only_mesh->GetNumNodes(), 543u);
            TS_ASSERT_EQUALS(p_nodes_only_mesh->GetNumElements(), 0u);

            // Check some node co-ordinates
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(1)->GetPoint()[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(1)->GetPoint()[1], 0.0, 1e-6);

            // Check some cell radii
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(0)->GetRadius(), 1.12, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(1)->GetRadius(), 2.34, 1e-6);

            // Tidy up
            delete p_mesh2;
        }
    }

    void TestArchiveNodesOnlyMeshWithNodeAttributes() throw (Exception)
    {
        FileFinder archive_dir("archive_mutable_mesh", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mutable_mesh_with_attributes.arch";

        {
            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            std::vector<Node<3> *> nodes;

            nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
            nodes.push_back(new Node<3>(1, true,  0.0,  1.0,  0.0));
            nodes.push_back(new Node<3>(2, true,  0.0,  0.0,  1.0));
            nodes.push_back(new Node<3>(3, true,  1.0,  0.0,  0.0));

            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);

            for (unsigned i=0; i<4; i++)
            {
                mesh.GetNode(i)->SetRadius(1.2);
                mesh.GetNode(i)->SetIsParticle(true);
            }

            AbstractTetrahedralMesh<3,3>* const p_mesh = &mesh;

            (*p_arch) << p_mesh;

            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }
        }

        {
            AbstractTetrahedralMesh<3,3>* p_mesh;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_mesh;

            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 4u);

            for (unsigned i=0; i<4; i++)
            {
                TS_ASSERT_DELTA(p_mesh->GetNode(i)->GetRadius(), 1.2, 1e-15);
                TS_ASSERT(p_mesh->GetNode(i)->IsParticle());
            }

            delete p_mesh;
        }
    }
};

#endif /*TESTNODESONLYMESH_HPP_*/
