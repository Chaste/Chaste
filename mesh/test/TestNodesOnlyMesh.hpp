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

#ifndef TESTNODESONLYMESH_HPP_
#define TESTNODESONLYMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

#include "SmartPointers.hpp"
#include "UblasCustomFunctions.hpp"
#include "NodesOnlyMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "ArchiveOpener.hpp"
#include "PetscSetupAndFinalize.hpp"

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

        // For coverage, give an attribute to one of the nodes
        nodes[0]->AddNodeAttribute(9.81);

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        if (PetscTools::IsSequential())
        {
            // Check that node 0 has the correct attribute
            TS_ASSERT_EQUALS(mesh.GetNode(0)->HasNodeAttributes(), true);
            TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetNodeAttributes().size(), 1u);
            TS_ASSERT_DELTA(mesh.GetNode(0)->rGetNodeAttributes()[0], 9.81, 1e-4);
        }

        unsigned num_nodes = PetscTools::AmMaster() ? 8 : 0;    // All nodes will lie on the master process.
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), num_nodes);

        // Check that the nodes mapping is set up correctly
        if (PetscTools::AmMaster())
        {
            for (unsigned i=0; i<nodes.size(); i+=PetscTools::GetNumProcs())
            {
                TS_ASSERT(!(mesh.mNodesMapping.find(i) == mesh.mNodesMapping.end()));
                TS_ASSERT_EQUALS(mesh.SolveNodeMapping(i), mesh.mNodesMapping[i]);
            }
        }

        if (PetscTools::IsSequential())
        {
            /*
             * Check that other nodes do not have attributes.
             *
             * Note: since the radius of each node is set to 0.5 in
             * NodesOnlyMesh::ConstructNodesWithoutMesh(), this means
             * that every node in a NodeBasedCellPopulaton has called
             * ConstructNodeAttributes(); however, rGetNodeAttributes()
             * will return an empty vector unless any attributes have
             * been set by the user.
             */
            for (unsigned i=1; i<7; i++)
            {
                TS_ASSERT_EQUALS(mesh.GetNode(i)->HasNodeAttributes(), true);
                TS_ASSERT_EQUALS(mesh.GetNode(i)->rGetNodeAttributes().size(), 0u);
            }
        }

        TS_ASSERT_THROWS_CONTAINS(mesh.SolveNodeMapping(8*PetscTools::GetNumProcs() + PetscTools::GetMyRank() + 1), " does not belong to process ");

        // Cover being able to set the maximum interation distance.
        mesh.SetMaximumInteractionDistance(5.0);
        TS_ASSERT_DELTA(mesh.GetMaximumInteractionDistance(), 5.0, 1e-4);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestConstructNodesWithoutMeshSharedPtr()
    {
        std::vector<boost::shared_ptr<Node<3> > > nodes;
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(0, true,  0.0, 0.0, 0.0)));
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(1, false, 1.0, 0.0, 0.0)));
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(2, false, 0.0, 1.0, 0.0)));
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(3, false, 1.0, 1.0, 0.0)));
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(4, false, 0.0, 0.0, 1.0)));
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(5, false, 1.0, 0.0, 1.0)));
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(6, false, 0.0, 1.0, 1.0)));
        nodes.push_back(boost::shared_ptr<Node<3> > (new Node<3>(7, false, 1.0, 1.0, 1.0)));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        unsigned num_nodes = PetscTools::AmMaster() ? 8 : 0;    // All nodes will lie on the master process.
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), num_nodes);

        // Check that the nodes mapping is set up correctly
        if (PetscTools::AmMaster())
        {
            for (unsigned i=0; i<nodes.size(); i+=PetscTools::GetNumProcs())
            {
                TS_ASSERT(!(mesh.mNodesMapping.find(i) == mesh.mNodesMapping.end()));
                TS_ASSERT_EQUALS(mesh.SolveNodeMapping(i), mesh.mNodesMapping[i]);
            }
        }

        TS_ASSERT_THROWS_CONTAINS(mesh.SolveNodeMapping(8*PetscTools::GetNumProcs() + PetscTools::GetMyRank() + 1), " does not belong to process ");

        // Cover being able to set the maximum interation distance.
        mesh.SetMaximumInteractionDistance(5.0);
        TS_ASSERT_DELTA(mesh.GetMaximumInteractionDistance(), 5.0, 1e-4);
    }

    void TestCalculateBoundingBox()
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

        // Make sure we get a zero-element box if there are no nodes.
        std::vector<Node<3>* > empty_nodes;
        NodesOnlyMesh<3> empty_mesh;

        ChasteCuboid<3> empty_bounding_cuboid = empty_mesh.CalculateBoundingBox(empty_nodes);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(empty_bounding_cuboid.rGetLowerCorner()[i], empty_bounding_cuboid.rGetUpperCorner()[i], 1e-4);
        }
    }

    void TestSetInitialBoxCollection()
    {
        double cut_off = 1.0;

        c_vector<double, 6> domain_size = zero_vector<double>(6);
        domain_size[1] = 2.9;
        domain_size[3] = 2.9;
        domain_size[5] = 2.9;

        NodesOnlyMesh<3> mesh;

        // Set the initial box collection.
        mesh.SetInitialBoxCollection(domain_size, cut_off);

        c_vector<double, 6> set_domain_size = mesh.GetBoxCollection()->rGetDomainSize();
        double box_width = mesh.GetBoxCollection()->GetBoxWidth();
        double expected_domain_size[6] = {0.0, 3.0, 0.0, 3.0, 0.0, 3.0};
        // Distributed boxes will extend the range of the z dimension when there are redundant processes
        expected_domain_size[5] = (double) std::max(3u, PetscTools::GetNumProcs());

        // Check the correct domain size has been set (swelled to be a multiple of the box width)
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(set_domain_size[2*i], 0.0, 1e-4);
            TS_ASSERT_DELTA(set_domain_size[2*i+1], expected_domain_size[2*i+1], 1e-4);
        }

        // Make sure the boxes are the width of the cut-off distance.
        TS_ASSERT_DELTA(box_width, cut_off, 1e-4);

        if (PetscTools::AmMaster())
        {
            c_vector<double,3> location = zero_vector<double>(3);
            location[2] = 0.1;    // This should be owned by process 0 in any space decomposition.
            TS_ASSERT(mesh.IsOwned(location));
        }
    }

    void TestConstuctingAndEnlargingInitialBoxCollection()
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

        DistributedBoxCollection<3>* p_box_collection = mesh.mpBoxCollection;

        TS_ASSERT(p_box_collection != NULL);

        // 5x5xnum_procs box collection
        unsigned num_boxes = 5*5*PetscTools::GetNumProcs();
        TS_ASSERT_EQUALS(p_box_collection->GetNumBoxes(), num_boxes);

        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(nodes[0]), 10u);
        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(nodes[1]), 4u);
        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(nodes[2]), 22u);
        TS_ASSERT_EQUALS(p_box_collection->CalculateContainingBox(nodes[3]), 24u);

        mesh.EnlargeBoxCollection();

        DistributedBoxCollection<3>* p_new_collection = mesh.mpBoxCollection;

        // 7x7xnum_procs box collection
        num_boxes = 7*7* (2 + PetscTools::GetNumProcs());
        TS_ASSERT_EQUALS(p_new_collection->GetNumBoxes(), num_boxes);

        // The "old" boxes should be in the same spatial position.
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(nodes[0]), 71u);
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(nodes[1]), 61u);
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(nodes[2]), 87u);
        TS_ASSERT_EQUALS(p_new_collection->CalculateContainingBox(nodes[3]), 89u);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestBoxCollectionSizeAndSwelling()
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

        TS_ASSERT_EQUALS(mesh.mpBoxCollection->GetNumBoxes(), 16*PetscTools::GetNumProcs());

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

    void TestMovingNodesInBoxCollection()
    {
        EXIT_IF_PARALLEL;    // This wont work until nodes can be re-assigned on re-mesh. #2260

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

        // 4 x 4 x 1 collection
        TS_ASSERT_EQUALS(mesh.mpBoxCollection->GetNumBoxes(), 16*PetscTools::GetNumProcs());

        // Test what happens if nodes move in initial collection.
        AbstractMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
        if (node_iter != mesh.GetNodeIteratorEnd())
        {
            c_vector<double, 3> initial_location = node_iter->rGetLocation();
            c_vector<double, 3> translation = scalar_vector<double>(3, -2.0*cut_off);
            ChastePoint<3> new_point(initial_location + translation);
            node_iter->SetPoint(new_point);
        }

        NodeMap map(4);
        TS_ASSERT_THROWS_NOTHING(mesh.ReMesh(map));

        mesh.ResizeBoxCollection();

        // 12 x 12 x 9
        TS_ASSERT_EQUALS(mesh.mpBoxCollection->GetNumBoxes(), 1296u);

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

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        unsigned num_nodes = PetscTools::AmMaster() ? 8 : 0;
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);

        mesh.Clear();

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), 0u);

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

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        unsigned num_boxes = PetscTools::AmMaster() ? 5 : 0;
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_boxes);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), num_boxes);

        // Coverage of rGetInitiallyOwnedNodes()
        std::vector<bool>& r_nodes = mesh.rGetInitiallyOwnedNodes();
        TS_ASSERT_EQUALS(r_nodes.size(), 5u);
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
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 4u);
        }
        else if (PetscTools::GetNumProcs() == 2)
        {
            if (PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 4u);
            }
            else if (PetscTools::GetMyRank() == 1)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 5u);
            }
        }
        else if (PetscTools::GetNumProcs() == 3)
        {
            if (PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 6u);
            }
            else if (PetscTools::GetMyRank() == 1)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 7u);
            }
            else if (PetscTools::GetMyRank() == 2)
            {
                TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), 2u);
            }
        }

        // Delete some nodes and make sure that their global indices become available.
        if (mesh.mpBoxCollection->IsOwned(nodes[0]))
        {
            mesh.DeleteNode(PetscTools::GetMyRank());
            TS_ASSERT_EQUALS(mesh.GetNextAvailableIndex(), PetscTools::GetMyRank());
        }
        else // Covers an exception
        {
            TS_ASSERT_THROWS_CONTAINS(mesh.GetNode(0), "Requested node 0 does not belong to process ");
        }

        // Clean up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestWriteNodesWithoutMeshUsingVtk()
    {
        EXIT_IF_PARALLEL;    // Cannot write to file yet in parallel.
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

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Note in my version of Paraview, you need data on points before you can view with Glyphs
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "just_nodes", false);

        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh.GetNode(i)->rGetLocation()));
        }
        writer.AddPointData("Distance from origin", distance);

        // Add boundary node "point" data
        std::vector<double> boundary;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            boundary.push_back(mesh.GetNode(i)->IsBoundaryNode());
        }
        writer.AddPointData("Boundary", boundary);

        // Add fibre type to "point" data
        std::vector< c_vector<double, 3> > location;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            location.push_back(mesh.GetNode(i)->rGetLocation());
        }
        writer.AddPointData("Location", location);

        writer.WriteFilesUsingMesh(mesh);

        // Avoid memory leak
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
        EXIT_IF_PARALLEL;    // Cannot write to file yet in parallel.
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        TrianglesMeshWriter<3,3> writer("TestMeshWriter", "3dNodesOnlyMesh");
        TS_ASSERT_THROWS_NOTHING(writer.WriteFilesUsingMesh(mesh));

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestGetSetRadiusMethods()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        if (PetscTools::AmMaster())
        {
            mesh.GetNode(0)->SetRadius(1.0);

            TS_ASSERT_DELTA(mesh.GetNode(0)->GetRadius(), 1.0, 1e-6);

        }

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestAddNode()
    {
        std::vector<Node<2>*> nodes;
        Node<2> node0(0, true, 0.0, 0.0);
        nodes.push_back(&node0);
        Node<2> node1(1, true, 0.0, 0.5);
        nodes.push_back(&node1);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        if (PetscTools::IsSequential())
        {
            // Test add node
            mesh.AddNode(new Node<2>(2, true, 0.0, 0.25)); // This node pointer is added to the mesh and deleted by the destructor

            unsigned num_nodes = 3;
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
            TS_ASSERT_DELTA(mesh.GetNode(2)->GetRadius(), 0.5, 1e-4);

            // Cover an exception.
            TS_ASSERT_THROWS_CONTAINS(mesh.GetNode(3)->SetRadius(1.0), " does not belong to process ");
        }
        else
        {
            if (PetscTools::GetMyRank() == 1)
            {
                mesh.AddNode(new Node<2>(2, true, 0.0, 2.0)); // This node pointer is added to the mesh and deleted by the destructor
            }

            unsigned num_nodes = PetscTools::AmMaster() ? 2 : ((PetscTools::GetMyRank() == 1) ? 1 : 0);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);

            if (PetscTools::GetMyRank() == 1)
            {
                TS_ASSERT_DELTA(mesh.GetNode(1)->GetRadius(), 0.5, 1e-4);
            }
        }
    }

    void TestDeleteNodesAndRemesh()
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

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Test that there are never any boundary nodes
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 0u);
        NodeMap node_map(8);
        mesh.ReMesh(node_map);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 0u);

        // Set radii of cells from 1 to 8
        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            node_iter->SetRadius(node_index+1);
        }

        // Delete from interior
        unsigned node_index;
        if (mesh.mpBoxCollection->IsOwned(nodes[0]))
        {
            AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
            node_index = node_iter->GetIndex();

            TS_ASSERT_DELTA(mesh.GetNode(node_index)->GetRadius(), (double)(node_index+1), 1e-4);
            mesh.DeleteNode(node_index);
            TS_ASSERT_THROWS_THIS(mesh.DeleteNode(node_index), "Trying to delete a deleted node");
        }

        // Global node indices should stay the same after remesh
        NodeMap map(mesh.GetMaximumNodeIndex());
        mesh.ReMesh(map);

        if (PetscTools::AmMaster())
        {
            TS_ASSERT(map.IsDeleted(node_index));
            TS_ASSERT_EQUALS(map.GetNewIndex(1), 1u);
            TS_ASSERT_EQUALS(map.GetNewIndex(2), 2u);
            TS_ASSERT_EQUALS(map.GetNewIndex(3), 3u);
            TS_ASSERT_EQUALS(map.GetNewIndex(4), 4u);
            TS_ASSERT_EQUALS(map.GetNewIndex(5), 5u);
            TS_ASSERT_EQUALS(map.GetNewIndex(7), 7u);

            // But local indices (the location in mNodes) should change).
            unsigned base_index = PetscTools::GetNumProcs();
            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(base_index), 0u);
            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(2*base_index), 1u);
            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(3*base_index), 2u);
            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(4*base_index), 3u);
            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(5*base_index), 4u);
            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(6*base_index), 5u);
        }

        TS_ASSERT_THROWS_CONTAINS(mesh.SolveNodeMapping(0), " does not belong to process ");

        if (PetscTools::AmMaster())
        {
            unsigned base_index = PetscTools::GetNumProcs();
            TS_ASSERT_DELTA(mesh.GetNode(base_index)->GetRadius(),   (double)base_index + 1.0, 1e-4);
            TS_ASSERT_DELTA(mesh.GetNode(2*base_index)->GetRadius(), (double)2*base_index + 1.0, 1e-4);
            TS_ASSERT_DELTA(mesh.GetNode(3*base_index)->GetRadius(), (double)3*base_index + 1.0, 1e-4);
            TS_ASSERT_DELTA(mesh.GetNode(4*base_index)->GetRadius(), (double)4*base_index + 1.0, 1e-4);
            TS_ASSERT_DELTA(mesh.GetNode(5*base_index)->GetRadius(), (double)5*base_index + 1.0, 1e-4);
            TS_ASSERT_DELTA(mesh.GetNode(6*base_index)->GetRadius(), (double)6*base_index + 1.0, 1e-4);
        }

        // Free memory - the constructor does a deep copy of its input
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCleanDeleteAndAddNode()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.0, 1.6));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        unsigned num_initial_nodes = mesh.GetNumNodes();

        if (PetscTools::AmMaster())
        {
            mesh.DeleteMovedNode(0);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_initial_nodes -1);
            TS_ASSERT(mesh.mDeletedGlobalNodeIndices.empty());
            TS_ASSERT(mesh.mDeletedNodeIndices.size() == 1u);

            boost::shared_ptr<Node<2> > p_node(new Node<2>(1));
            p_node->AddNodeAttribute(1.0);
            mesh.AddMovedNode(p_node);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_initial_nodes);
            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(1), 0u);
            TS_ASSERT(mesh.mDeletedNodeIndices.size() == 0u);
        }

        mesh.ClearBoxCollection();

        TS_ASSERT(!mesh.mpBoxCollection);

        // Clean up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestHaloNodes()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.0, 1.6));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        boost::shared_ptr<Node<2> > p_node(new Node<2>(3, false));
        mesh.AddHaloNode(p_node);

        TS_ASSERT_EQUALS(mesh.mHaloNodes.size(), 1u);
        TS_ASSERT_EQUALS(mesh.mHaloNodesMapping[3], 0u);

        mesh.ClearHaloNodes();

        TS_ASSERT_EQUALS(mesh.mHaloNodes.size(), 0u);

        delete nodes[0];
        delete nodes[1];
    }

    void TestGetNodesOutsideLocalDomain()
    {
        if (PetscTools::GetNumProcs() > 1)
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.0, 1.6));

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);

            if (PetscTools::AmMaster())
            {
                c_vector<double, 2> new_location = zero_vector<double>(2);
                new_location[1] = 1.7;
                ChastePoint<2> new_point(new_location);

                mesh.GetNode(0)->SetPoint(new_point);
            }
            if (PetscTools::GetMyRank() == 1)
            {
                c_vector<double, 2> new_location = zero_vector<double>(2);
                new_location[1] = 0.5;
                ChastePoint<2> new_point(new_location);

                mesh.GetNode(1)->SetPoint(new_point);
            }

            mesh.CalculateNodesOutsideLocalDomain();

            std::vector<unsigned> nodes_left = mesh.rGetNodesToSendLeft();
            std::vector<unsigned> nodes_right = mesh.rGetNodesToSendRight();

            if (PetscTools::AmMaster())
            {
                TS_ASSERT(nodes_left.empty());
                TS_ASSERT_EQUALS(nodes_right.size(), 1u);
                TS_ASSERT_EQUALS(nodes_right[0], 0u);
            }
            if (PetscTools::GetMyRank()==1)
            {
                TS_ASSERT(nodes_right.empty());
                TS_ASSERT_EQUALS(nodes_left.size(), 1u);
                TS_ASSERT_EQUALS(nodes_left[0], 1u);
            }
        }
    }

    void TestAddAndGetHaloNodes()
    {
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.0);

        unsigned num_original_nodes = mesh.GetNumNodes();
        boost::shared_ptr<Node<2> > p_halo_node(new Node<2>(10, false, 0.0, 1.0));
        mesh.AddHaloNode(p_halo_node);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_original_nodes);

        if (num_original_nodes > 0)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(0), mesh.GetNodeOrHaloNode(0));
        }
        TS_ASSERT_EQUALS(mesh.mHaloNodes.size(), 1u);
        TS_ASSERT_THROWS_NOTHING(mesh.GetNodeOrHaloNode(10));
        TS_ASSERT_EQUALS(mesh.GetNodeOrHaloNode(10)->GetIndex(), 10u);
        TS_ASSERT_DELTA(mesh.GetNodeOrHaloNode(10)->rGetLocation()[1], 1.0, 1e-4);

        delete nodes[0];
    }

    void TestArchiving()
    {
        EXIT_IF_PARALLEL;    ///\todo parallel archiving not yet possible.

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "nodes_only_mesh.arch";
        ArchiveLocationInfo::SetMeshFilename("nodes_only_mesh");

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
            TetrahedralMesh<2,2> generating_mesh;
            generating_mesh.ConstructFromMeshReader(mesh_reader);

            // Convert this to a NodesOnlyMesh
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543u);
            TS_ASSERT_EQUALS(mesh.GetAllNodeIndices().size(), 543u);
            TS_ASSERT_EQUALS(mesh.GetAllNodeIndices()[2], 2u);
            mesh.GetNode(0)->SetRadius(1.12);
            mesh.GetNode(1)->SetRadius(2.34);

            // Delete a node in order to test re-indexing (nodes 3,4,5... should still be named 3,4,5...)
            mesh.DeleteNode(2);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 542u);
            TS_ASSERT_EQUALS(mesh.GetAllNodeIndices().size(), 542u);
            TS_ASSERT_EQUALS(mesh.GetAllNodeIndices()[2], 3u);

            TS_ASSERT_DELTA(mesh.GetNode(0)->GetRadius(), 1.12, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1)->GetRadius(), 2.34, 1e-6);

            AbstractTetrahedralMesh<2,2>* const p_const_mesh = &mesh;

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) << p_const_mesh;
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
            TS_ASSERT_EQUALS(p_nodes_only_mesh->GetNumNodes(), 543u - 1u);
            TS_ASSERT_EQUALS(p_nodes_only_mesh->GetNumElements(), 0u);

            // Check some node co-ordinates
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(1)->GetPoint()[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(1)->GetPoint()[1], 0.0, 1e-6);

            // Node with index 2 is not available
            TS_ASSERT_THROWS_CONTAINS(p_nodes_only_mesh->GetNode(2), "Requested node 2 does not belong to process ");
            // Test re-indexing - Node with index 3 still has index 3
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(3)->GetPoint()[0], 0.992114, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(3)->GetPoint()[1], 0.125333, 1e-6);

            // Check some cell radii
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(0)->GetRadius(), 1.12, 1e-6);
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(1)->GetRadius(), 2.34, 1e-6);

            // Test that a newly added node gets the correct index (543)
            TS_ASSERT_THROWS_CONTAINS(p_nodes_only_mesh->GetNode(543), "Requested node 543 does not belong to process ");
            // Node that the index "2" here is fake.  The new node gets a fresh index
            p_nodes_only_mesh->AddNode(new Node<2>(2, true, 0.12345678, 2.0)); // This node pointer is added to the mesh and deleted by the destructor
            // Node with index 2 is still not available
            TS_ASSERT_THROWS_CONTAINS(p_nodes_only_mesh->GetNode(2), "Requested node 2 does not belong to process ");
            TS_ASSERT_DELTA(p_nodes_only_mesh->GetNode(543)->GetPoint()[0], 0.12345678, 1e-8);
            // Tidy up
            delete p_mesh2;
        }
    }

    void TestArchiveNodesOnlyMeshWithNodeAttributes()
    {
        EXIT_IF_PARALLEL;    ///\todo parallel archiving not yet possible.

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

    void TestLoadBalanceMesh()
    {
        // Test designed for np=2
        if (PetscTools::GetNumProcs() == 2)
        {
            std::vector<Node<1>*> nodes;
            nodes.push_back(new Node<1>(0, true,  0.0));
            nodes.push_back(new Node<1>(1, false, 1.0));
            nodes.push_back(new Node<1>(2, false, 2.0));
            nodes.push_back(new Node<1>(3, false, 3.0));
            nodes.push_back(new Node<1>(4, false, 4.0));
            nodes.push_back(new Node<1>(5, false, 5.0));
            nodes.push_back(new Node<1>(6, false, 5.5));
            nodes.push_back(new Node<1>(7, false, 11.0));

            NodesOnlyMesh<1> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);

            unsigned old_num_rows = mesh.mpBoxCollection->CalculateNumberOfNodesInEachStrip().size();

            mesh.AddNodesToBoxes();
            mesh.LoadBalanceMesh();

            unsigned new_num_rows = mesh.mpBoxCollection->CalculateNumberOfNodesInEachStrip().size();

            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(new_num_rows, old_num_rows - 1);
            }
            else
            {
                TS_ASSERT_EQUALS(new_num_rows, old_num_rows + 1);
            }

            // Tidy up
            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }
        }
    }
};

#endif /*TESTNODESONLYMESH_HPP_*/
