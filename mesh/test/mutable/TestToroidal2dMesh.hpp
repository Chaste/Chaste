/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTTOROIDAL2DMESH_HPP_
#define TESTTOROIDAL2DMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include <algorithm>

#include "UblasCustomFunctions.hpp"
#include "VertexMesh.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"
#include "ArchiveOpener.hpp"
#include "TrianglesMeshWriter.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestToroidal2dMesh : public CxxTest::TestSuite
{
public:



    void TestCreateMirrorCellsAndAlignmentTester()
    {
        // Note that elements are not created (and boundary elements are not changed),
        // this just creates a set of new nodes

        unsigned cells_across = 4;
        unsigned cells_up = 4;
        double domain_width = (double) cells_across;
        double domain_height = (double) cells_up*0.5*sqrt(3.0);


        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1,1);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        // Reset the mesh
        p_mesh = generator.GetToroidalMesh();

        p_mesh->CreateMirrorNodes();

        // Check the vectors are the right size...
        TS_ASSERT_EQUALS(p_mesh->mLeftOriginals.size(), 3u*cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->mRightOriginals.size(), 3u*cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->mLeftOriginals.size(), p_mesh->mLeftImages.size());
        TS_ASSERT_EQUALS(p_mesh->mRightOriginals.size(), p_mesh->mRightImages.size());

        TS_ASSERT_EQUALS(p_mesh->mTopOriginals.size(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->mBottomOriginals.size(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->mTopOriginals.size(), p_mesh->mTopImages.size());
        TS_ASSERT_EQUALS(p_mesh->mBottomOriginals.size(), p_mesh->mBottomImages.size());

        // Check that the image nodes are where they should be.
        for (unsigned i=0; i<p_mesh->mLeftOriginals.size(); i++)
        {
            c_vector<double,2> original_location;
            original_location = p_mesh->GetNode(p_mesh->mLeftOriginals[i])->rGetLocation();
            c_vector<double,2> image_location;
            image_location = p_mesh->GetNode(p_mesh->mLeftImages[i])->rGetLocation();

            TS_ASSERT_DELTA(original_location[0] + domain_width, image_location[0], 1e-7);
            TS_ASSERT_DELTA(original_location[1], image_location[1], 1e-7);
        }

        for (unsigned i=0; i<p_mesh->mRightOriginals.size(); i++)
        {
            c_vector<double,2> original_location;
            original_location = p_mesh->GetNode(p_mesh->mRightOriginals[i])->rGetLocation();
            c_vector<double,2> image_location;
            image_location = p_mesh->GetNode(p_mesh->mRightImages[i])->rGetLocation();

            TS_ASSERT_DELTA(original_location[0] - domain_width, image_location[0], 1e-7);
            TS_ASSERT_DELTA(original_location[1], image_location[1], 1e-7);
        }

        for (unsigned i=0; i<p_mesh->mBottomOriginals.size(); i++)
        {
            c_vector<double,2> original_location;
            original_location = p_mesh->GetNode(p_mesh->mBottomOriginals[i])->rGetLocation();
            c_vector<double,2> image_location;
            image_location = p_mesh->GetNode(p_mesh->mBottomImages[i])->rGetLocation();

            TS_ASSERT_DELTA(original_location[1] + domain_height, image_location[1], 1e-7);
            TS_ASSERT_DELTA(original_location[0], image_location[0], 1e-7);
        }

        for (unsigned i=0; i<p_mesh->mTopOriginals.size(); i++)
        {
            c_vector<double,2> original_location;
            original_location = p_mesh->GetNode(p_mesh->mTopOriginals[i])->rGetLocation();
            c_vector<double,2> image_location;
            image_location = p_mesh->GetNode(p_mesh->mTopImages[i])->rGetLocation();

            TS_ASSERT_DELTA(original_location[1] - domain_height, image_location[1], 1e-7);
            TS_ASSERT_DELTA(original_location[0], image_location[0], 1e-7);
        }


        // Check that we've got the correct number
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 9*cells_across*cells_up);
    }

    // void TestReconstructToroidalMesh()
    // {
    //     // This test takes in a new mesh created using the mirror function above
    //     // and a ReMesh call, then removes nodes, elements and boundary elements

    //     unsigned cells_across = 6;
    //     unsigned cells_up = 12;

    //     // Set up a mesh which can be mirrored (no scaling in this case)
    //     ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1 , 1);
    //     boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

    //     // Create a mirrored load of nodes for the normal remesher to work with
    //     p_mesh->CreateMirrorNodes();

    //     // Call the normal re-mesh
    //     NodeMap map(p_mesh->GetNumNodes());
    //     p_mesh->MutableMesh<2,2>::ReMesh(map);

    //     // Re-Index the vectors regarding left/right nodes with the node map
    //     for (unsigned i = 0; i<p_mesh->mLeftOriginals.size(); i++)
    //     {
    //          p_mesh->mLeftOriginals[i] = map.GetNewIndex(p_mesh->mLeftOriginals[i]);
    //          p_mesh->mLeftImages[i] = map.GetNewIndex(p_mesh->mLeftImages[i]);
    //     }
    //     for (unsigned i = 0; i<p_mesh->mRightOriginals.size(); i++)
    //     {
    //          p_mesh->mRightOriginals[i] = map.GetNewIndex(p_mesh->mRightOriginals[i]);
    //          p_mesh->mRightImages[i] = map.GetNewIndex(p_mesh->mRightImages[i]);
    //     }

    //     p_mesh->ReconstructToroidalMesh();

    //     unsigned elements_for_node_0 = 0;
    //     unsigned elements_for_node_11 = 0;
    //     unsigned elements_for_node_12 = 0;
    //     unsigned elements_for_node_18 = 0;

    //     unsigned checksum_for_node_0 = 0;
    //     unsigned checksum_for_node_11 = 0;
    //     unsigned checksum_for_node_12 = 0;
    //     unsigned checksum_for_node_18 = 0;

    //     for (Toroidal2dMesh::ElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
    //          elem_iter != p_mesh->GetElementIteratorEnd();
    //          ++elem_iter)
    //     {
    //         for (unsigned i=0; i<3; i++)
    //         {
    //             unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);

    //             if (this_node_index==0)
    //             {
    //                 elements_for_node_0++;
    //                 for (unsigned j=0; j<3; j++)
    //                 {
    //                     checksum_for_node_0 += elem_iter->GetNodeGlobalIndex(j);
    //                 }
    //             }
    //             if (this_node_index==11)
    //             {
    //                 elements_for_node_11++;
    //                 for (unsigned j=0; j<3; j++)
    //                 {
    //                     checksum_for_node_11 += elem_iter->GetNodeGlobalIndex(j);
    //                 }
    //             }
    //             if (this_node_index==12)
    //             {
    //                 elements_for_node_12++;
    //                 for (unsigned j=0; j<3; j++)
    //                 {
    //                     checksum_for_node_12 += elem_iter->GetNodeGlobalIndex(j);
    //                 }
    //             }
    //             if (this_node_index==18)
    //             {
    //                 elements_for_node_18++;
    //                 for (unsigned j=0; j<3; j++)
    //                 {
    //                     checksum_for_node_18 += elem_iter->GetNodeGlobalIndex(j);
    //                 }
    //             }
    //         }
    //     }

    //     TS_ASSERT_EQUALS(elements_for_node_0, 3u);
    //     TS_ASSERT_EQUALS(elements_for_node_11, 6u);
    //     TS_ASSERT_EQUALS(elements_for_node_12, 6u);
    //     TS_ASSERT_EQUALS(elements_for_node_18, 6u);

    //     // These are nodes on the edge of the cylindrical region.
    //     // If the mesh is periodic they will be joined to the following other nodes.
    //     unsigned checksum_target_for_node_0 = (0 + 5 + 11) + (0 + 6 + 11) + (0 + 1 + 6);
    //     TS_ASSERT_EQUALS(checksum_for_node_0, checksum_target_for_node_0);

    //     unsigned checksum_target_for_node_11 = (0 + 5 + 11) + (0 + 6 + 11) + (11 + 12 + 6)
    //                                           +(11 + 17 + 12) + (10 + 11 + 17) + (5 + 10 + 11);
    //     TS_ASSERT_EQUALS(checksum_for_node_11, checksum_target_for_node_11);

    //     unsigned checksum_target_for_node_12 = (12 + 13 + 18) + (12 + 18 + 23) + (12 + 17 + 23)
    //                                           +(12 + 11 + 17) + (12 + 6 + 11) + (12 + 6 + 13);
    //     TS_ASSERT_EQUALS(checksum_for_node_12, checksum_target_for_node_12);

    //     unsigned checksum_target_for_node_18 = (18 + 19 + 25) + (18 + 24 + 25) + (18 + 23 + 24)
    //                                             +(18 + 12 + 23) + (18 + 12 + 13) + (18 + 13 + 19);
    //     TS_ASSERT_EQUALS(checksum_for_node_18, checksum_target_for_node_18);

    //     // Check that there are the correct number of everything
    //     TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
    //     TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*(cells_up-1));
    //     TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 12u);
    // }

    void TestToroidalReMesh()
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;

        // Set up a mesh which can be mirrored (no scaling in this case)
        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1 , 1);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u);  // No boundary elements so one added to avoid errors.
    }

    void TestToroidalReMeshAfterDelete()
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;

        // Set up a mesh which can be mirrored (no scaling in this case)
        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1 , 1);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        unsigned num_old_nodes = p_mesh->GetNumNodes();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up)

        p_mesh->DeleteNode(15);

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up-1)

        NodeMap map(p_mesh->GetNumNodes());

        p_mesh->ReMesh(map);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up-1);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u);  // No boundary elements so one addd to avoid errors

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

    void TestToroidalReMeshOnSmallMesh()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;

        // Set up a mesh which can be mirrored (no scaling in this case)
        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1 , 1);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        NodeMap map(p_mesh->GetNumNodes());
        p_mesh->ReMesh(map);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u); // boundary elements removed so one add to avoid errors
    }

    void TestGetVectorBetweenPeriodicPoints()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;

        // Set up a mesh which can be mirrored (no scaling in this case)
        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1 , 1);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        c_vector<double, 2> location1 = p_mesh->GetNode(1)->rGetLocation();
        c_vector<double, 2> location2 = p_mesh->GetNode(4)->rGetLocation();

        // Test a normal distance calculation...
        c_vector<double, 2> vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], sqrt(3.0)/2.0, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetDistanceBetweenNodes(1, 4), 1.0, 1e-7);

        // Periodic at x=3, so closest node should be node 0 at (0,0).
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
        location1[1] = 2.5*sqrt(3)/2.0;
        location2[0] = 0.5;
        location2[1] = 0.5*sqrt(3)/2.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], 1.0*sqrt(3)/2.0, 1e-7);

        // We also want GetVectorFromAtoB to work when the x coord of A and B is not
        // between 0 and domain width, by first normalizing the x coord (by taking modulus by crypt width)
        location1[0] = -0.5;
        location1[1] = -0.5*sqrt(3)/2.0;
        location2[0] = 2.0;
        location2[1] = 2.0*sqrt(3)/2.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], -0.5, 1e-7);
        TS_ASSERT_DELTA(vector[1], -0.5*sqrt(3)/2.0, 1e-7);

        // Some tests where the location[0] is not between -0.25*crypt_width and 1.25*crypt_width
        location1[0] = -2.5;
        location1[1] = 0;
        location2[0] = 4.5;
        location2[1] = 1.0;
        vector = p_mesh->GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], 1.0, 1e-7);
    }

    void TestSetNodeLocationForToroidalMesh()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;

        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1, 1);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        // Move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = 2.999999999;
        new_point_location[1] = 0.0001;
        ChastePoint<2> new_point(new_point_location);
        p_mesh->GetNode(5)->SetPoint(new_point);

        new_point.SetCoordinate(0, -0.0001);
        new_point.SetCoordinate(1, -0.0001);

        // This node was on bottom right and is now near the top right
        p_mesh->SetNode(0, new_point, false);
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 2.9999, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[1], 1.5*sqrt(3)-0.0001, 1e-4);

        new_point.SetCoordinate(0, 3.0001);
        p_mesh->SetNode(0, new_point, false);

        // This node was on right and is now on the left
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 0.0001, 1e-4);
    }

    void TestAddNodeAndReMesh()
    {
        unsigned cells_across = 3;
        unsigned cells_up = 3;

        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1, 1);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*cells_up);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u); // boundary elements removed so one add to avoid errors

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
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2*cells_across*cells_up+2);
        TS_ASSERT_EQUALS(p_mesh->GetNumBoundaryElements(), 1u); // boundary elements removed so one add to avoid errors+

        // Check that we have moved the new node across
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 3.0+point[0], 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], point[1], 1e-7);

        // Test GetWidth
        TS_ASSERT_DELTA(p_mesh->GetWidth(0u), 3.0, 1e-9);
        TS_ASSERT_DELTA(p_mesh->GetWidth(1u), 3.0*0.5*sqrt(3.0), 1e-6);
    }

    // NB This checks that periodicity is maintained through archiving...
    void TestArchiving()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "toriodal_mesh_base.arch";
        ArchiveLocationInfo::SetMeshFilename("toroidal_mesh");

        // Set up a mesh
        unsigned cells_across = 5;
        unsigned cells_up = 3;

        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1, 1);
        boost::shared_ptr<AbstractTetrahedralMesh<2,2> > const p_mesh = boost::static_pointer_cast<AbstractTetrahedralMesh<2,2> >(generator.GetToroidalMesh());

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
            TS_ASSERT_DELTA(width, 5.0, 1e-7);
            double height = p_mesh->GetWidth(1);
            TS_ASSERT_DELTA(height, 1.5*sqrt(3.0), 1e-7);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost.
            (*p_arch) << p_mesh;
        }

        {
            // De-serialize and compare
            boost::shared_ptr<AbstractTetrahedralMesh<2,2> > p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            // Toroidal2dMesh now remeshes itself on load (convert from TrianglesMeshReader to normal format)

            TS_ASSERT_DELTA(p_mesh2->GetWidth(0), 5.0, 1e-7);
            TS_ASSERT_DELTA(p_mesh2->GetWidth(1), 1.5*sqrt(3.0), 1e-7);

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
        }
    }

    void TestConstructFromNodeList()
    {
        std::vector<Node<2>*> nodes;

        nodes.push_back(new Node<2>(0, true, 0.1, 0.01));
        nodes.push_back(new Node<2>(1, true, 0.5, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.9, 0.01));
        nodes.push_back(new Node<2>(3, true, 0.1, 0.99));
        nodes.push_back(new Node<2>(4, true, 0.5, 0.9));
        nodes.push_back(new Node<2>(5, true, 0.9, 0.99));

        const double width = 1.0;
        const double height = 1.0;

        Toroidal2dMesh mesh(width, height, nodes);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 12u);

        // Find the element with node indices 2,0,3 (which stradles both periodic boundaries)
        unsigned element_index;
        std::set<unsigned> target_element_node_indices;
        for (element_index=0; element_index<mesh.GetNumElements(); element_index++)
        {
            target_element_node_indices.clear();
            target_element_node_indices.insert(2);
            target_element_node_indices.insert(0);
            target_element_node_indices.insert(3);

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
        c_matrix<double, 2, 2> jacobian;
        c_matrix<double, 2, 2> inverse_jacobian;
        double jacobian_det;

        mesh.GetInverseJacobianForElement(element_index, jacobian, jacobian_det, inverse_jacobian);
        c_vector<double, 3> circumsphere = mesh.GetElement(element_index)->CalculateCircumsphere(jacobian,inverse_jacobian);

        TS_ASSERT_DELTA(circumsphere[0], 1.0, 1e-3);
        TS_ASSERT_DELTA(circumsphere[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(circumsphere[2], 0.25, 1e-3);

        /* The reason that the circumsphere is calculated correctly for a periodic boundary
         * stradling element is somewhat obscure.
         * The Jacobian of the element is calculated when the element has mirror nodes
         * The mirror nodes are then replaced with the nodes within the periodic mesh
         * The circumsphere is calculated based on the Jacobian and the replaced node within the periodic mesh
         *
         * uncommenting the following line of code causes an error:
         * Jacobian determinant is non-positive
         *
         * mesh.GetElement(element_index)->RefreshJacobianDeterminant();
         */
    }

    void TestVoronoiTessellationUsesOverriddenMetric()
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;

        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up, 1, 1);
        boost::shared_ptr<Toroidal2dMesh> const p_mesh = generator.GetToroidalMesh();

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

    void TestRefreshMesh()
    {
        // Create a simple Toroidal2dMesh
        unsigned cells_across = 4;
        unsigned cells_up = 4;
        ToroidalHoneycombMeshGenerator generator(cells_across, cells_up);
        boost::shared_ptr<Toroidal2dMesh> p_mesh = generator.GetToroidalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 32u);

        //Translate mesh which calls RefreshMesh()
        p_mesh->Translate(-1,-1);

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 32u);

        //Translate mesh which calls RefreshMesh()
        p_mesh->Translate(2,2);

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 32u);
    }

};

#endif /*TESTTOROIDAL2DMESH_HPP_*/
