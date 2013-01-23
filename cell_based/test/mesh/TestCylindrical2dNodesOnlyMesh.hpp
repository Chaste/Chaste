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

#ifndef TESTCYLINDRICAL2DNODESONLYMESH_HPP_
#define TESTCYLINDRICAL2DNODESONLYMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "HoneycombMeshGenerator.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "ArchiveOpener.hpp"

class TestCylindrical2dNodesOnlyMesh : public CxxTest::TestSuite
{
public:

    void TestMeshGetWidth()
    {
        // Create generating mesh
        HoneycombMeshGenerator generator(4, 4);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        double periodic_width = 4.0;
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Test CalculateBoundingBox() method
        ChasteCuboid<2> bounds = p_mesh->CalculateBoundingBox();

        ///\todo this should really be 4 as mesh is periodic
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], 3.5, 1e-4);

        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], 3.0*0.5*sqrt(3), 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], 0.0,1e-4);

        // Test GetWidth() method
        double width = p_mesh->GetWidth(0);
        double height = p_mesh->GetWidth(1);

        TS_ASSERT_DELTA(width, 4, 1e-4);
        TS_ASSERT_DELTA(height, 3.0*0.5*sqrt(3), 1e-4);

        // Avoid memory leak
        delete p_mesh;
    }

    void TestGetVectorFromAtoB() throw (Exception)
    {
        // Create generating mesh
        HoneycombMeshGenerator generator(4, 4);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        double periodic_width = 4.0;
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        c_vector<double, 2> node10_location = p_mesh->GetNode(10)->rGetLocation();
        c_vector<double, 2> node11_location = p_mesh->GetNode(11)->rGetLocation();

        // Test a normal vector and distance calculation
        c_vector<double, 2> vector = p_mesh->GetVectorFromAtoB(node10_location, node11_location);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetDistanceBetweenNodes(10, 11), 1.0, 1e-4);

        // Test the opposite vector
        vector = p_mesh->GetVectorFromAtoB(node11_location, node10_location);
        TS_ASSERT_DELTA(vector[0], -1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);

        // Test a periodic calculation
        c_vector<double, 2> node12_location = p_mesh->GetNode(12)->rGetLocation();
        vector = p_mesh->GetVectorFromAtoB(node11_location, node12_location);

        TS_ASSERT_DELTA(vector[0], 1.5, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.5*sqrt(3.0), 1e-4);

        // Test the opposite vector
        vector = p_mesh->GetVectorFromAtoB(node12_location, node11_location);
        TS_ASSERT_DELTA(vector[0], -1.5, 1e-4);
        TS_ASSERT_DELTA(vector[1], -0.5*sqrt(3.0), 1e-4);

        // Avoid memory leak
        delete p_mesh;
    }

    void TestSetNodeLocationForCylindricalMesh() throw (Exception)
    {
        // Create generating mesh
        HoneycombMeshGenerator generator(4, 4);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        double periodic_width = 4.0;
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Move one of the nodes to past the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = -0.01;
        new_point_location[1] = 0.5*sqrt(3);
        ChastePoint<2> new_point(new_point_location);

        // This node was on left and is now near the right
        p_mesh->SetNode(4, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(4u)->rGetLocation()[0], 3.99, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(4u)->rGetLocation()[1], 3.0*0.5/sqrt(3), 1e-4);

        // This node has stayed close to where it was
        TS_ASSERT_DELTA(p_mesh->GetNode(5u)->rGetLocation()[0], 1.5, 1e-4);
        new_point.SetCoordinate(0, 1.4);
        p_mesh->SetNode(5, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(5u)->rGetLocation()[0], 1.4, 1e-4);

        // This node was on right and is now near the left
        new_point.SetCoordinate(0, 4.1);
        p_mesh->SetNode(7, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(7u)->rGetLocation()[0], 0.1, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(7u)->rGetLocation()[1], 3.0*0.5/sqrt(3), 1e-4);

        // Avoid memory leak
        delete p_mesh;
    }

    void TestAddNode() throw (Exception)
    {
        // Create generating mesh
        HoneycombMeshGenerator generator(4, 4);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        double periodic_width = 4.0;
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Choose a node on the left boundary
        ChastePoint<2> point = p_mesh->GetNode(4)->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(point[1], 0.5*sqrt(3), 1e-4);

        // Create a new node close to this node
        point.SetCoordinate(0, -0.01);
        Node<2>* p_node = new Node<2>(p_mesh->GetNumNodes(), point);

        unsigned old_num_nodes = p_mesh->GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = p_mesh->AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 17u);

        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 3.99, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], 0.5*sqrt(3), 1e-4);

        // Now test AddNode() when mDeletedNodeIndices is populated

        // Label node 14 as deleted
        p_mesh->mDeletedNodeIndices.push_back(14);

        // Create a new node
        ChastePoint<2> point2;
        point2.SetCoordinate(0, 2.0);
        point2.SetCoordinate(1, 2.1);
        Node<2>* p_node2 = new Node<2>(p_mesh->GetNumNodes(), point);

        // Add this new node to the mesh
        new_index = p_mesh->AddNode(p_node2);
        TS_ASSERT_EQUALS(new_index, 14u);

        // Avoid memory leak
        delete p_mesh;
    }

    // NB This checks that periodicity is maintained through archiving...
    void TestArchiving() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "cylindrical_nodes_only_mesh_base.arch";
        ArchiveLocationInfo::SetMeshFilename("cylindrical_nodes_only_mesh");

        // Create generating mesh
        unsigned num_cells_across = 4;
        unsigned num_cells_up = 7;
        HoneycombMeshGenerator generator(num_cells_across,num_cells_up);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        double periodic_width = 4.0;
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        AbstractMesh<2,2>* const p_saved_mesh = p_mesh;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save
         * tracked objects while the compiler considers them const, to prevent
         * the objects changing during the save, and so object tracking leading
         * to wrong results. For example, A is saved once via pointer, then
         * changed, then saved again.  The second save notes that A was saved
         * before, so doesn't write its data again, and the change is lost.
         */
        {
            // Serialize the mesh
            TS_ASSERT_DELTA((static_cast<Cylindrical2dNodesOnlyMesh*>(p_saved_mesh))->GetWidth(0), periodic_width, 1e-7);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost.
            (*p_arch) << p_saved_mesh;
        }

        {
            // De-serialize and compare
            AbstractMesh<2,2>* p_loaded_mesh;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_loaded_mesh;

            // Compare the loaded mesh against the original
            Cylindrical2dNodesOnlyMesh* p_static_cast_loaded_mesh = static_cast<Cylindrical2dNodesOnlyMesh*>(p_loaded_mesh);
            Cylindrical2dNodesOnlyMesh* p_static_cast_saved_mesh = static_cast<Cylindrical2dNodesOnlyMesh*>(p_saved_mesh);

            // Compare width
            TS_ASSERT_DELTA(p_static_cast_loaded_mesh->GetWidth(0), periodic_width, 1e-7);
            TS_ASSERT_DELTA(p_static_cast_saved_mesh->GetWidth(0), periodic_width, 1e-7);

            // Compare nodes
            TS_ASSERT_EQUALS(p_static_cast_saved_mesh->GetNumNodes(), p_static_cast_loaded_mesh->GetNumNodes());
            TS_ASSERT_EQUALS(p_static_cast_saved_mesh->GetNumNodes(), 28u);

            for (unsigned i=0; i<p_static_cast_saved_mesh->GetNumNodes(); i++)
            {
                Node<2>* p_node = p_static_cast_saved_mesh->GetNode(i);
                Node<2>* p_node2 = p_static_cast_loaded_mesh->GetNode(i);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-4);
                }
            }
            // Avoid memory leak
            delete p_loaded_mesh;
        }

        // Avoid memory leak
        delete p_mesh;
    }
};

#endif /*TESTCYLINDRICAL2DNODESONLYMESH_HPP_*/
