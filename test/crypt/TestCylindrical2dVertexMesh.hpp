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
#ifndef TESTCYLINDRICAL2DVERTEXMESH_HPP_
#define TESTCYLINDRICAL2DVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "VertexMeshWriter.hpp"
#include "ArchiveOpener.hpp"

class TestCylindrical2dVertexMesh : public CxxTest::TestSuite
{
public:

    void TestEachNodeIsContainedInAtLeastOneElement()
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(18, 25, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            TS_ASSERT_LESS_THAN(0u, num_containing_elements);
        }
    }

    void TestMeshGetWidth()
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(4, 4);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Test CalculateBoundingBox() method
        ChasteCuboid<2> bounds=p_mesh->CalculateBoundingBox();
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], 3.5,              1e-4);// \todo this should really be 4 as mesh is periodic
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], 13.0*0.5/sqrt(3), 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], 0.0,    1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], 0.0,    1e-4);

        // Test GetWidth() method
        double width = p_mesh->GetWidth(0);
        double height = p_mesh->GetWidth(1);

        TS_ASSERT_DELTA(width, 4, 1e-4);
        TS_ASSERT_DELTA(height, 13.0*0.5/sqrt(3), 1e-4);
    }

    void TestGetVectorFromAtoB() throw (Exception)
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(4, 4);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        c_vector<double, 2> node18_location = p_mesh->GetNode(18)->rGetLocation();
        c_vector<double, 2> node19_location = p_mesh->GetNode(19)->rGetLocation();

        // Test a normal vector and distance calculation
        c_vector<double, 2> vector = p_mesh->GetVectorFromAtoB(node18_location, node19_location);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetDistanceBetweenNodes(18, 19), 1.0, 1e-4);

        // Test the opposite vector
        vector = p_mesh->GetVectorFromAtoB(node19_location, node18_location);
        TS_ASSERT_DELTA(vector[0], -1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);

        // Test a periodic calculation
        c_vector<double, 2> node16_location = p_mesh->GetNode(16)->rGetLocation();
        vector = p_mesh->GetVectorFromAtoB(node16_location, node19_location);

        TS_ASSERT_DELTA(vector[0], -1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);
    }

    void TestSetNodeLocationForCylindricalMesh() throw (Exception)
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(4, 4);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = -0.01;
        new_point_location[1] = 3.0*0.5/sqrt(3);
        ChastePoint<2> new_point(new_point_location);

        // This node was on left and is now near the right
        p_mesh->SetNode(12, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(12u)->rGetLocation()[0], 3.99, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(12u)->rGetLocation()[1], 3.0*0.5/sqrt(3), 1e-4);

        // This node has stayed close to where it was
        new_point.SetCoordinate(0, 0.2);
        p_mesh->SetNode(0, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 0.2, 1e-4);

        // This node was on right and is now near the left
        new_point.SetCoordinate(0, 4.1);
        p_mesh->SetNode(8, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(8u)->rGetLocation()[0], 0.1, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(8u)->rGetLocation()[1], 3.0*0.5/sqrt(3), 1e-4);
    }

    void TestAddNodeAndReMesh() throw (Exception)
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(6, 6);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 84u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 36u);

        // Choose a node on the left boundary
        ChastePoint<2> point = p_mesh->GetNode(18)->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(point[1], 4.0*0.5/sqrt(3), 1e-4);

        // Create a new node close to this node
        point.SetCoordinate(0, -0.01);
        point.SetCoordinate(1, 4.5);
        Node<2>* p_node = new Node<2>(p_mesh->GetNumNodes(), point);

        unsigned old_num_nodes = p_mesh->GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = p_mesh->AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Remesh to update correspondences
        VertexElementMap map(p_mesh->GetNumElements());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 85u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 36u);

        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 5.99, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], 4.5000, 1e-4);

        // Now tet AddNode() when mDeletedNodeIndices is populated

        // Label node 29 as deleted
        p_mesh->mDeletedNodeIndices.push_back(29);

        // Create a new node close to this node
        ChastePoint<2> point2;
        point2.SetCoordinate(0, 2.0);
        point2.SetCoordinate(1, 2.1);
        Node<2>* p_node2 = new Node<2>(p_mesh->GetNumNodes(), point);

        // Add this new node to the mesh
        new_index = p_mesh->AddNode(p_node2);
        TS_ASSERT_EQUALS(new_index, 29u);
    }

    void TestElementAreaPerimeterCentroidAndMoments()
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(4, 4);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 40u);

        // Test area and perimeter calculations for all elements
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(i), 0.8660, 1e-4);
            TS_ASSERT_DELTA(p_mesh->GetSurfaceAreaOfElement(i), 3.4641, 1e-4);
        }

        // Test centroid calculations for non-periodic element
        c_vector<double, 2> centroid = p_mesh->GetCentroidOfElement(5);
        TS_ASSERT_DELTA(centroid(0), 2.0, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 5.0*0.5/sqrt(3), 1e-4);

        // Test centroid calculations for periodic element
        centroid = p_mesh->GetCentroidOfElement(7);
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 5.0*0.5/sqrt(3), 1e-4);

        // Test CalculateMomentOfElement() for all elements
        // all elements are regular hexagons with edge 1/sqrt(3)
        c_vector<double, 3> moments;
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            moments = p_mesh->CalculateMomentsOfElement(i);

            TS_ASSERT_DELTA(moments(0), 5*sqrt(3)/16/9, 1e-6);    // Ixx
            TS_ASSERT_DELTA(moments(1), 5*sqrt(3)/16/9, 1e-6);    // Iyy
            TS_ASSERT_DELTA(moments(2), 0.0, 1e-6);    // Ixy
        }
    }

    void TestDivideElementAlongGivenAxis()
    {
        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(4, 4);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 40u);

        c_vector<double, 2> axis_of_division;
        axis_of_division(0) = 1.0/sqrt(2.0);
        axis_of_division(1) = 1.0/sqrt(2.0);

        //Divide non periodic element
        unsigned new_element_index = p_mesh->DivideElementAlongGivenAxis(p_mesh->GetElement(2), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 17u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 42u);

        TS_ASSERT_DELTA(p_mesh->GetNode(40)->rGetLocation()[0], 2.8660, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(40)->rGetLocation()[1], 0.9433, 1e-4);

        TS_ASSERT_DELTA(p_mesh->GetNode(41)->rGetLocation()[0], 2.1339, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(41)->rGetLocation()[1], 0.2113, 1e-4);

        // Test new elements have correct nodes
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNode(1)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNode(2)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNode(3)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNode(4)->GetIndex(), 41u);

        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNode(0)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNode(1)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNode(2)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNode(3)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNode(4)->GetIndex(), 41u);

        // Divide periodic element
        new_element_index = p_mesh->DivideElementAlongGivenAxis(p_mesh->GetElement(3), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, 17u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 18u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 44u);

        TS_ASSERT_DELTA(p_mesh->GetNode(42)->rGetLocation()[0], -0.1339, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(42)->rGetLocation()[1], 0.9433, 1e-4);

        TS_ASSERT_DELTA(p_mesh->GetNode(43)->rGetLocation()[0], 3.1339, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(43)->rGetLocation()[1], 0.2113, 1e-4);

        // Test new elements have correct nodes
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNode(0)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNode(2)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNode(3)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNode(4)->GetIndex(), 43u);

        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNode(0)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNode(1)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNode(2)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNode(3)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNode(4)->GetIndex(), 43u);
    }

    void TestArchiving() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "cylindrical_vertex_mesh_base.arch";
        ArchiveLocationInfo::SetMeshFilename("cylindrical_vertex_mesh");

        // Create mesh
        unsigned num_cells_across = 4;
        unsigned num_cells_up = 7;
        CylindricalHoneycombVertexMeshGenerator generator(num_cells_across, num_cells_up);
        AbstractMesh<2,2>* const p_saved_mesh = generator.GetCylindricalMesh();

        double crypt_width = num_cells_across;

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
            TS_ASSERT_DELTA((static_cast<Cylindrical2dVertexMesh*>(p_saved_mesh))->GetWidth(0), crypt_width, 1e-7);

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
            Cylindrical2dVertexMesh* p_mesh2 = static_cast<Cylindrical2dVertexMesh*>(p_loaded_mesh);
            Cylindrical2dVertexMesh* p_mesh = static_cast<Cylindrical2dVertexMesh*>(p_saved_mesh);

            // Compare width
            TS_ASSERT_DELTA(p_mesh2->GetWidth(0), crypt_width, 1e-7);
            TS_ASSERT_DELTA(p_mesh->GetWidth(0), crypt_width, 1e-7);

            // Compare nodes
            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), p_mesh2->GetNumNodes());

            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                Node<2>* p_node = p_mesh->GetNode(i);
                Node<2>* p_node2 = p_mesh2->GetNode(i);
                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-4);
                }
            }

            // Compare elements
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), p_mesh2->GetNumElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh2->GetNumAllElements());

            for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
            {
                VertexElement<2,2>* p_elt = p_mesh->GetElement(i);
                VertexElement<2,2>* p_elt2 = p_mesh2->GetElement(i);
                TS_ASSERT_EQUALS(p_elt->GetNumNodes(), p_elt2->GetNumNodes());
                for (unsigned i=0; i<p_elt->GetNumNodes(); i++)
                {
                    TS_ASSERT_EQUALS(p_elt->GetNodeGlobalIndex(i), p_elt2->GetNodeGlobalIndex(i));
                }
            }

            // Tidy up
            delete p_mesh2;
        }
    }

    void TestCylindricalReMesh() throw (Exception)
    {
        // Create mesh
        unsigned num_cells_across = 6;
        unsigned num_cells_up = 12;
        CylindricalHoneycombVertexMeshGenerator generator(num_cells_across, num_cells_up);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Remesh
        VertexElementMap map(p_mesh->GetNumElements());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that there are the correct number of elements
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_cells_across*num_cells_up);
    }

    void TestCylindricalReMeshAfterDelete() throw (Exception)
    {
        // Create mesh
        unsigned num_cells_across = 6;
        unsigned num_cells_up = 12;
        CylindricalHoneycombVertexMeshGenerator generator(num_cells_across, num_cells_up);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        unsigned num_old_nodes = p_mesh->GetNumNodes();
        unsigned num_old_elements = num_cells_across*num_cells_up;

        // Delete a node
        p_mesh->DeleteElementPriorToReMesh(8);

        // Remesh
        VertexElementMap map(p_mesh->GetNumElements());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);
        TS_ASSERT_EQUALS(map.Size(), num_old_elements);

        // Check that there are the correct number of elements and nodes
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), num_old_nodes);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_old_elements-1);
    }

    void TestCylindricalElementIncludesPointAndGetLocalIndexForElementEdgeClosestToPoint()
    {
        // Set up a simple cylindrical mesh with one triangular element

        // Make 3 nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 9.0, 2.0));
        nodes.push_back(new Node<2>(1, false, 9.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        // Make element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Make mesh
        Cylindrical2dVertexMesh mesh(10.0, nodes, elements);

        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 4.0, 1e-10);

        // Make some test points and test ElementIncludesPoint()

        // A point far outside the element
        c_vector<double, 2> test_point1;
        test_point1[0] = -1.0;
        test_point1[1] = -1.0;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point1, 0), false);

        // A point outside the element due to periodicity
        c_vector<double, 2> test_point2;
        test_point1[0] = 3.0;
        test_point1[1] = 1.0;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point1, 0), false);


        // A point far inside the element
        c_vector<double, 2> test_point3;
        test_point3[0] = 9.5;
        test_point3[1] = 1.0;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point3, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point3, 0), 0u);

        // A point far inside the element as periodic
        c_vector<double, 2> test_point4;
        test_point4[0] = 0.5;
        test_point4[1] = 1.0;

        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point4, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point4, 0), 2u);
    }
};

#endif /*TESTCYLINDRICAL2DVERTEXMESH_HPP_*/
