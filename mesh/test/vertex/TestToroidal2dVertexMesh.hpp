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

#ifndef TESTTOROIDAL2DVERTEXMESH_HPP_
#define TESTTOROIDAL2DVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexMeshReader.hpp"
#include "ArchiveOpener.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestToroidal2dVertexMesh : public CxxTest::TestSuite
{
public:

    void TestEachNodeIsContainedInAtLeastOneElement()
    {
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(18, 24);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 864u); // 2*18*24
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 432u); // 18*24

        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            TS_ASSERT_LESS_THAN(0u, num_containing_elements);
        }
    }

    void TestMeshGetAndSetWidthAndHeight()
    {
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(4, 4);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        // Test CalculateBoundingBox() method
        ChasteCuboid<2> bounds = p_mesh->CalculateBoundingBox();

        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[0], 3.5, 1e-4);
        TS_ASSERT_DELTA(bounds.rGetUpperCorner()[1], 10.0*0.5/sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounds.rGetLowerCorner()[1], 0.0, 1e-4);

        // Test GetWidth() method
        double width = p_mesh->GetWidth(0);
        double height = p_mesh->GetWidth(1);
        TS_ASSERT_DELTA(width, 4.0, 1e-4);
        TS_ASSERT_DELTA(height, 6.0/sqrt(3.0), 1e-4);

        // Test SetWidth() method
        p_mesh->SetWidth(61.47);
        width = p_mesh->GetWidth(0);
        TS_ASSERT_DELTA(width, 61.47, 1e-4);

        // Test SetHeight() method
        p_mesh->SetHeight(123.89);
        height = p_mesh->GetWidth(1);
        TS_ASSERT_DELTA(height, 123.89, 1e-4);
    }

    void TestGetVectorFromAtoB()
    {
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(4, 4);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        // Store the locations of some nodes
        c_vector<double, 2> node0_location = p_mesh->GetNode(0)->rGetLocation();
        c_vector<double, 2> node12_location = p_mesh->GetNode(12)->rGetLocation();
        c_vector<double, 2> node16_location = p_mesh->GetNode(16)->rGetLocation();
        c_vector<double, 2> node18_location = p_mesh->GetNode(18)->rGetLocation();
        c_vector<double, 2> node19_location = p_mesh->GetNode(19)->rGetLocation();
        c_vector<double, 2> node28_location = p_mesh->GetNode(28)->rGetLocation();
        c_vector<double, 2> node31_location = p_mesh->GetNode(31)->rGetLocation();

        // Test a normal vector and distance calculation
        c_vector<double, 2> vector = p_mesh->GetVectorFromAtoB(node18_location, node19_location);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1.0, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetDistanceBetweenNodes(18, 19), 1.0, 1e-4);

        // Test the opposite vector
        vector = p_mesh->GetVectorFromAtoB(node19_location, node18_location);
        TS_ASSERT_DELTA(vector[0], -1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0, 1e-4);

        // Test another normal calculation
        vector = p_mesh->GetVectorFromAtoB(node12_location, node16_location);
        TS_ASSERT_DELTA(vector[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 1/sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 1/sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetDistanceBetweenNodes(12, 16), 1/sqrt(3.0), 1e-4);

        vector = p_mesh->GetVectorFromAtoB(node16_location, node12_location);
        TS_ASSERT_DELTA(vector[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], -1/sqrt(3.0), 1e-4);

        // Test a left-right periodic calculation
        vector = p_mesh->GetVectorFromAtoB(node16_location, node19_location);
        TS_ASSERT_DELTA(vector[0], -1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1],  0.0, 1e-4);

        // Test a top-bottom periodic calculation
        vector = p_mesh->GetVectorFromAtoB(node28_location, node0_location);
        TS_ASSERT_DELTA(vector[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 1/sqrt(3.0), 1e-4);

        vector = p_mesh->GetVectorFromAtoB(node0_location, node28_location);
        TS_ASSERT_DELTA(vector[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], -1/sqrt(3.0), 1e-4);

        // Test a 'diagonal' periodic calculation
        vector = p_mesh->GetVectorFromAtoB(node31_location, node0_location);
        TS_ASSERT_DELTA(vector[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(vector[1], 1/sqrt(3.0), 1e-4);
    }

    void TestSetNodeLocationForToroidalMesh()
    {
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(4, 4);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        // Move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = -0.01;
        new_point_location[1] = 3.0*0.5/sqrt(3.0);
        ChastePoint<2> new_point(new_point_location);

        // This node was on left and is now near the right
        p_mesh->SetNode(12, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(12)->rGetLocation()[0], 3.99, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(12)->rGetLocation()[1], 3.0*0.5/sqrt(3.0), 1e-4);

        // This node has stayed close to where it was
        new_point.SetCoordinate(0, 0.2);
        p_mesh->SetNode(0, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetLocation()[0], 0.2, 1e-4);

        // This node was on the right and is now near the left
        new_point.SetCoordinate(0, 4.1);
        p_mesh->SetNode(8, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(8)->rGetLocation()[0], 0.1, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(8)->rGetLocation()[1], 3.0*0.5/sqrt(3.0), 1e-4);

        // This node was on the top and is now near the bottom
        new_point.SetCoordinate(0, 3.5);
        new_point.SetCoordinate(1, 6.0/sqrt(3.0)+0.1);
        p_mesh->SetNode(31, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(31)->rGetLocation()[0], 3.5, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(31)->rGetLocation()[1], 0.1, 1e-4);

        // This node was on the bottom and is now near the top
        new_point.SetCoordinate(0, 2.5);
        new_point.SetCoordinate(1, -0.1);
        p_mesh->SetNode(2, new_point);
        TS_ASSERT_DELTA(p_mesh->GetNode(2)->rGetLocation()[0], 2.5, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(2)->rGetLocation()[1], 6.0/sqrt(3.0) - 0.1, 1e-4);
    }

    void TestAddNodeAndReMesh()
    {
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(6, 6);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 72u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 36u);

        // Choose a node on the left boundary
        ChastePoint<2> point = p_mesh->GetNode(18)->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(point[1], 4.0*0.5/sqrt(3.0), 1e-4);

        // Create a new node close to this node
        point.SetCoordinate(0, -0.01);
        point.SetCoordinate(1, 4.0*0.5/sqrt(3.0));
        Node<2>* p_node = new Node<2>(p_mesh->GetNumNodes(), point);

        unsigned old_num_nodes = p_mesh->GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = p_mesh->AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Call ReMesh() to update correspondences
        VertexElementMap map(p_mesh->GetNumElements());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 73u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 36u);

        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 5.99, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], 4.0*0.5/sqrt(3.0), 1e-4);

        // Repeat this with a node on the top boundary
        point = p_mesh->GetNode(70)->GetPoint();
        TS_ASSERT_DELTA(point[0], 4.5, 1e-4);
        TS_ASSERT_DELTA(point[1], 8.0/sqrt(3.0), 1e-4);

        // Create a new node close to this node
        point.SetCoordinate(0, 4.5);
        point.SetCoordinate(1, 9.0/sqrt(3.0)+0.01);
        p_node = new Node<2>(p_mesh->GetNumNodes(), point);

        old_num_nodes = p_mesh->GetNumNodes();
        new_index = p_mesh->AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Call ReMesh() to update correspondences
        VertexElementMap map2(p_mesh->GetNumElements());
        p_mesh->ReMesh(map2);

        TS_ASSERT_EQUALS(map2.Size(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(map2.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 74u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 36u);

        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 4.5, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], 0.01, 1e-4);

        // Now test AddNode() when mDeletedNodeIndices is populated

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
        ToroidalHoneycombVertexMeshGenerator generator(4, 4);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u); // 2*4*4

        // Test area and perimeter calculations for all elements
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(i), 0.8660, 1e-4);
            TS_ASSERT_DELTA(p_mesh->GetSurfaceAreaOfElement(i), 3.4641, 1e-4);
        }

        // Test centroid calculations for non-periodic element
        c_vector<double, 2> centroid = p_mesh->GetCentroidOfElement(5);
        TS_ASSERT_DELTA(centroid(0), 2.0, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 5.0*0.5/sqrt(3.0), 1e-4);

        // Test centroid calculations for periodic element (left to right)
        centroid = p_mesh->GetCentroidOfElement(7);
        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 5.0*0.5/sqrt(3.0), 1e-4);

        // Test centroid calculations for periodic element (top to bottom)
        centroid = p_mesh->GetCentroidOfElement(12);
        TS_ASSERT_DELTA(centroid(0), 1.0, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 11.0*0.5/sqrt(3.0), 1e-4);

        // Test CalculateMomentOfElement() for all elements;
        // all elements are regular hexagons with edge length 1/sqrt(3.0)
        c_vector<double, 3> moments;
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            moments = p_mesh->CalculateMomentsOfElement(i);

            TS_ASSERT_DELTA(moments(0), 5*sqrt(3.0)/16/9, 1e-6); // Ixx
            TS_ASSERT_DELTA(moments(1), 5*sqrt(3.0)/16/9, 1e-6); // Iyy
            TS_ASSERT_DELTA(moments(2), 0.0, 1e-6);            // Ixy = 0 by symmetry
        }

        // Test methods with a toroidal mesh comprising a single rectangular element
        std::vector<Node<2>*> rectangle_nodes;
        rectangle_nodes.push_back(new Node<2>(0, false, 8.0, 2.0));
        rectangle_nodes.push_back(new Node<2>(1, false, 8.0, 0.0));
        rectangle_nodes.push_back(new Node<2>(2, false, 2.0, 0.0));
        rectangle_nodes.push_back(new Node<2>(3, false, 2.0, 2.0));
        std::vector<VertexElement<2,2>*> rectangle_elements;
        rectangle_elements.push_back(new VertexElement<2,2>(0, rectangle_nodes));

        Toroidal2dVertexMesh rectangle_mesh(10.0, 3.0, rectangle_nodes, rectangle_elements);

        TS_ASSERT_DELTA(rectangle_mesh.GetVolumeOfElement(0), 4.0, 1e-10);
        TS_ASSERT_DELTA(rectangle_mesh.GetSurfaceAreaOfElement(0), 10.0, 1e-4);

        ///\todo #2393 - for consistency, the centroid should be at (0, 2.5)
        c_vector<double, 2> rectangle_centroid = rectangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(rectangle_centroid(0), 10.0, 1e-4);
        TS_ASSERT_DELTA(rectangle_centroid(1), 2.5, 1e-4);

        c_vector<double, 3> rectangle_moments = rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(rectangle_moments(0), 1.0/3.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(rectangle_moments(1), 16.0/3.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(rectangle_moments(2), 0.0, 1e-6);      // Ixy = 0 by symmetry

        c_vector<double, 2> rectangle_short_axis = rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(rectangle_short_axis(0), 0.0, 1e-4);
        TS_ASSERT_DELTA(rectangle_short_axis(1), 1.0, 1e-4);
    }

    void TestDivideElementAlongGivenAxis()
    {
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(4, 4);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u);

        c_vector<double, 2> axis_of_division;
        axis_of_division(0) = 1.0/sqrt(2.0);
        axis_of_division(1) = 1.0/sqrt(2.0);

        // Divide non-periodic element
        unsigned new_element_index = p_mesh->DivideElementAlongGivenAxis(p_mesh->GetElement(2), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 17u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 34u);

        TS_ASSERT_DELTA(p_mesh->GetNode(32)->rGetLocation()[0], 2.8660, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(32)->rGetLocation()[1], 0.9433, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(33)->rGetLocation()[0], 2.1339, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(33)->rGetLocation()[1], 0.2113, 1e-4);

        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(1), 7u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(2), 11u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(3), 32u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(4), 33u);

        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNodeGlobalIndex(0), 32u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNodeGlobalIndex(1), 14u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNodeGlobalIndex(2), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNodeGlobalIndex(3), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(16)->GetNodeGlobalIndex(4), 33u);

        // Divide periodic element (left to right)
        new_element_index = p_mesh->DivideElementAlongGivenAxis(p_mesh->GetElement(3), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, 17u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 18u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 36u);

        TS_ASSERT_DELTA(p_mesh->GetNode(34)->rGetLocation()[0], -0.1339, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(34)->rGetLocation()[1], 0.9433, 1e-4);

        TS_ASSERT_DELTA(p_mesh->GetNode(35)->rGetLocation()[0], 3.1339, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(35)->rGetLocation()[1], 0.2113, 1e-4);

        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(2), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(3), 34u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(4), 35u);

        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNodeGlobalIndex(0), 34u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNodeGlobalIndex(1), 15u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNodeGlobalIndex(2), 11u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNodeGlobalIndex(3), 7u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(17)->GetNodeGlobalIndex(4), 35u);

        // Divide periodic element (top to bottom)
        new_element_index = p_mesh->DivideElementAlongGivenAxis(p_mesh->GetElement(13), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, 18u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 19u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 38u);

        TS_ASSERT_DELTA(p_mesh->GetNode(36)->rGetLocation()[0], 2.3660, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(36)->rGetLocation()[1], 0.0773, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(37)->rGetLocation()[0], 1.6339, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(37)->rGetLocation()[1], 2.8094, 1e-4);

        TS_ASSERT_EQUALS(p_mesh->GetElement(13)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(13)->GetNodeGlobalIndex(0), 36u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(13)->GetNodeGlobalIndex(1), 33u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(13)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(13)->GetNodeGlobalIndex(3), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(13)->GetNodeGlobalIndex(4), 29u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(13)->GetNodeGlobalIndex(5), 37u);

        TS_ASSERT_EQUALS(p_mesh->GetElement(18)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(18)->GetNodeGlobalIndex(0), 26u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(18)->GetNodeGlobalIndex(1), 30u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(18)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(18)->GetNodeGlobalIndex(3), 36u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(18)->GetNodeGlobalIndex(4), 37u);
    }

    void TestArchiving()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "toroidal_vertex_mesh_base.arch";
        ArchiveLocationInfo::SetMeshFilename("toroidal_vertex_mesh");

        // Create mesh
        unsigned num_cells_across = 4;
        unsigned num_cells_up = 6;
        ToroidalHoneycombVertexMeshGenerator generator(num_cells_across, num_cells_up);
        AbstractMesh<2,2>* const p_saved_mesh = generator.GetToroidalMesh();

        double mesh_width = num_cells_across;
        double mesh_height = num_cells_up*1.5/sqrt(3.0);

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
            TS_ASSERT_DELTA((static_cast<Toroidal2dVertexMesh*>(p_saved_mesh))->GetWidth(0), mesh_width, 1e-7);

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
            Toroidal2dVertexMesh* p_mesh2 = static_cast<Toroidal2dVertexMesh*>(p_loaded_mesh);
            Toroidal2dVertexMesh* p_mesh = static_cast<Toroidal2dVertexMesh*>(p_saved_mesh);

            // Compare width
            TS_ASSERT_DELTA(p_mesh2->GetWidth(0), mesh_width, 1e-7);
            TS_ASSERT_DELTA(p_mesh->GetWidth(0), mesh_width, 1e-7);
            TS_ASSERT_DELTA(p_mesh2->GetWidth(1), mesh_height, 1e-7);
            TS_ASSERT_DELTA(p_mesh->GetWidth(1), mesh_height, 1e-7);

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
                for (unsigned j=0; j<p_elt->GetNumNodes(); j++)
                {
                    TS_ASSERT_EQUALS(p_elt->GetNodeGlobalIndex(j), p_elt2->GetNodeGlobalIndex(j));
                }
            }

            // Tidy up
            delete p_mesh2;
        }
    }

    void TestToroidalReMesh()
    {
        // Create mesh
        unsigned num_cells_across = 6;
        unsigned num_cells_up = 12;
        ToroidalHoneycombVertexMeshGenerator generator(num_cells_across, num_cells_up);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        // Call Remesh()
        VertexElementMap map(p_mesh->GetNumElements());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that there are the correct number of elements
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_cells_across*num_cells_up);
    }

    void TestToroidalReMeshAfterDelete()
    {
        // Create mesh
        unsigned num_cells_across = 6;
        unsigned num_cells_up = 12;
        ToroidalHoneycombVertexMeshGenerator generator(num_cells_across, num_cells_up);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        unsigned num_old_nodes = p_mesh->GetNumNodes();
        unsigned num_old_elements = num_cells_across*num_cells_up;

        // Delete a node
        p_mesh->DeleteElementPriorToReMesh(8);

        // Call Remesh()
        VertexElementMap map(p_mesh->GetNumElements());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);
        TS_ASSERT_EQUALS(map.Size(), num_old_elements);

        // Check that there are the correct number of elements and nodes
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), num_old_nodes);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_old_elements-1);
    }

    void TestToroidalElementIncludesPointAndGetLocalIndexForElementEdgeClosestToPoint()
    {
        // Set up a simple toroidal mesh with one rectangular element

        // Make four nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 9.0, 2.0));
        nodes.push_back(new Node<2>(1, false, 9.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 2.0));

        // Make element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Make mesh
        Toroidal2dVertexMesh mesh(10.0, 3.0, nodes, elements);

        /*
         * The area of the element, taking periodicity into account,
         * should be equal to |((9-1) - 10)|*|((2-0) - 3)| = 2*1 = 2.
         */
        TS_ASSERT_DELTA(mesh.GetVolumeOfElement(0), 2.0, 1e-10);

        // Make some test points and test ElementIncludesPoint()

        // A point far outside the element
        c_vector<double, 2> test_point1;
        test_point1[0] = -1.1;
        test_point1[1] = -0.9;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point1, 0), false);

        // A point outside the element due to periodicity
        c_vector<double, 2> test_point2;
        test_point2[0] = 3.0;
        test_point2[1] = 1.0;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point2, 0), false);

        // Another point outside the element due to periodicity
        // (note that this point would be inside the element if we had a Cylindrical2dVertexMesh)
        c_vector<double, 2> test_point3;
        test_point3[0] = 9.5;
        test_point3[1] = 1.0;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point3, 0), false);

        // A point far inside the element as periodic
        c_vector<double, 2> test_point4;
        test_point4[0] = 0.9;
        test_point4[1] = -0.2;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point4, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point4, 0), 2u);

        // Another point far inside the element as periodic
        c_vector<double, 2> test_point5;
        test_point5[0] = 9.99;
        test_point5[1] = 2.1;
        TS_ASSERT_EQUALS(mesh.ElementIncludesPoint(test_point5, 0), true);
        TS_ASSERT_EQUALS(mesh.GetLocalIndexForElementEdgeClosestToPoint(test_point5, 0), 3u);
    }

    void TestGetMeshForVtk()
    {
        // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(4, 4);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u); // 2*4*4
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u); // 4*4

        // Test GetMeshForVtk() method
        VertexMesh<2, 2>* p_mesh_for_vtk = p_mesh->GetMeshForVtk();

        // The mesh for VTK should have the same number of elements, but 16 extra nodes
        TS_ASSERT_EQUALS(p_mesh_for_vtk->GetNumElements(), 16u);
        TS_ASSERT_EQUALS(p_mesh_for_vtk->GetNumNodes(), 48u);

        // Every element in the mesh for VTK should have 6 nodes
        for (unsigned elem_index=0; elem_index<p_mesh_for_vtk->GetNumElements(); elem_index++)
        {
            TS_ASSERT_EQUALS(p_mesh_for_vtk->GetElement(elem_index)->GetNumNodes(), 6u);
        }

        VertexElement<2, 2>* p_element0 = p_mesh_for_vtk->GetElement(0);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(2), 9u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(3), 12u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(4), 8u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(5), 4u);

        VertexElement<2, 2>* p_element3 = p_mesh_for_vtk->GetElement(3);
        TS_ASSERT_EQUALS(p_element3->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(p_element3->GetNodeGlobalIndex(1), 32u);
        TS_ASSERT_EQUALS(p_element3->GetNodeGlobalIndex(2), 33u);
        TS_ASSERT_EQUALS(p_element3->GetNodeGlobalIndex(3), 15u);
        TS_ASSERT_EQUALS(p_element3->GetNodeGlobalIndex(4), 11u);
        TS_ASSERT_EQUALS(p_element3->GetNodeGlobalIndex(5), 7u);

        VertexElement<2, 2>* p_element7 = p_mesh_for_vtk->GetElement(7);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(0), 33u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(1), 34u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(2), 35u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(3), 36u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(4), 19u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(5), 15u);

        VertexElement<2, 2>* p_element12 = p_mesh_for_vtk->GetElement(12);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(0), 25u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(1), 29u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(2), 40u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(3), 43u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(4), 39u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(5), 28u);

        VertexElement<2, 2>* p_element15 = p_mesh_for_vtk->GetElement(15);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(0), 37u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(1), 38u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(2), 46u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(3), 47u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(4), 42u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(5), 31u);
    }

    void TestConstructFromMeshReader()
    {
        // Create mesh reader and mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/toroidal_2d_mesh");
        Toroidal2dVertexMesh mesh;

        // Check nothing is thrown when width and height are input correctly
        TS_ASSERT_THROWS_NOTHING(mesh.ConstructFromMeshReader(mesh_reader, 4.0, 3.0));

        // Check exception is tripped when width and height are input incorrectly
        TS_ASSERT_THROWS_THIS(mesh.ConstructFromMeshReader(mesh_reader, 4.0, 2.0),
                              "Mesh width and height do not match sheet surface area.");
    }
};

#endif /*TESTTOROIDAL2DVERTEXMESH_HPP_*/
