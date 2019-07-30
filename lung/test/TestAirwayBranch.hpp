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

#ifndef _TESTAIRWAYBRANCH_HPP_
#define _TESTAIRWAYBRANCH_HPP_

#include <cxxtest/TestSuite.h>

#include "AirwayBranch.hpp"

class TestAirwayBranch : public CxxTest::TestSuite
{
public:
    void TestMakeBranch()
    {
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/three_generation_branch_mesh_refined");
        mesh.ConstructFromMeshReader(mesh_reader);

        AirwayBranch branch_one;
        branch_one.AddElement(mesh.GetElement(0));

        TS_ASSERT_DELTA(branch_one.GetLength(), 0.05, 1e-6);
        TS_ASSERT_DELTA(branch_one.GetAverageRadius(), (75+62.5)/2.0, 1e-3);
        TS_ASSERT_DELTA(branch_one.GetDirection()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(branch_one.GetDirection()[1], 0.0, 1e-3);
        TS_ASSERT_DELTA(branch_one.GetDirection()[2], -1.0, 1e-3);

        AirwayBranch branch_two;
        branch_two.AddElement(mesh.GetElement(0));
        branch_two.AddElement(mesh.GetElement(1));
        branch_two.AddElement(mesh.GetElement(2));

        TS_ASSERT_DELTA(branch_two.GetLength(), 0.1*(1.0 + 0.5*sqrt(2.0)), 1e-6);
        TS_ASSERT_DELTA(branch_two.GetAverageRadius(), 54.7334, 1e-3);

        c_vector<double, 3> dir;
        dir[0] = 0.05;
        dir[1] = 0.05;
        dir[2] = -0.1;
        dir = dir/norm_2(dir);
        TS_ASSERT_DELTA(branch_two.GetDirection()[0], dir[0], 1e-3);
        TS_ASSERT_DELTA(branch_two.GetDirection()[1], dir[1], 1e-3);
        TS_ASSERT_DELTA(branch_two.GetDirection()[2], dir[2], 1e-3);

        branch_one.SetChildOne(&branch_two);
        branch_one.SetChildTwo(&branch_one);

        TS_ASSERT_EQUALS(branch_one.GetChildOne(), &branch_two);
        TS_ASSERT_EQUALS(branch_one.GetChildTwo(), &branch_one);

        branch_two.SetParent(&branch_one);
        TS_ASSERT_EQUALS(branch_two.GetParent(), &branch_one);

        branch_one.SetSibling(&branch_two);
        TS_ASSERT_EQUALS(branch_one.GetSibling(), &branch_two);

        TS_ASSERT_EQUALS(branch_one.GetElements().size(), 1u);
        TS_ASSERT_EQUALS(branch_two.GetElements().size(), 3u);

        TS_ASSERT_EQUALS(branch_two.IsMajor(), true);
        branch_two.SetSibling(&branch_one);
        TS_ASSERT_EQUALS(branch_two.IsMajor(), false);

        TS_ASSERT_EQUALS(branch_one.IsMajor(), true);

        //Check radius on edge calculations
        AirwayBranch branch_three(true);
        mesh.GetElement(0)->SetAttribute(0.5);
        mesh.GetElement(1)->SetAttribute(0.5);
        mesh.GetElement(2)->SetAttribute(0.5);
        branch_three.AddElement(mesh.GetElement(0));
        branch_three.AddElement(mesh.GetElement(1));
        branch_three.AddElement(mesh.GetElement(2));

        TS_ASSERT_DELTA(branch_three.GetAverageRadius(), 0.5, 1e-3);
        TS_ASSERT_DELTA(branch_three.GetPoiseuilleResistance(), 2.7313, 1e-3);
    }

    void TestAngleCalculations()
    {
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/branched_1d_in_3d_mesh");
        mesh.ConstructFromMeshReader(mesh_reader);

        AirwayBranch branch_one;
        AirwayBranch branch_two;
        AirwayBranch branch_three;

        for (unsigned i = 0; i < 10; ++i)
        {
            branch_one.AddElement(mesh.GetElement(i));
            branch_two.AddElement(mesh.GetElement(10 + i));
            branch_three.AddElement(mesh.GetElement(20 + i));
        }

        TS_ASSERT_THROWS_CONTAINS(branch_one.GetBranchAngle(), "Insufficient airway tree structure to calculate branch angle");
        TS_ASSERT_THROWS_CONTAINS(branch_one.GetRotationAngle(), "Insufficient airway tree structure to calculate rotation angle");

        branch_one.SetChildOne(&branch_two);
        branch_one.SetChildTwo(&branch_three);
        branch_two.SetParent(&branch_one);
        branch_three.SetParent(&branch_one);
        branch_two.SetSibling(&branch_three);
        branch_three.SetSibling(&branch_two);
        branch_one.SetSibling(&branch_two);

        TS_ASSERT_DELTA(branch_two.GetBranchAngle(), M_PI/2.0, 1e-6);
        TS_ASSERT_DELTA(branch_two.GetRotationAngle(), M_PI/2.0, 1e-6);

        TS_ASSERT(!branch_one.IsTerminal());
        TS_ASSERT(branch_two.IsTerminal());
        TS_ASSERT(branch_three.IsTerminal());

        //Check special case of a zero rotation angle.
        branch_one.SetChildOne(&branch_two);
        branch_one.SetChildTwo(&branch_one);
        branch_one.SetSibling(&branch_two);
        branch_two.SetSibling(&branch_one);

        TS_ASSERT_DELTA(branch_two.GetRotationAngle(), 0.0, 1e-6);
    }

    void TestBranchIndexing()
    {
        AirwayBranch branch;

        // Branch index should default to UINT_MAX
        TS_ASSERT_EQUALS(branch.GetIndex(), UINT_MAX);

        // Test index set method
        branch.SetIndex(5u);
        TS_ASSERT_EQUALS(branch.GetIndex(), 5u);
    }

    void TestGetProximalAndDistalNodes()
    {
        // Load mesh
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/three_generation_branch_mesh_refined");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Add elements to branch
        AirwayBranch branch;
        branch.AddElement(mesh.GetElement(0));
        branch.AddElement(mesh.GetElement(1));
        branch.AddElement(mesh.GetElement(2));

        // Get proximal and distal nodes using member functions
        Node<3>* p_proximal_node = branch.GetProximalNode();
        Node<3>* p_distal_node = branch.GetDistalNode();

        // Get first and last elements from branch
        Element<1, 3>* p_first_element = branch.GetElements().front();
        Element<1, 3>* p_last_element = branch.GetElements().back();

        // Check each element has precisely two nodes
        TS_ASSERT_EQUALS(p_first_element->GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(p_last_element->GetNumNodes(), 2u);

        // Get nodes from first element
        Node<3>* p_first_elem_node_0 = p_first_element->GetNode(0);
        Node<3>* p_first_elem_node_1 = p_first_element->GetNode(1);

        // Get nodes from last element
        Node<3>* p_last_elem_node_0 = p_last_element->GetNode(0);
        Node<3>* p_last_elem_node_1 = p_last_element->GetNode(1);

        // Set up non-proximal and non-distal nodes
        Node<3>* p_first_elem_non_proximal;
        Node<3>* p_last_elem_non_distal;

        // Check that proximal node is indeed in the first element
        bool proximal_is_node_0 = (p_proximal_node->GetIndex() == p_first_elem_node_0->GetIndex());
        bool proximal_is_node_1 = (p_proximal_node->GetIndex() == p_first_elem_node_1->GetIndex());
        bool proximal_in_first_elem = proximal_is_node_0 || proximal_is_node_1;
        TS_ASSERT_EQUALS(proximal_in_first_elem, true);

        if (proximal_is_node_0)
        {
            p_first_elem_non_proximal = p_first_elem_node_1;
        }
        else
        {
            p_first_elem_non_proximal = p_first_elem_node_0;
        }

        // Check that distal node is indeed in the last element
        bool distsal_is_node_0 = (p_distal_node->GetIndex() == p_last_elem_node_0->GetIndex());
        bool distsal_is_node_1 = (p_distal_node->GetIndex() == p_last_elem_node_1->GetIndex());
        bool distal_in_last_elem = distsal_is_node_0 || distsal_is_node_1;
        TS_ASSERT_EQUALS(distal_in_last_elem, true);

        if (distsal_is_node_0)
        {
            p_last_elem_non_distal = p_last_elem_node_1;
        }
        else
        {
            p_last_elem_non_distal = p_last_elem_node_0;
        }

        // Get distance from proximal to distal node
        double proximal_to_distal = norm_2(p_proximal_node->rGetLocation() - p_distal_node->rGetLocation());

        // Get distance from non-proximal to distal
        double non_proximal_to_distal = norm_2(p_first_elem_non_proximal->rGetLocation() - p_distal_node->rGetLocation());

        // Get distance from non-distal to proximal
        double non_distal_to_proximal = norm_2(p_last_elem_non_distal->rGetLocation() - p_proximal_node->rGetLocation());

        // Check proximal to distal distance is bigger than other possibilities
        TS_ASSERT_LESS_THAN(non_proximal_to_distal, proximal_to_distal);
        TS_ASSERT_LESS_THAN(non_distal_to_proximal, proximal_to_distal);
    }

    void TestBranchProperties()
    {
        Node<3> node_a(0u, true, 0.0, 0.0, 0.0);
        Node<3> node_b(1u, true, 1.0, 0.0, 0.0);
        Node<3> node_c(2u, true, 4.0, 0.0, 0.0);
        Node<3> node_d(3u, true, 6.0, 0.0, 0.0);

        node_a.AddNodeAttribute(0.5);
        node_b.AddNodeAttribute(0.5);
        node_c.AddNodeAttribute(0.5);
        node_d.AddNodeAttribute(0.5);

        std::vector<Node<3>*> elem_a_nodes;
        elem_a_nodes.push_back(&node_a);
        elem_a_nodes.push_back(&node_b);

        std::vector<Node<3>*> elem_b_nodes;
        elem_b_nodes.push_back(&node_b);
        elem_b_nodes.push_back(&node_c);

        std::vector<Node<3>*> elem_c_nodes;
        elem_c_nodes.push_back(&node_c);
        elem_c_nodes.push_back(&node_d);

        Element<1,3> element_a(0u, elem_a_nodes, true);
        Element<1,3> element_b(1u, elem_b_nodes, true);
        Element<1,3> element_c(2u, elem_c_nodes, true);

        AirwayBranch branch;
        branch.AddElement(&element_a);
        branch.AddElement(&element_b);
        branch.AddElement(&element_c);

        double branch_length = branch.GetLength();
        double branch_resistance = branch.GetPoiseuilleResistance();
        double branch_volume = branch.GetBranchVolume();
        double branch_surface = branch.GetBranchLateralSurfaceArea();
        c_vector<double, 3> branch_centroid = branch.GetBranchCentroid();

        TS_ASSERT_DELTA(branch_length, 6.0, 1e-6);
        TS_ASSERT_DELTA(branch_resistance, 6.0 / (0.5 * 0.5 * 0.5 * 0.5), 1e-6);
        TS_ASSERT_DELTA(branch_volume, 6.0 * M_PI * 0.5 * 0.5, 1e-6);
        TS_ASSERT_DELTA(branch_surface, 6.0 * 2 * M_PI * 0.5, 1e-6);
        TS_ASSERT_DELTA(norm_2(branch_centroid), 3.0, 1e-6);


    }
};

#endif /*_TESTAIRWAYBRANCH_HPP_*/
