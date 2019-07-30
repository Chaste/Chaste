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


#ifndef _TESTELEMENT_HPP_
#define _TESTELEMENT_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>

#include "MutableMesh.hpp"
#include "Exception.hpp"
#include "TrianglesMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

typedef BoundaryElement<2,3> BOUNDARY_ELEMENT_2_3;

class TestElement : public CxxTest::TestSuite
{
    /**
     * Create a node with given index that has a point at the origin and is
     * not on the boundary.
     *
     * @param index The global index of the created node.
     * @return A pointer to a new 3d node.
     */
    Node<3>* CreateZeroPointNode(int index)
    {
        c_vector<double,3> zero_point;
        zero_point.clear();

        Node<3>* p_node = new Node<3>(index, zero_point, false);
        return p_node;
    }

public:

    void TestConstructionForLinearBasisFunctions()
    {
        std::vector<Node<3>*> corner_nodes;
        corner_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        corner_nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(INDEX_IS_NOT_USED, corner_nodes);

        // Check nodes on the new element have the right indices
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(i), i);
        }

        c_matrix<double, 3, 3> jacob;
        double det;

        element.CalculateJacobian(jacob, det);

        TS_ASSERT_DELTA(det, 1.0, 1e-5);
        TS_ASSERT_DELTA(element.GetVolume(det), 1.0/6.0, 1e-5);

        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
        }

        element.SetIndex(27);
        TS_ASSERT_EQUALS(element.GetIndex(), 27u);
    }

    void TestEquals()
    {
        std::vector<Node<3>*> corner_nodes;
        corner_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        corner_nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(INDEX_IS_NOT_USED, corner_nodes);

        std::vector<Node<3>*> more_nodes;
        more_nodes.push_back(new Node<3>(0, false, 10.0, 10.0, 10.0));
        more_nodes.push_back(new Node<3>(1, false, 11.0, 10.0, 10.0));
        more_nodes.push_back(new Node<3>(2, false, 10.0, 11.0, 10.0));
        more_nodes.push_back(new Node<3>(3, false, 10.0, 10.0, 11.0));
        Element<3,3> another_element(INDEX_IS_NOT_USED, more_nodes);

        // Test (and cover) equals operator
        another_element = element;

        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(another_element.GetNode(i)->GetPoint()[j], element.GetNode(i)->GetPoint()[j], 1e-10);
            }
        }

        // Coverage of GetVolume()
        element.MarkAsDeleted();
        TS_ASSERT_DELTA(element.GetVolume(1.0), 0.0, 1e-6);

        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
            delete more_nodes[i];
        }
    }

    void TestGetSetAbstractTetrahedralElementMethods()
    {
        std::vector<Node<3>*> corner_nodes;
        corner_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        corner_nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(INDEX_IS_NOT_USED, corner_nodes);

        element.SetOwnership(true);

        TS_ASSERT_EQUALS(element.GetOwnership(), true);

        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
        }
    }

    void TestJacobian()
    {
        // 1D case
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);
        c_matrix<double, 1, 1> J1d;
        double det_1d;
        c_matrix<double, 1, 1> J1dinv;
        element1d.CalculateInverseJacobian(J1d, det_1d, J1dinv);

        TS_ASSERT_DELTA(J1d(0,0), 0.5, 1e-12);
        TS_ASSERT_DELTA(element1d.GetVolume(det_1d), 0.5, 1e-5);
        TS_ASSERT_DELTA(det_1d, 0.5, 1e-12);
        TS_ASSERT_DELTA(J1dinv(0,0), 2.0, 1e-12);

        delete nodes1d[0];
        delete nodes1d[1];

        // Easy 2D case
        std::vector<Node<2>*> nodes2d;
        nodes2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes2d.push_back(new Node<2>(2, false, 0.0, 1.0));
        Element<2,2> element2d(INDEX_IS_NOT_USED, nodes2d);
        c_matrix<double, 2, 2> J2d;
        double det;
        element2d.CalculateJacobian(J2d, det);

        TS_ASSERT_DELTA(J2d(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(J2d(0,1), 0.0, 1e-12);
        TS_ASSERT_DELTA(J2d(1,0), 0.0, 1e-12);
        TS_ASSERT_DELTA(J2d(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(det, 1.0, 1e-5);
        TS_ASSERT_DELTA(element2d.GetVolume(det), 0.5, 1e-5);

        delete nodes2d[0];
        delete nodes2d[1];
        delete nodes2d[2];

        // General 2D case
        std::vector<Node<2>*> nodes2d2;
        nodes2d2.push_back(new Node<2>(0, false, 1.0, -2.0));
        nodes2d2.push_back(new Node<2>(1, false, 4.0, -3.0));
        nodes2d2.push_back(new Node<2>(2, false, 2.0, -1.0));
        Element<2,2> element2d2(INDEX_IS_NOT_USED, nodes2d2);
        c_matrix<double, 2, 2> J2d2;
        c_matrix<double, 2, 2> J2d2inv;
        double det_2d;
        element2d2.CalculateInverseJacobian(J2d2, det_2d,J2d2inv);

        TS_ASSERT_DELTA(J2d2(0,0), 3.0, 1e-12);
        TS_ASSERT_DELTA(J2d2(0,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(J2d2(1,0), -1.0, 1e-12);
        TS_ASSERT_DELTA(J2d2(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(det_2d, 4.0, 1e-12);
        TS_ASSERT_DELTA(J2d2inv(0,0), 0.25, 1e-12);
        TS_ASSERT_DELTA(J2d2inv(0,1), -0.25, 1e-12);
        TS_ASSERT_DELTA(J2d2inv(1,0), 0.25, 1e-12);
        TS_ASSERT_DELTA(J2d2inv(1,1), 0.75, 1e-12);

        delete nodes2d2[0];
        delete nodes2d2[1];
        delete nodes2d2[2];

        // Easy 3D case
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element3d(INDEX_IS_NOT_USED, nodes3d);
        c_matrix<double, 3, 3> J3d;
        double det_3d;
        element3d.CalculateJacobian(J3d,det_3d);

        TS_ASSERT_DELTA(J3d(0,0), 1.0, 1e-12);
        TS_ASSERT_DELTA(J3d(0,1), 0.0, 1e-12);
        TS_ASSERT_DELTA(J3d(0,2), 0.0, 1e-12);
        TS_ASSERT_DELTA(J3d(1,0), 0.0, 1e-12);
        TS_ASSERT_DELTA(J3d(1,1), 1.0, 1e-12);
        TS_ASSERT_DELTA(J3d(1,2), 0.0, 1e-12);
        TS_ASSERT_DELTA(J3d(2,0), 0.0, 1e-12);
        TS_ASSERT_DELTA(J3d(2,1), 0.0, 1e-12);
        TS_ASSERT_DELTA(J3d(2,2), 1.0, 1e-12);

        delete nodes3d[0];
        delete nodes3d[1];
        delete nodes3d[2];
        delete nodes3d[3];

        // General 3D case
        std::vector<Node<3>*> nodes3d2;
        nodes3d2.push_back(new Node<3>(0, false, 1.0, 2.0, 3.0));
        nodes3d2.push_back(new Node<3>(1, false, 2.0, 1.0, 3.0));
        nodes3d2.push_back(new Node<3>(2, false, 5.0, 5.0, 5.0));
        nodes3d2.push_back(new Node<3>(3, false, 0.0, 3.0, 4.0));
        Element<3,3> element3d2(INDEX_IS_NOT_USED, nodes3d2);
        c_matrix<double, 3, 3> J3d2;
        c_matrix<double, 3, 3> J3d2inv;
        element3d2.CalculateInverseJacobian(J3d2,det_3d,J3d2inv);

        TS_ASSERT_DELTA(J3d2(0,0), 1.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(0,1), 4.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(0,2), -1.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(1,0), -1.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(1,1), 3.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(1,2), 1.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(2,0), 0.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(2,1), 2.0, 1e-4);
        TS_ASSERT_DELTA(J3d2(2,2), 1.0, 1e-4);

        // Check Jacobian determinant
        TS_ASSERT_DELTA(det_3d, 7.0, 1e-4);

        // Check inverse Jacobian
        TS_ASSERT_DELTA(J3d2inv(0,0), 1.0/7.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(0,1), -6.0/7.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(0,2), 1.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(1,0), 1.0/7.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(1,1), 1.0/7.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(1,2), 0.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(2,0), -2.0/7.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(2,1), -2.0/7.0, 1e-4);
        TS_ASSERT_DELTA(J3d2inv(2,2), 1.0, 1e-4);

        delete nodes3d2[0];
        delete nodes3d2[1];
        delete nodes3d2[2];
        delete nodes3d2[3];
    }

    void TestNodeToElementConversion()
    {
        ChastePoint<1> point1(1.0);
        ChastePoint<2> point2(2.0, -1.0);

        Node<1> node1(0, point1);
        Node<2> node2(0, point2);

        BoundaryElement<0,1> element1(INDEX_IS_NOT_USED, &node1);
        BoundaryElement<0,2> element2(INDEX_IS_NOT_USED, &node2);

        TS_ASSERT_EQUALS(element1.GetNode(0)->GetPoint()[0], point1[0]);

        TS_ASSERT_EQUALS(element2.GetNode(0)->GetPoint()[0], point2[0]);
        TS_ASSERT_EQUALS(element2.GetNode(0)->GetPoint()[1], point2[1]);
    }

    void TestGetNodeLocation()
    {
        ChastePoint<2> point1(0.0, 1.0);
        ChastePoint<2> point2(4.0, 6.0);
        ChastePoint<2> point3(2.0, 3.0);

        std::vector<Node<2>*> element_nodes;
        element_nodes.push_back(new Node<2>(0, point1));
        element_nodes.push_back(new Node<2>(0, point2));
        element_nodes.push_back(new Node<2>(0, point3));

        Element<2,2> element(INDEX_IS_NOT_USED, element_nodes);

        // Note that nodes 2 and 3 are swapped by the element constructor
        // to ensure that the Jacobian determinant is positive
        TS_ASSERT_EQUALS(element.GetNodeLocation(0,0), 0.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(0)(0), 0.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(1,0), 2.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(1)(0), 2.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(2,0), 4.0);
        TS_ASSERT_EQUALS(element.GetNodeLocation(2)(0), 4.0);

        delete element_nodes[0];
        delete element_nodes[1];
        delete element_nodes[2];
    }

    void TestElementSwapsNodesIfJacobianIsNegative()
    {
        ChastePoint<1> a0(0),     a1(1);
        ChastePoint<2> b0(0,0),   b1(1,0),   b2(0,1);
        ChastePoint<3> c0(0,0,0), c1(1,0,0), c2(0,1,0), c3(0,0,1);

        Node<1>  na0(0,a0), na1(1,a1);
        Node<2>  nb0(0,b0), nb1(1,b1), nb2(2,b2);
        Node<3>  nc0(0,c0), nc1(1,c1), nc2(2,c2), nc3(3,c3);

        ////////////////////////////////////////////
        // 1D case
        ////////////////////////////////////////////
        std::vector<Node<1>*> nodes_1d_correct;
        nodes_1d_correct.push_back(&na0);
        nodes_1d_correct.push_back(&na1);

        std::vector<Node<1>*> nodes_1d_incorrect;
        nodes_1d_incorrect.push_back(&na1);
        nodes_1d_incorrect.push_back(&na0);

        Element<1,1> e_1d_correct_orientation(INDEX_IS_NOT_USED, nodes_1d_correct);
        Element<1,1> e_1d_incorrect_orientation(INDEX_IS_NOT_USED, nodes_1d_incorrect);

        // Index of second node should be 1
        TS_ASSERT_EQUALS( e_1d_correct_orientation.GetNode(1)->GetIndex(), 1u);

        // Index of second node for incorrect orientation element should also be 1
        // because the element should have swapped the nodes around
        TS_ASSERT_EQUALS( e_1d_incorrect_orientation.GetNode(1)->GetIndex(), 1u);

        ////////////////////////////////////////////
        // 2D case
        ////////////////////////////////////////////
        std::vector<Node<2>*> nodes_2d_correct;
        nodes_2d_correct.push_back(&nb0);
        nodes_2d_correct.push_back(&nb1);
        nodes_2d_correct.push_back(&nb2);

        std::vector<Node<2>*> nodes_2d_incorrect;
        nodes_2d_incorrect.push_back(&nb1);
        nodes_2d_incorrect.push_back(&nb0);
        nodes_2d_incorrect.push_back(&nb2);

        Element<2,2>   e_2d_correct_orientation(INDEX_IS_NOT_USED, nodes_2d_correct);
        Element<2,2> e_2d_incorrect_orientation(INDEX_IS_NOT_USED, nodes_2d_incorrect);

        // Index of last node should be 2
        TS_ASSERT_EQUALS(e_2d_correct_orientation.GetNode(2)->GetIndex(), 2u);

        // Index of last node for incorrect orientation element should be 0
        // because the element should have swapped the last two nodes around
        TS_ASSERT_EQUALS(e_2d_incorrect_orientation.GetNode(2)->GetIndex(), 0u);

        ////////////////////////////////////////////
        // 3D case
        ////////////////////////////////////////////
        std::vector<Node<3>*> nodes_3d_correct;
        nodes_3d_correct.push_back(&nc0);
        nodes_3d_correct.push_back(&nc1);
        nodes_3d_correct.push_back(&nc2);
        nodes_3d_correct.push_back(&nc3);

        std::vector<Node<3>*> nodes_3d_incorrect;
        nodes_3d_incorrect.push_back(&nc0);
        nodes_3d_incorrect.push_back(&nc1);
        nodes_3d_incorrect.push_back(&nc3);
        nodes_3d_incorrect.push_back(&nc2);

        Element<3,3> e_3d_correct_orientation(INDEX_IS_NOT_USED, nodes_3d_correct);
        Element<3,3> e_3d_incorrect_orientation(INDEX_IS_NOT_USED, nodes_3d_incorrect);

        // Index of last node should be 3
        TS_ASSERT_EQUALS( e_3d_correct_orientation.GetNode(3)->GetIndex(), 3u);

        // Index of last node for incorrect orientation element should be 3
        // because the element should have swapped the last two nodes around
        TS_ASSERT_EQUALS( e_3d_incorrect_orientation.GetNode(3)->GetIndex(), 3u);
    }

    void TestElementCopyConstructor()
    {
        std::vector<Node<3>*> corner_nodes;
        corner_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        corner_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 3.0));
        corner_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        corner_nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element(31415, corner_nodes);

        // Create a copy of the element and test that it's the same as the original one
        Element<3,3> copied_element(element);

        TS_ASSERT_EQUALS(element.GetNumNodes(), copied_element.GetNumNodes());

        // Check that this version of the copy constructor gives the copied
        // element the same index number
        TS_ASSERT_EQUALS(copied_element.GetIndex(), 31415u);

        for (int node=0; node<4; node++)
        {
            for (int dimension=0; dimension<3; dimension++)
            {
                TS_ASSERT_EQUALS(element.GetNodeLocation(node, dimension), copied_element.GetNodeLocation(node, dimension));
            }
        }

        Element<3,3> another_copied_element(element, 2345);
        TS_ASSERT_EQUALS(another_copied_element.GetIndex(), 2345u);

        // Update a node of the element
        Node<3>* update_for_node= new Node<3>(4, false, 0.0, 0.0, 2.0);
        another_copied_element.UpdateNode(1, update_for_node);
        TS_ASSERT_EQUALS(another_copied_element.GetNodeLocation(1, 2), 2.0);

        for (unsigned i=0; i<corner_nodes.size(); i++)
        {
            delete corner_nodes[i];
        }
        delete update_for_node;
    }

    void TestBoundaryElement()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));

        BoundaryElement<2,3> element(INDEX_IS_NOT_USED, nodes);
        c_vector<double, 3> weighted_direction;
        double det;
        element.CalculateWeightedDirection(weighted_direction, det);
        TS_ASSERT_DELTA(det, 1.0, 1e-6);

        // Alter to be collinear (for coverage)
        nodes[2]->rGetModifiableLocation()[1] = 0.0;
        TS_ASSERT_THROWS_THIS(element.CalculateWeightedDirection(weighted_direction, det),"Jacobian determinant is zero");

        // Alter to be deleted (for coverage)
        element.MarkAsDeleted();
        TS_ASSERT_THROWS_THIS(element.CalculateWeightedDirection(weighted_direction, det),"Attempting to Refresh a deleted element");

        //Construction of face with zero area
        std::vector<Node<3>*> nodes_the_same;
        nodes_the_same.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes_the_same.push_back(new Node<3>(1, false, 0.0, 0.0, 0.0));
        nodes_the_same.push_back(new Node<3>(2, false, 0.0, 0.0, 0.0));
        TS_ASSERT_THROWS_THIS(new BOUNDARY_ELEMENT_2_3(INDEX_IS_NOT_USED, nodes_the_same), "Jacobian determinant is zero");


        Node<3>* p_fake_node = new Node<3>(0, false, 0.0, 0.0, 0.0);
        Node<3>* p_fake_node2 = new Node<3>(0, false, 0.0, 0.0, 0.0);
        TS_ASSERT_THROWS_THIS(element.ReplaceNode(p_fake_node, p_fake_node2),"You didn\'t have that node to start with.");
        delete p_fake_node;
        delete p_fake_node2;

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
            delete nodes_the_same[i];
        }

        // Test attribute when none is set

//        TS_ASSERT_THROWS_THIS(element.GetAttribute(),
//                              "Attempting to get element attribute when there are none defined");
        TS_ASSERT_DELTA(element.GetAttribute(), 0.00, 1e-9);


        // Test attribute setting and conversion to unsigned integers

        element.SetAttribute(3);
        TS_ASSERT_DELTA(element.GetAttribute(), 3.00, 1e-9);

        element.SetAttribute(3.0);
        TS_ASSERT_EQUALS(element.GetUnsignedAttribute(), 3u);

        //Check that rounding happens correctly
        element.SetAttribute(2.99999999999999999);
        TS_ASSERT_EQUALS(element.GetUnsignedAttribute(), 3u);

        element.SetAttribute(3.1);
        TS_ASSERT_THROWS_THIS(element.GetUnsignedAttribute(),
                              "Element attribute '3.1' cannot be converted to an unsigned.");

        element.SetAttribute(0);
        TS_ASSERT_EQUALS(element.GetUnsignedAttribute(), 0u);

        element.SetAttribute(-1.0);
        TS_ASSERT_THROWS_THIS(element.GetUnsignedAttribute(),
                              "Element attribute '-1' cannot be converted to an unsigned.");

        element.SetAttribute(-1.2);
        TS_ASSERT_THROWS_THIS(element.GetUnsignedAttribute(),
                              "Element attribute '-1.2' cannot be converted to an unsigned.");

    }

    void TestCircum1d()
    {
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, false, 10.0));
        nodes.push_back(new Node<1>(1, false, 15.0));

        Element<1,1> element(0, nodes);

        c_matrix<double, 1, 1> jacobian, inverse_jacobian;
        double det;
        element.CalculateInverseJacobian(jacobian, det, inverse_jacobian);
        c_vector<double, 2> circum = element.CalculateCircumsphere(jacobian, inverse_jacobian);
        TS_ASSERT_DELTA(circum[0], 12.5, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum[1]), 2.5, 1e-7);

        TS_ASSERT_DELTA(element.CalculateQuality(), 1.0, 1e-7);

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCircum2d()
    {
        std::vector<Node<2>*> equilateral_nodes;
        equilateral_nodes.push_back(new Node<2>(0, false, 2.0,   0.0));
        equilateral_nodes.push_back(new Node<2>(1, false, -1.0,  sqrt(3.0)));
        equilateral_nodes.push_back(new Node<2>(2, false, -1.0, -sqrt(3.0)));

        Element<2,2> equilateral_element(0, equilateral_nodes);

        c_matrix<double, 2, 2> jacobian, inverse_jacobian;
        double det;

        equilateral_element.CalculateInverseJacobian(jacobian, det, inverse_jacobian);
        c_vector<double, 3> circum = equilateral_element.CalculateCircumsphere(jacobian, inverse_jacobian);
        TS_ASSERT_DELTA(circum[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(circum[1], 0.0, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum[2]), 2.0, 1e-7);

        TS_ASSERT_DELTA(equilateral_element.CalculateQuality(), 1.0, 1e-7);

        std::vector<Node<2>*> right_angle_nodes;
        right_angle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        right_angle_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        right_angle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        Element<2,2> right_angle_element(0, right_angle_nodes);

        right_angle_element.CalculateInverseJacobian(jacobian, det, inverse_jacobian);
        c_vector<double, 3> circum2 = right_angle_element.CalculateCircumsphere(jacobian, inverse_jacobian);
        TS_ASSERT_DELTA(circum2[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(circum2[1], 0.5, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum2[2]), 1.0/sqrt(2.0), 1e-7);

        TS_ASSERT_DELTA(right_angle_element.CalculateQuality(), 4.0*sqrt(3.0)/9.0, 1e-7);


        c_vector<double, 2> equilateral_min_max = equilateral_element.CalculateMinMaxEdgeLengths();
        TS_ASSERT_DELTA(equilateral_min_max[0], 2*sqrt(3.0), 1e-5);
        TS_ASSERT_DELTA(equilateral_min_max[1], 2*sqrt(3.0), 1e-5);
        c_vector<double, 2> right_angle_min_max = right_angle_element.CalculateMinMaxEdgeLengths();
        TS_ASSERT_DELTA(right_angle_min_max[0], 1.0, 1e-5);
        TS_ASSERT_DELTA(right_angle_min_max[1], sqrt(2.0), 1e-5);


        for (unsigned i=0; i<equilateral_nodes.size(); i++)
        {
            delete equilateral_nodes[i];
        }
        for (unsigned i=0; i<right_angle_nodes.size(); i++)
        {
            delete right_angle_nodes[i];
        }
    }

    void TestCircum3d()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  1.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(1, false, -1.0, -1.0,  1.0));
        nodes.push_back(new Node<3>(2, false, -1.0,  1.0, -1.0));
        nodes.push_back(new Node<3>(3, false,  1.0, -1.0, -1.0));

        Element<3,3> element(0, nodes);

        c_matrix<double, 3, 3> jacobian, inverse_jacobian;
        double det;
        element.CalculateInverseJacobian(jacobian, det, inverse_jacobian);
        c_vector<double, 4> circum = element.CalculateCircumsphere(jacobian, inverse_jacobian);
        TS_ASSERT_DELTA(circum[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(circum[1], 0.0, 1e-7);
        TS_ASSERT_DELTA(circum[2], 0.0, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum[3]), sqrt(3.0), 1e-7);

        TS_ASSERT_DELTA(element.CalculateQuality(), 1.0, 1e-7);

        std::vector<Node<3>*> right_angle_nodes;
        right_angle_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        right_angle_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        right_angle_nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        right_angle_nodes.push_back(new Node<3>(2, false, 0.0, 0.0, 1.0));

        Element<3,3> right_angle_element(0, right_angle_nodes);

        right_angle_element.CalculateInverseJacobian(jacobian, det, inverse_jacobian);
        c_vector<double, 4> circum2 = right_angle_element.CalculateCircumsphere(jacobian, inverse_jacobian);
        TS_ASSERT_DELTA(circum2[0], 0.5, 1e-7);
        TS_ASSERT_DELTA(circum2[1], 0.5, 1e-7);
        TS_ASSERT_DELTA(circum2[2], 0.5, 1e-7);
        TS_ASSERT_DELTA(sqrt(circum2[3]), sqrt(3.0)/2.0, 1e-7);

        TS_ASSERT_DELTA(right_angle_element.CalculateQuality(), 0.5, 1e-7);

        c_vector<double, 2> min_max = element.CalculateMinMaxEdgeLengths();
        TS_ASSERT_DELTA(min_max[0], 2*sqrt(2.0), 1e-5);
        TS_ASSERT_DELTA(min_max[1], 2*sqrt(2.0), 1e-5);
        c_vector<double, 2> right_angle_min_max = right_angle_element.CalculateMinMaxEdgeLengths();
        TS_ASSERT_DELTA(right_angle_min_max[0], 1.0, 1e-5);
        TS_ASSERT_DELTA(right_angle_min_max[1], sqrt(2.0), 1e-5);

        for (unsigned i=0; i<right_angle_nodes.size(); i++)
        {
            delete right_angle_nodes[i];
        }
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestNormal1dIn2d()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.5, 0.0));
        BoundaryElement<1,2> element_1d(0, nodes);

        c_vector<double, 2> normal = element_1d.CalculateNormal();

        TS_ASSERT_EQUALS(normal[0], 0.0);
        TS_ASSERT_EQUALS(normal[1], -1.0);

        delete nodes[0];
        delete nodes[1];
    }

    void TestCentroidAndDirection()
    {
        c_vector<double,3> direction;
        c_vector<double,3> centroid;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.5, 0.0, 0.0));
        BoundaryElement<1,3> element_1d(0, nodes);

        double det;
        element_1d.CalculateWeightedDirection(direction, det);

        // 1D element in higher space is orientated by vector between endpoints
        TS_ASSERT_EQUALS(direction[0], 0.5);
        TS_ASSERT_EQUALS(direction[1], 0.0);
        TS_ASSERT_EQUALS(direction[2], 0.0);
        TS_ASSERT_EQUALS(det, 0.5);

        // Can't compute sensible normal for 1d-in-3d
        TS_ASSERT_THROWS_THIS(element_1d.CalculateNormal(),
                "Don't have enough information to calculate a normal vector");

        centroid = element_1d.CalculateCentroid();
        TS_ASSERT_EQUALS(centroid[0], 0.25);
        TS_ASSERT_EQUALS(centroid[1], 0.0);
        TS_ASSERT_EQUALS(centroid[2], 0.0);

        nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        BoundaryElement<2,3> element_2d(0, nodes);

        element_2d.CalculateWeightedDirection(direction, det);

        // 2D element in higher space is oriented by a normal
        TS_ASSERT_EQUALS(direction[0], 0.0);
        TS_ASSERT_EQUALS(direction[1], 0.0);
        TS_ASSERT_EQUALS(direction[2], -0.5);
        TS_ASSERT_EQUALS(det, 0.5);

        c_vector<double,3> normal = element_2d.CalculateNormal();
        TS_ASSERT_EQUALS(normal[0], 0.0);
        TS_ASSERT_EQUALS(normal[1], 0.0);
        TS_ASSERT_EQUALS(normal[2], -1.0);

        // Note it's negative z-axis:
        // Vertices are *clockwise* when we look at them from outward facing normal
        centroid = element_2d.CalculateCentroid();
        TS_ASSERT_DELTA(centroid[0], 0.5/3.0, 1e-8);
        TS_ASSERT_DELTA(centroid[1], 1.0/3.0, 1e-8);
        TS_ASSERT_EQUALS(centroid[2], 0.0);

        nodes.push_back(new Node<3>(2, false, 0.0, 0.0, 1.0));
        Element<3,3> element_3d(0, nodes);
        // Weighted direction only makes sense in a subspace element
        TS_ASSERT_THROWS_THIS(element_3d.CalculateWeightedDirection(direction, det),
                "WeightedDirection undefined for fully dimensional element");

        // 3D element in 3D space has no orientation (other than JacobianDeterminant)
        centroid = element_3d.CalculateCentroid();
        TS_ASSERT_DELTA(centroid[0], 0.125, 1e-8);
        TS_ASSERT_DELTA(centroid[1], 0.25, 1e-8);
        TS_ASSERT_EQUALS(centroid[2], 0.25 );

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void Test1DRefineElement()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        BoundaryElement<0,1>* p_end_boundary = mesh.GetBoundaryElement(0);
        c_vector<double, 1> dir;
        double unit_det;
        p_end_boundary->CalculateWeightedDirection(dir, unit_det);
        TS_ASSERT_EQUALS(dir[0], 1.0);
        TS_ASSERT_EQUALS(unit_det, 1.0);
        dir = p_end_boundary->CalculateNormal();
        TS_ASSERT_EQUALS(dir[0], 0.0);

        ChastePoint<1> new_point(0.01);
        Element<1,1>* p_first_element = mesh.GetElement(0);


        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_first_element, new_point));

        // Instead of a element with nodes at 0 and 0.1
        // there should be an element with nodes at 0 and 0.01 and
        // an element with nodes at 0.01 and 0.1. Other elements should stay the same

        const Node<1>* p_first_node = p_first_element->GetNode(0);
        const Node<1>* p_second_node = p_first_element->GetNode(1);

        TS_ASSERT_EQUALS(p_first_node->GetPoint().rGetLocation()(0), 0);
        TS_ASSERT_EQUALS(p_second_node->GetPoint().rGetLocation()(0), 0.01);

        // Test second element
        Element<1,1>* p_second_element = mesh.GetElement(1);
        p_first_node = p_second_element->GetNode(0);
        p_second_node = p_second_element->GetNode(1);

        TS_ASSERT_EQUALS(p_first_node->GetPoint().rGetLocation()(0), 0.1);
        TS_ASSERT_EQUALS(p_second_node->GetPoint().rGetLocation()(0), 0.2);

        // Test last element
        Element<1,1>* p_last_element = mesh.GetElement(10);
        p_first_node = p_last_element->GetNode(0);
        p_second_node = p_last_element->GetNode(1);

        TS_ASSERT_EQUALS(p_first_node->GetPoint().rGetLocation()(0), 0.01);
        TS_ASSERT_EQUALS(p_second_node->GetPoint().rGetLocation()(0), 0.1);

        // Test Jacobians
        c_matrix<double,1,1> jacobian;
        double det;
        p_first_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_EQUALS(det, 0.01);
        p_second_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_EQUALS(det, 0.1);
        p_last_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.09, 1e-6);

        mesh.RefreshMesh();
        // Test mesh length
        TS_ASSERT_DELTA(mesh.GetVolume(), 1.0, 1e-6);
    }

    void Test1DRefineElementWithPointTooFarRightFails()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<1> new_point(0.11);
        Element<1,1>* p_first_element = mesh.GetElement(0);

        // Trying to add Point(0.11) to Element(0)
        // This point is contained in Element(1)
        TS_ASSERT_THROWS_THIS(mesh.RefineElement(p_first_element, new_point),
                "RefineElement could not be started (point is not in element)");
    }

    void Test1DRefineElementWithPointTooFarLeftFails()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<1> new_point(-0.1);
        Element<1,1>* p_first_element = mesh.GetElement(0);

        // Trying to add Point(-0.1) to Element(0)
        // This point is to the left of Element(0)
        TS_ASSERT_THROWS_THIS(mesh.RefineElement(p_first_element, new_point),
                "RefineElement could not be started (point is not in element)");
    }

    void Test1DRefineElementManyNodes()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Element<1,1>* p_first_element = mesh.GetElement(0);

        // There's space on the node vector for 10 new points
        // but more than 10 should still work
        for (int i=1; i<=20; i++)
        {
            ChastePoint<1> new_point(0.1 - i*0.0005);
            mesh.RefineElement(p_first_element, new_point);
        }
    }

    void Test2DRefineElement()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA(mesh.GetVolume(), 0.01,1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 0.4,1e-6);

        // Refine an element in the bottom right corner
        Element<2,2>* p_corner_element = mesh.GetElement(18);

        c_matrix<double,2,2> jacobian;
        double det;
        p_corner_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.0001, 1e-6);

        // Point to be inserted in the bottom right corner
        ChastePoint<2> new_point(0.095, 0.003);

        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_corner_element, new_point));

        // Testing the JacobianDeterminant of the element that has changed
        p_corner_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 3e-05, 1e-6);

        // Testing invariants
        mesh.RefreshMesh();
        TS_ASSERT_DELTA(mesh.GetVolume(), 0.01, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 0.4, 1e-6);

        // Refine an element in the middle of the mesh
        Element<2,2>* p_middle_element = mesh.GetElement(108);

        p_middle_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.0001, 1e-6);

        // Point to be inserted in the middle of the mesh
        ChastePoint<2> new_point1(0.045,0.053);

        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_middle_element, new_point1));

        // Testing the JacobianDeterminant of the element that has changed
        p_middle_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 3e-05, 1e-6);

        // Testing invariants
        mesh.RefreshMesh();
        TS_ASSERT_DELTA(mesh.GetVolume(), 0.01, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 0.4, 1e-6);
    }

    void Test2DRefineElementFails()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Refine an element in the bottom right corner
        Element<2,2>* p_corner_element = mesh.GetElement(18);

        // Point to be inserted on the edge of the element
        ChastePoint<2> new_point(0.095, 0.005);

        // Shouldn't be able to insert this point at the edge of the element
        TS_ASSERT_THROWS_THIS(mesh.RefineElement(p_corner_element, new_point),
                "RefineElement could not be started (point is not in element)");
    }

    void Test3DRefineElement()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA(mesh.GetVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 6.0, 1e-6);

        // Refine an element in the top corner (1, 1, 1)
        Element<3,3>* p_corner_element = mesh.GetElement(64);
        c_matrix<double,3,3> jacobian;
        double det;
        p_corner_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.03125, 1e-6);

        // Point to be inserted in the top corner
        ChastePoint<3> new_point(0.9, 0.75, 0.9);

        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_corner_element, new_point));

        // Testing the JacobianDeterminant of the element that has changed
        p_corner_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.0125, 1e-6);

        // Testing invariants
        mesh.RefreshMesh();
        TS_ASSERT_DELTA(mesh.GetVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 6.0, 1e-6);

        // Refine an element which includes the middle node
        Element<3,3>* p_middle_element = mesh.GetElement(49);
        p_middle_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.0625, 1e-6);

        // Point to be inserted near node 22 (middle of cube)
        ChastePoint<3> new_point1(0.49, 0.47, 0.6);

        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_middle_element, new_point1));

        // Testing the JacobianDeterminant of the element that has changed
        p_middle_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.01125, 1e-6);

        // Testing invariants
        mesh.RefreshMesh();
        TS_ASSERT_DELTA(mesh.GetVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 6.0, 1e-6);
    }

    void Test3DRefineElementFails()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA(mesh.GetVolume(), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetSurfaceArea(), 6.0, 1e-6);

        // Refine an element in the top corner (1, 1, 1)
        Element<3,3>* p_corner_element = mesh.GetElement(64);
        c_matrix<double,3,3> jacobian;
        double det;
        p_corner_element->CalculateJacobian(jacobian, det);
        TS_ASSERT_DELTA(det, 0.03125, 1e-6);

        // Point to be inserted in wrong place
        ChastePoint<3> new_point(0.9, 0.75, 1.0);

        TS_ASSERT_THROWS_THIS(mesh.RefineElement(p_corner_element, new_point),
                "RefineElement could not be started (point is not in element)");
    }

    void TestGetStiffnessMatrixGlobalIndices()
    {
        unsigned PROBLEM_DIM = 1;

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Element<2,2>* p_element = mesh.GetElement(2);

        unsigned indices_1[3];

        // Element 2 in the mesh has global node indices: 1,2,4
        unsigned expected_indices_1[] = {1u, 2u, 4u};
        p_element->GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, indices_1);
        TS_ASSERT_SAME_DATA(indices_1, expected_indices_1, 3*sizeof(unsigned));

        PROBLEM_DIM = 2;
        unsigned indices_2[6];
        unsigned expected_indices_2[] = {2u, 3u, 4u, 5u, 8u, 9u};

        p_element->GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, indices_2);

        TS_ASSERT_SAME_DATA(indices_2, expected_indices_2, 6*sizeof(unsigned));

        p_element->SetAttribute(3);
        TS_ASSERT_EQUALS(p_element->GetUnsignedAttribute(), 3u);
    }
};

#endif //_TESTELEMENT_HPP_
