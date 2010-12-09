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

#ifndef TESTPOTTSELEMENT_HPP_
#define TESTPOTTSELEMENT_HPP_

#include <cxxtest/TestSuite.h>

#include "PottsElement.hpp"

class TestPottsElement : public CxxTest::TestSuite
{
public:

    void TestSimple()
    {
        TS_ASSERT(true);
    }

    // These are the sort of tests we need

//    void Test2dPottsElement()
//    {
//        std::vector<Node<2>*> nodes;
//        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
//
//        VertexElement<2,2> element(0, nodes);
//
//        // Test RegisterWithNodes()
//        element.RegisterWithNodes();
//
//        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
//        {
//            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 1u);
//        }
//
//        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
//        TS_ASSERT_EQUALS(element.GetNode(0)->GetIndex(), 0u);
//        TS_ASSERT_EQUALS(element.GetNode(1)->GetIndex(), 1u);
//
//        // Test GetNodeLocalIndex()
//        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), 0u);
//        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), 1u);
//
//        // Test UpdateNode()
//        Node<2>* p_node_2 = new Node<2>(2, false, 1.2, 1.3);
//        element.UpdateNode(0, p_node_2);
//
//        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.2, 1e-12);
//        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 1.3, 1e-12);
//
//        // Test ResetIndex()
//        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
//        element.ResetIndex(5);
//        TS_ASSERT_EQUALS(element.GetIndex(), 5u);
//
//        // Test DeleteNode() and AddNode()
//        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
//        element.DeleteNode(1);
//        TS_ASSERT_EQUALS(element.GetNumNodes(), 1u);
//
//        Node<2>* p_node_3 = new Node<2>(3, false, 0.1, 0.4);
//        element.AddNode(0, p_node_3);
//        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
//
//        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.2, 1e-12);
//        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 1.3, 1e-12);
//        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[0], 0.1, 1e-12);
//        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[1], 0.4, 1e-12);
//
//        // Test GetNodeLocalIndex()
//        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), UINT_MAX);
//        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), UINT_MAX);
//
//        // Test MarkAsDeleted()
//        element.MarkAsDeleted();
//
//        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
//        {
//            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 0u);
//        }
//
//        // Tidy up
//        for (unsigned i=0; i<nodes.size(); i++)
//        {
//            delete nodes[i];
//        }
//        delete p_node_2;
//        delete p_node_3;
//    }
//
//    void TestMarkAsDeleted()
//    {
//        // Create nodes
//        std::vector<Node<2>*> nodes;
//        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        // Create element
//        VertexElement<2,2> vertex_element(0, nodes);
//        vertex_element.RegisterWithNodes();
//
//        for (unsigned i=0; i<nodes.size(); i++)
//        {
//            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 1u);
//        }
//
//        vertex_element.MarkAsDeleted();
//
//        for (unsigned i=0; i<nodes.size(); i++)
//        {
//            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 0u);
//        }
//
//        // Tidy up
//        for (unsigned i=0; i<nodes.size(); i++)
//        {
//            delete nodes[i];
//        }
//    }
//
//    void TestUpdateNode()
//    {
//        // Create nodes
//        std::vector<Node<2>*> nodes;
//        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        // Create element
//        VertexElement<2,2> vertex_element(0, nodes);
//        vertex_element.RegisterWithNodes();
//
//        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);
//        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);
//
//        // Update location of node 2
//        Node<2>* p_node = new Node<2>(4, false, 1.2, 1.3);
//        vertex_element.UpdateNode(2, p_node);
//
//        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.2, 1e-12);
//        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[1], 1.3, 1e-12);
//
//        // Tidy up
//        for (unsigned i=0; i<nodes.size(); ++i)
//        {
//            delete nodes[i];
//        }
//        delete p_node;
//    }
//
//    void TestGetNodeLocalIndex()
//    {
//        // Create nodes
//        std::vector<Node<2>*> nodes;
//
//        // This is a square
//        nodes.push_back(new Node<2>(3, false, 0.0, 0.0));
//        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
//        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
//        nodes.push_back(new Node<2>(0, false, 0.0, 1.0));
//
//        // Create element
//        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);
//
//        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0), 3u);
//        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(1), 2u);
//        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(2), 1u);
//        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(3), 0u);
//
//        vertex_element.DeleteNode(3); // Removes (1,1) node
//
//        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0), UINT_MAX);
//
//        // Tidy up
//        for (unsigned i=0; i<nodes.size(); i++)
//        {
//            delete nodes[i];
//        }
//    }
};

#endif /*TESTPOTTSELEMENT_HPP_*/
