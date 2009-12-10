/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef TESTVERTEXELEMENT_HPP_
#define TESTVERTEXELEMENT_HPP_

#include <cxxtest/TestSuite.h>

#include "VertexElement.hpp"


class TestVertexElement : public CxxTest::TestSuite
{
public:

    void TestVertexElementDeleteAndAddNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);

        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 6u);

        vertex_element.DeleteNode(3); // Removes (-1,0) node
        vertex_element.DeleteNode(0); // Removes (1,0) node

        // Test node is removed
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);

        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        // Add new node
        Node<2>* p_new_node = new Node<2>(4, false, 0.0, 0.0);
        vertex_element.AddNode(3, p_new_node); // Add node at (0,0) between nodes 3 and 0

        // Test node is added
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 5u);

        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[1], 0.0, 1e-9);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_new_node;
    }


    void TestMarkAsDeleted()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);
        vertex_element.RegisterWithNodes();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 1u);
        }

        vertex_element.MarkAsDeleted();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }


    void TestUpdateNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);
        vertex_element.RegisterWithNodes();

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);

        // Update location of node 2
        Node<2>* p_node = new Node<2>(4, false, 1.2, 1.3);
        vertex_element.UpdateNode(2, p_node);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.2, 1e-12);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[1], 1.3, 1e-12);

         // Tidy up
        for (unsigned i=0; i<nodes.size(); ++i)
        {
            delete nodes[i];
        }
        delete p_node;
    }


//     void xTestAnticlockwisenessOfNodes() throw(Exception)
//     {
//        // Tests to check that the nodes are anticlockwise when we create element
//        std::vector<Node<2>*> corner_nodes;
//        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
//
//        std::vector<Node<2>*> corner_nodes2;
//        corner_nodes2.push_back(new Node<2>(0, false, 0.0, 0.0));
//        corner_nodes2.push_back(new Node<2>(2, false, 1.0, 1.0));
//        corner_nodes2.push_back(new Node<2>(1, false, 1.0, 0.0));
//        corner_nodes2.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        VertexElement<2,2> vertex_element2(INDEX_IS_NOT_USED, corner_nodes2);
//
//        for (unsigned i=0; i<corner_nodes.size(); ++i)
//        {
//            delete corner_nodes[i];
//        }
//        for (unsigned i=0; i<corner_nodes2.size(); ++i)
//        {
//            delete corner_nodes2[i];
//        }
//     }


    void TestGetNodeLocalIndex()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;

        // This is a square
        nodes.push_back(new Node<2>(3, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(0, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);

        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0),3u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(1),2u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(2),1u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(3),0u);

        vertex_element.DeleteNode(3); // Removes (1,1) node

        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0), UINT_MAX);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

};
#endif /*TESTVERTEXELEMENT_HPP_*/
