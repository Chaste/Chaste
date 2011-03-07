/*

Copyright (C) University of Oxford, 2005-2011

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

    void Test2dPottsElement()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.0, 1.0));

        PottsElement element(0, nodes);

        // Test RegisterWithNodes()
        element.RegisterWithNodes();
        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 1u);
        }

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), 1u);

        // Test IsElementOnBoundary()
        TS_ASSERT_EQUALS(element.IsElementOnBoundary(),true);

        // Test UpdateNode()
        Node<2>* p_node_2 = new Node<2>(2, false, 1.0, 0.0);
        element.UpdateNode(0, p_node_2);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 0.0, 1e-12);

        // Test ResetIndex()
        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
        element.ResetIndex(5);
        TS_ASSERT_EQUALS(element.GetIndex(), 5u);

        // Test DeleteNode() and AddNode()
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        element.DeleteNode(1);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 1u);

        Node<2>* p_node_3 = new Node<2>(3, false, 1.0, 1.0);
        element.AddNode(p_node_3);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[1], 1.0, 1e-12);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(2), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(3), 1u);

        // Test MarkAsDeleted()
        element.MarkAsDeleted();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_node_2;
        delete p_node_3;
    }
};

#endif /*TESTPOTTSELEMENT_HPP_*/
