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


#ifndef _TESTNODE_HPP_
#define _TESTNODE_HPP_

#include "Node.hpp"
#include "TetrahedralMesh.hpp"
#include <cxxtest/TestSuite.h>

class TestNode : public CxxTest::TestSuite
{
public:

    /**
     * Test that get and set methods work as expected.
     * Also check constructors.
     * We only test 1 dimensional nodes. Nothing much changes in higher
     * dimensions.
     */
    void TestNodeMethod()
    {
        ChastePoint<1> point1(1.0);
        ChastePoint<1> point2(2.0);

        Node<1> node1(0, point1);
        TS_ASSERT_EQUALS(node1.GetIndex(), 0u);
        TS_ASSERT_DELTA(1.0, node1.GetPoint()[0], 1e-12);

        node1.SetIndex(1);
        TS_ASSERT_EQUALS(node1.GetIndex(), 1u);

        node1.SetPoint(point2);
        TS_ASSERT_DELTA(2.0, node1.GetPoint()[0], 1e-12);

        node1.SetAsBoundaryNode();
        TS_ASSERT(node1.IsBoundaryNode());

        node1.SetAsBoundaryNode(false);
        TS_ASSERT(!node1.IsBoundaryNode());

        Node<2> node2(1, false, 1.0, 2.0);
        TS_ASSERT_EQUALS(node2.GetIndex(), 1u);
        TS_ASSERT_DELTA(node2.GetPoint()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(node2.GetPoint()[1], 2.0, 1e-12);
        TS_ASSERT_EQUALS(node2.IsBoundaryNode(), false);

        // This shouldn't give an error, even though we specify too many
        // coordinates: the 3rd coord should be ignored.
        Node<2> node3(2, true, 1.0, 2.0, 3.0);

        double location[3];
        location[0]=100.0;
        location[1]=101.0;
        location[2]=102.0;
        Node<1> node4(2, location, true);
        Node<2> node5(2, location, true);
        Node<3> node6(2, location, true);
        TS_ASSERT_EQUALS(node4.IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(node5.IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(node6.IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(node4.GetIndex(), 2u);
        TS_ASSERT_EQUALS(node5.GetIndex(), 2u);
        TS_ASSERT_EQUALS(node6.GetIndex(), 2u);
        TS_ASSERT_DELTA(node4.GetPoint()[0], location[0], 1e-12);
        TS_ASSERT_DELTA(node5.GetPoint()[0], location[0], 1e-12);
        TS_ASSERT_DELTA(node6.GetPoint()[0], location[0], 1e-12);

        TS_ASSERT_DELTA(node5.GetPoint()[1], location[1], 1e-12);
        TS_ASSERT_DELTA(node6.GetPoint()[1], location[1], 1e-12);

        TS_ASSERT_DELTA(node6.GetPoint()[2], location[2], 1e-12);

        // Test the node attributes
        TS_ASSERT_EQUALS(node6.rGetNodeAttributes().size(),0u);
        double attribute = 54.98;
        node6.AddNodeAttribute(attribute);
        TS_ASSERT_EQUALS(node6.rGetNodeAttributes().size(),1u);
        TS_ASSERT_DELTA(node6.rGetNodeAttributes()[0],attribute, 1e-12);

        // Add another one
        node6.AddNodeAttribute(attribute*2);
        TS_ASSERT_EQUALS(node6.rGetNodeAttributes().size(),2u);
        TS_ASSERT_DELTA(node6.rGetNodeAttributes()[1],attribute*2, 1e-12);

        // Test node deletion (from a mesh) methods
        TS_ASSERT_EQUALS(node1.IsDeleted(), false);
        node1.MarkAsDeleted();
        TS_ASSERT_EQUALS(node1.IsDeleted(), true);

        TS_ASSERT_EQUALS(node1.GetRegion(), 0u);
        node1.SetRegion(4);
        TS_ASSERT_EQUALS(node1.GetRegion(), 4u);

        TS_ASSERT_EQUALS(node1.IsInternal(), false);
        node1.MarkAsInternal();
        TS_ASSERT_EQUALS(node1.IsInternal(), true);
    }

    /**
     * Test that we can use both the new and old interfaces
     */
    void TestNodeNewAndOld()
    {
        ChastePoint<1> point1(1.0);

        // Create a node with old interface
        Node<1> node1(0, point1);

        // check location with new interface
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], point1[0]);

        c_vector<double, 1> new_location;
        new_location[0] = 10.0;

        // Update location with new interface
        node1.rGetModifiableLocation() = new_location;

        // Check location with old interface
        TS_ASSERT_EQUALS(node1.GetPoint()[0], 10.0);

        // Update location with old interface
        node1.SetPoint(point1);

        // Check location with new interface
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], point1[0]);
    }

    void TestNodeNew()
    {
        c_vector<double, 1> point3;
        point3[0] = 23.0;

        // Create node
        Node<1> node1(0, point3);

        // Read back location
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], point3[0]);

        c_vector<double, 1> new_location;
        new_location[0] = 10.0;

        // Update location
        node1.rGetModifiableLocation() = new_location;

        // Read back location
        TS_ASSERT_EQUALS(node1.rGetLocation()[0], 10.0);
    }

    void TestNodeFromStdVector()
    {
        std::vector<double> coords(3u);
        coords[0] = 1.5;
        coords[1] = 15.9;
        coords[2] = 777.7;
        Node<3> node(0, coords);

        for (int i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(node.rGetLocation()[i], coords[i], 1e-12);
        }

        TS_ASSERT_THROWS_THIS(node.RemoveElement(256),"Tried to remove an index which was not in the set");
        TS_ASSERT_THROWS_THIS(node.RemoveBoundaryElement(256),"Tried to remove an index which was not in the set");
    }

    void TestFlagging()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1, 1);
        unsigned num_nodes = 4;

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsFlagged(mesh), false);
        }

        mesh.GetElement(0)->Flag();

        TS_ASSERT_EQUALS(mesh.GetNode(0)->IsFlagged(mesh), false);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->IsFlagged(mesh), true);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->IsFlagged(mesh), true);
        TS_ASSERT_EQUALS(mesh.GetNode(3)->IsFlagged(mesh), true);

        mesh.GetElement(0)->Unflag();

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetNode(i)->IsFlagged(mesh), false);
        }

        mesh.GetElement(1)->Flag();
        TS_ASSERT_EQUALS(mesh.GetNode(0)->IsFlagged(mesh), true);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->IsFlagged(mesh), true);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->IsFlagged(mesh), true);
        TS_ASSERT_EQUALS(mesh.GetNode(3)->IsFlagged(mesh), false);
    }

};

#endif //_TESTNODE_HPP_
