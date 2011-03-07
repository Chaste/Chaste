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

#ifndef TESTPOTTSMESHGENERATOR_HPP_
#define TESTPOTTSMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "PottsMeshGenerator.hpp"

class TestPottsMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestPottsMeshGeneratorSimple() throw(Exception)
    {
        PottsMeshGenerator generator(9, 5, 3, 2, 3, 2);

        // Create mesh
        PottsMesh* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 45u);

        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[0], 5.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[1], 0.0, 1e-3);

        TS_ASSERT_DELTA(p_mesh->GetNode(15)->rGetLocation()[0], 6.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(15)->rGetLocation()[1], 1.0, 1e-3);

        // Check that each node under 36 is contained in only one element and the rest aren't
        // in any elements
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            if (node_index < 36)
            {
                TS_ASSERT_EQUALS(num_containing_elements,1u);
            }
            else
            {
                TS_ASSERT_EQUALS(num_containing_elements,0u);
            }
        }

        // Check that some elements contain the correct nodes

        // Element 0 contains nodes 0, 1, 2, 9, 10, and 11
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(3), 9u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(4), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(5), 11u);

        // Element 2 contains nodes 6, 7, 8,  and 7
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(0), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(1), 7u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(2), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(3), 15u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(4), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(5), 17u);

        // Element 4 contains nodes 19, 20, 21, 28, 29, and 30
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(0), 21u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(1), 22u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(2), 23u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(3), 30u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(4), 31u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(5), 32u);

        // Check that some nodes are contained in the correct elements

        // Node 0 is only in element 0
        std::set<unsigned> temp_list_1;
        temp_list_1.insert(0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(0)->rGetContainingElementIndices(), temp_list_1);

        // Node 13 is in element 1
        std::set<unsigned> temp_list_2;
        temp_list_2.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(13)->rGetContainingElementIndices(), temp_list_2);

        // Node 38 is in no element
        std::set<unsigned> temp_list_3;
        TS_ASSERT_EQUALS(p_mesh->GetNode(38)->rGetContainingElementIndices(), temp_list_3);

        // The edges of the mesh should be boundary nodes
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            bool is_boundary_node = (node_index <= 9 || node_index==17 || node_index==18 || node_index==26 || node_index==27 || node_index>=35) ? true : false;
            TS_ASSERT_EQUALS(p_mesh->GetNode(node_index)->IsBoundaryNode(), is_boundary_node);
        }
    }
};

#endif /*TESTPOTTSMESHGENERATOR_HPP_*/
