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

#ifndef TESTPOTTSMESHGENERATOR_HPP_
#define TESTPOTTSMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "PottsMeshGenerator.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestPottsMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestPottsMeshGeneratorIn2dFromBottomLeft()
    {
        // Coverage
        //TS_ASSERT_THROWS_NOTHING(PottsMeshGenerator<2> empty_generator());
    ///\todo Line above is not valid code.  What does it cover?

        PottsMeshGenerator<2> generator(9, 3, 3, 5, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left

        // Create mesh
        PottsMesh<2>* p_mesh = generator.GetMesh();

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

    void TestPottsMeshGenerator2dInCentre()
    {
        PottsMeshGenerator<2> generator(6, 2, 2, 7, 2, 2); //should have a gap of one on the left right and bottom and 2 on the top

        // Create mesh
        PottsMesh<2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 42u);

        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[0], 5.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[1], 0.0, 1e-3);

        TS_ASSERT_DELTA(p_mesh->GetNode(15)->rGetLocation()[0], 3.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(15)->rGetLocation()[1], 2.0, 1e-3);

        // Check that each node on the left right and top aren't in any elements and the rest
        // are in only one element
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            if ((node_index < 6) || (node_index > 29) || (node_index % 6 == 0) || (node_index % 6 == 5))
            {
                TS_ASSERT_EQUALS(num_containing_elements,0u);
            }
            else
            {
                TS_ASSERT_EQUALS(num_containing_elements,1u);
            }
        }

        // Check that some elements contain the correct nodes

        // Element 0 contains nodes 7, 8, 13, and 14,
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(0), 7u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(1), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(2), 13u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(3), 14u);

        // Element 1 contains nodes 9, 10, 15,  and 16
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(0), 9u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(1), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(2), 15u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(3), 16u);

        // Element 2 contains nodes 19, 20, 25, and 26
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(0), 19u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(1), 20u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(2), 25u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(3), 26u);

        // Element 3 contains nodes 21, 22, 27, and 28
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(0), 21u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(1), 22u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(2), 27u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(3)->GetNodeGlobalIndex(3), 28u);

        // Check that some nodes are contained in the correct elements

        // Node 7 is only in element 0
        std::set<unsigned> temp_list_1;
        temp_list_1.insert(0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(7)->rGetContainingElementIndices(), temp_list_1);

        // Node 15 is in element 1
        std::set<unsigned> temp_list_2;
        temp_list_2.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(15)->rGetContainingElementIndices(), temp_list_2);

        // Node 30 is in no element
        std::set<unsigned> temp_list_3;
        TS_ASSERT_EQUALS(p_mesh->GetNode(30)->rGetContainingElementIndices(), temp_list_3);

        // The edges of the mesh should be boundary nodes
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            bool is_boundary_node = (node_index <= 5 || node_index>=36 || node_index % 6 == 0 || (node_index % 6 == 5)) ? true : false;
            TS_ASSERT_EQUALS(p_mesh->GetNode(node_index)->IsBoundaryNode(), is_boundary_node);
        }
    }

    void TestPottsMeshGenerator3dFromBottomLeft()
    {
        PottsMeshGenerator<3> generator(4, 2, 2, 4, 2, 2, 6, 2, 2, true); // last bool makes elements start in bottom left

        // Create mesh
        PottsMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 96u);

        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[0], 1.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[1], 1.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(p_mesh->GetNode(16)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(16)->rGetLocation()[1], 0.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(16)->rGetLocation()[2], 1.0, 1e-3);

        // Check that each node under 36 is contained in only one element and the rest aren't
        // in any elements
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            if (node_index < 64)
            {
                TS_ASSERT_EQUALS(num_containing_elements,1u);
            }
            else
            {
                TS_ASSERT_EQUALS(num_containing_elements,0u);
            }
        }

        // Check that some elements contain the correct nodes

        // Element 0 contains nodes 0, 1, 4, 5, 16, 17, 20, 21
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(2), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(3), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(4), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(5), 17u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(6), 20u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(7), 21u);

        // Element 2 contains nodes 8, 9, 12 ,13, 24 ,25, 28, 29
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(0), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(1), 9u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(2), 12u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(3), 13u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(4), 24u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(5), 25u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(6), 28u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(7), 29u);

        // Element 4 contains nodes 32, 33, 36, 37, 48, 49, 52, 53
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(0), 32u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(1), 33u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(2), 36u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(3), 37u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(4), 48u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(5), 49u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(6), 52u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(4)->GetNodeGlobalIndex(7), 53u);

        // Check that some nodes are contained in the correct elements

        // Node 0 is only in element 0
        std::set<unsigned> temp_list_1;
        temp_list_1.insert(0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(0)->rGetContainingElementIndices(), temp_list_1);

        // Node 13 is in element 2
        std::set<unsigned> temp_list_2;
        temp_list_2.insert(2u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(13)->rGetContainingElementIndices(), temp_list_2);

        // Node 65 is in no element
        std::set<unsigned> temp_list_3;
        TS_ASSERT_EQUALS(p_mesh->GetNode(65)->rGetContainingElementIndices(), temp_list_3);

        // The edges of the mesh should be boundary nodes
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            bool is_boundary_node = (node_index%16==0 || node_index%16==1 || node_index%16==2 || node_index%16==3 ||
                                     node_index%16==4 || node_index%16==7 || node_index==9 || node_index==10 || node_index==5 || node_index==6 ||
                                     node_index%16==8 || node_index%16==11 || node_index==85 || node_index==86 || node_index==89 || node_index==90 ||
                                     node_index%16==12 || node_index%16==13 || node_index%16==14 || node_index%16==15) ? true : false;
            TS_ASSERT_EQUALS(p_mesh->GetNode(node_index)->IsBoundaryNode(), is_boundary_node);
        }
    }

    void TestPottsMeshGenerator3dInCentre()
    {
        PottsMeshGenerator<3> generator(6, 2, 2, 4, 1, 2, 4, 1, 2); //should have a gap of one on all sides

        // Create mesh
        PottsMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 96u);

        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[0], 5.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[1], 0.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(p_mesh->GetNode(16)->rGetLocation()[0], 4.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(16)->rGetLocation()[1], 2.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(16)->rGetLocation()[2], 0.0, 1e-3);

        // Check that each node under 36 is contained in only one element and the rest aren't
        // in any elements
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            if ((p_mesh->GetNode(node_index)->rGetLocation()[0] <= 0.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[0] >= 5.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[1] <= 0.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[1] >= 3.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[2] <= 0.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[2] >= 3.0))
            {
                TS_ASSERT_EQUALS(num_containing_elements,0u);
            }
            else
            {
                TS_ASSERT_EQUALS(num_containing_elements,1u);
            }
        }

        // Check that the elements contain the correct nodes

        // Element 0 contains nodes 0, 1, 4, 5, 16, 17, 20, 21
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(0), 31u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(1), 32u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(2), 37u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(3), 38u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(4), 55u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(5), 56u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(6), 61u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNodeGlobalIndex(7), 62u);

        // Element 2 contains nodes 8, 9, 12 ,13, 24 ,25, 28, 29
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(0), 33u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(1), 34u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(2), 39u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(3), 40u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(4), 57u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(5), 58u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(6), 63u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(1)->GetNodeGlobalIndex(7), 64u);

        // Check that some nodes are contained in the correct elements

        // Node 31 is only in element 0
        std::set<unsigned> temp_list_1;
        temp_list_1.insert(0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(31)->rGetContainingElementIndices(), temp_list_1);

        // Node 34 is in element 1
        std::set<unsigned> temp_list_2;
        temp_list_2.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(34)->rGetContainingElementIndices(), temp_list_2);

        // Node 72 is in no element
        std::set<unsigned> temp_list_3;
        TS_ASSERT_EQUALS(p_mesh->GetNode(72)->rGetContainingElementIndices(), temp_list_3);

        // The edges of the mesh should be boundary nodes
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            bool is_boundary_node = ( (p_mesh->GetNode(node_index)->rGetLocation()[0] <= 0.0) ||
                                      (p_mesh->GetNode(node_index)->rGetLocation()[0] >= 5.0) ||
                                      (p_mesh->GetNode(node_index)->rGetLocation()[1] <= 0.0) ||
                                      (p_mesh->GetNode(node_index)->rGetLocation()[1] >= 3.0) ||
                                      (p_mesh->GetNode(node_index)->rGetLocation()[2] <= 0.0) ||
                                      (p_mesh->GetNode(node_index)->rGetLocation()[2] >= 3.0) ) ? true : false;

            TS_ASSERT_EQUALS(p_mesh->GetNode(node_index)->IsBoundaryNode(), is_boundary_node);
        }
    }

    void TestGenerator3dLarge()
    {
        // Create a simple 3D PottsMesh
        unsigned domain_size = 10;
        unsigned element_number = 4;
        unsigned element_size = 2;

        PottsMeshGenerator<3> generator(domain_size, element_number, element_size, domain_size, element_number, element_size, domain_size, element_number, element_size);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 64u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 1000u);

        // Check that each node on the boundary is in no element and the rest
        // are in one elements.
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            bool is_boundary_node;

            if ((p_mesh->GetNode(node_index)->rGetLocation()[0] <= 0.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[0] >= domain_size - 1.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[1] <= 0.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[1] >= domain_size - 1.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[2] <= 0.0) ||
                (p_mesh->GetNode(node_index)->rGetLocation()[2] >= domain_size - 1.0))
            {
                TS_ASSERT_EQUALS(num_containing_elements,0u);
                is_boundary_node = true;
            }
            else
            {
                TS_ASSERT_EQUALS(num_containing_elements,1u);
                is_boundary_node = false;
            }
            TS_ASSERT_EQUALS(p_mesh->GetNode(node_index)->IsBoundaryNode(), is_boundary_node);
        }
    }
};

#endif /*TESTPOTTSMESHGENERATOR_HPP_*/
