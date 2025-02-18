/*

Copyright (c) 2005-2025, University of Oxford.
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


#ifndef _TESTNODE_HPP_
#define _TESTNODE_HPP_

#include <cxxtest/TestSuite.h>

#include <new>

#include "CheckpointArchiveTypes.hpp"

#include "Node.hpp"
#include "TetrahedralMesh.hpp"

#include "PetscSetupAndFinalize.hpp"

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
        TS_ASSERT_THROWS_THIS(node6.rGetNodeAttributes(), "Node has no attributes associated with it. Construct attributes first");
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

    void TestNodeWithAttributes()
    {
        Node<3> node(0, false, 0.0, 1.0, 2.0);

        TS_ASSERT(!node.HasNodeAttributes());
        TS_ASSERT_EQUALS(node.GetNumNodeAttributes(), 0u);

        // Region should default to 0 if not attributes are set up.
        TS_ASSERT_EQUALS(node.GetRegion(), 0u);

        TS_ASSERT_THROWS_THIS(node.CheckForNodeAttributes(), "Node has no attributes associated with it. Construct attributes first");

        node.ConstructNodeAttributes();
        TS_ASSERT(node.HasNodeAttributes());

        // Check defaults are all returned.
        TS_ASSERT_EQUALS(node.rGetNodeAttributes().size(), 0u);
        TS_ASSERT_EQUALS(node.GetRegion(), 0u);
        TS_ASSERT_DELTA(node.GetRadius(), 0.0, 1e-4);

        TS_ASSERT_DELTA(node.rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node.rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(node.rGetAppliedForce()[2], 0.0, 1e-4);

        TS_ASSERT_EQUALS(node.IsParticle(), false);

        TS_ASSERT_EQUALS(node.NeighboursIsEmpty(), true);

        // Check that we can correctly set each of the attributes.
        node.AddNodeAttribute(1.0);
        TS_ASSERT_EQUALS(node.GetNumNodeAttributes(), 1u);
        TS_ASSERT_EQUALS(node.rGetNodeAttributes().size(), 1u);
        TS_ASSERT_DELTA(node.rGetNodeAttributes()[0], 1.0, 1e-4);

        node.SetRegion(1u);
        TS_ASSERT_EQUALS(node.GetRegion(), 1u);

        c_vector<double, 3> force_contribution = scalar_vector<double>(3, 1.0);
        node.AddAppliedForceContribution(force_contribution);

        TS_ASSERT_DELTA(node.rGetAppliedForce()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(node.rGetAppliedForce()[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(node.rGetAppliedForce()[2], 1.0, 1e-4);

        node.ClearAppliedForce();
        TS_ASSERT_DELTA(node.rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node.rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(node.rGetAppliedForce()[2], 0.0, 1e-4);

        TS_ASSERT_EQUALS(node.IsParticle(), false);
        node.SetIsParticle(true);
        TS_ASSERT_EQUALS(node.IsParticle(), true);

        TS_ASSERT_EQUALS(node.NeighboursIsEmpty(), true);
        node.AddNeighbour(1u);
        node.AddNeighbour(1u);
        TS_ASSERT_EQUALS(node.NeighboursIsEmpty(), false);
        TS_ASSERT_EQUALS(node.rGetNeighbours().size(), 2u);
        node.RemoveDuplicateNeighbours();
        TS_ASSERT_EQUALS(node.rGetNeighbours().size(), 1u);
        node.ClearNeighbours();
        TS_ASSERT_EQUALS(node.NeighboursIsEmpty(), true);


        node.SetRadius(1.6);
        TS_ASSERT_DELTA(node.GetRadius(), 1.6, 1e-4);
    }

    void TestArchiveNode()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestNode", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "node.arch";

        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            Node<3>* const p_node = new Node<3>(0, true, 0.0, 1.0, 2.0);
            Node<1>* const p_node_1d = new Node<1>(0, false, 100.0);

//            p_node->AddNodeAttribute(5.0);
            p_node->SetRegion(7);

            // Write the nodes to file
            output_arch << p_node;
            output_arch << p_node_1d;
            delete p_node;
            delete p_node_1d;
        }

        {
            // Restore the nodes
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            Node<3>* p_node;
            input_arch >> p_node;
            TS_ASSERT_EQUALS(p_node->rGetLocation()[0], 0.0);
            TS_ASSERT_EQUALS(p_node->rGetLocation()[1], 1.0);
            TS_ASSERT_EQUALS(p_node->rGetLocation()[2], 2.0);
            TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_node->GetRegion(), 7u);

            Node<1>* p_node_1d;
            input_arch >> p_node_1d;
            TS_ASSERT_EQUALS(p_node_1d->rGetLocation()[0], 100.0);
            TS_ASSERT_EQUALS(p_node_1d->IsBoundaryNode(), false);
            TS_ASSERT_EQUALS(p_node_1d->GetRegion(), 0u);

//
//
//            TS_ASSERT_EQUALS(p_node->rGetNodeAttributes()[0], 5.0);
//
            delete p_node;
            delete p_node_1d;
        }
    }


    /*
     * The following tests relate to some spurious GCC warnings about potentially uninitialised c_vectors
     * in Release mode. See https://github.com/Chaste/Chaste/issues/231 for full details. These aim to
     * verify that the warnings are indeed false positives.
     *
     * Note that the warnings can be reproduced by:
     *  1. Compiling in Release mode, with GCC >= 9.4 (reproduced with GCC == 13.2.0 on 2024-10-02)
     *  2. Turning warnings back on by adding the following to the Chaste_ADD_TEST macro in ChasteMacros.cmake,
     *     after line 147:
     *       if ("${_testname}" STREQUAL "TestNode")
     *         set_source_files_properties("${_test_real_output_filename}" PROPERTIES COMPILE_FLAGS
     *           "-Wmaybe-uninitialized -Warray-bounds -Wstringop-overflow -Wstringop-overread")
     *       endif ()
     */

    void TestGccReleaseModeBehaviour1dSingleNode()
    {
        // When the spurious warnings are not ignored, this test should result in a compilation error along these lines:

        // inlined from ‘virtual void TestDescription_TestNode_TestGccReleaseModeBehaviour1dSingleNode::runTest()’ at Chaste/Chaste/build/mesh/test/TestNode.cpp:67:73:
        // Chaste/mesh/test/TestNode.hpp:350:9: error: ‘location_1d_from_node.boost::numeric::ublas::c_vector<double, 1>::data_[0]’ may be used uninitialized [-Werror=maybe-uninitialized]
        // 350 |         TS_ASSERT(location_1d_from_node[0] == 1.23);
        //     |         ^~~~~~~~~

        Node<1> my_node_1d(0u, false, 1.23);
        c_vector<double, 1> location_1d_from_node = my_node_1d.rGetLocation();

        TS_ASSERT(location_1d_from_node[0] == 1.23);
    }

    void TestGccReleaseModeBehaviour2dSingleNode()
    {
        // When the spurious warnings are not ignored, this test should result in a compilation error along these lines:

        // inlined from ‘virtual void TestDescription_TestNode_TestGccReleaseModeBehaviour2dSingleNode::runTest()’ at Chaste/Chaste/build/mesh/test/TestNode.cpp:67:73:
        // Chaste/Chaste/mesh/test/TestNode.hpp:372:9: error: ‘location_2d_from_node.boost::numeric::ublas::c_vector<double, 2>::data_[0]’ may be used uninitialized [-Werror=maybe-uninitialized]
        // 372 |         TS_ASSERT(location_2d_from_node[0] == 2.34);
        //     |         ^~~~~~~~~

        Node<2> my_node_2d(0u, false, 2.34, 3.45);
        c_vector<double, 2> location_2d_from_node = my_node_2d.rGetLocation();

        TS_ASSERT(location_2d_from_node[0] == 2.34);
        TS_ASSERT(location_2d_from_node[1] == 3.45);
    }

    void TestGccReleaseModeBehaviour1dMultipleNodes()
    {
        // This test relates to a class of warnings that indicated a built-in copy operation was copying too much data,
        // i.e. copying beyond the end of a c_vector. This test copy-constructs multiple c_vectors contiguously, to
        // verify that no data is getting overwritten. We used placement new in an attempt to ensure no default
        // initialization was occurring e.g. when creating a std::array or std::vector. This did not reproduce the
        // warning, though, but this test still demonstrates that copy construction of c_vectors is working as intended.

        Node<1> node_1(0u, false, 4.56);
        Node<1> node_2(1u, false, 5.67);
        Node<1> node_3(2u, false, 6.78);

        // Allocate some raw memory so we get three contiguous c_vectors, copy-constructed from the node locations
        void* raw_memory = operator new[](3 * sizeof(c_vector<double, 1>));
        c_vector<double, 1>* loc1 = new (raw_memory) c_vector<double, 1>(node_1.rGetLocation());
        c_vector<double, 1>* loc2 = new (loc1 + 1) c_vector<double, 1>(node_2.rGetLocation());
        c_vector<double, 1>* loc3 = new (loc2 + 1) c_vector<double, 1>(node_3.rGetLocation());

        // Test that none of these values has been clobbered by copy construction of the other elements
        TS_ASSERT((*loc1)[0] == 4.56);
        TS_ASSERT((*loc2)[0] == 5.67);
        TS_ASSERT((*loc3)[0] == 6.78);

        // Just for good measure, do some more copy construction and check again
        c_vector<double, 1> loc1_2 = *loc1;
        c_vector<double, 1> loc2_2 = *loc2;
        c_vector<double, 1> loc3_2 = *loc3;

        TS_ASSERT(loc1_2[0] == 4.56);
        TS_ASSERT(loc2_2[0] == 5.67);
        TS_ASSERT(loc3_2[0] == 6.78);

        // Tidy up memory
        loc1->~c_vector();
        loc2->~c_vector();
        loc3->~c_vector();
        operator delete[](raw_memory);
    }

    void TestGccReleaseModeBehaviour3dMultipleNodes()
    {
        // This test relates to a class of warnings that indicated a built-in copy operation was copying too much data,
        // i.e. copying beyond the end of a c_vector. This test copy-constructs multiple c_vectors contiguously, to
        // verify that no data is getting overwritten. We used placement new in an attempt to ensure no default
        // initialization was occurring e.g. when creating a std::array or std::vector. This did not reproduce the
        // warning, though, but this test still demonstrates that copy construction of c_vectors is working as intended.

        Node<3> node_1(0u, false, 12.3, 13.4, 14.5);
        Node<3> node_2(1u, false, 15.6, 16.7, 17.8);
        Node<3> node_3(2u, false, 21.2, 22.3, 23.4);

        // Allocate some raw memory so we get three contiguous c_vectors, copy-constructed from the node locations
        void* raw_memory = operator new[](3 * sizeof(c_vector<double, 3>));
        c_vector<double, 3>* loc1 = new (raw_memory) c_vector<double, 3>(node_1.rGetLocation());
        c_vector<double, 3>* loc2 = new (loc1 + 1) c_vector<double, 3>(node_2.rGetLocation());
        c_vector<double, 3>* loc3 = new (loc2 + 1) c_vector<double, 3>(node_3.rGetLocation());

        // Test that none of these values has been clobbered by copy construction of the other elements
        TS_ASSERT((*loc1)[0] == 12.3);
        TS_ASSERT((*loc1)[1] == 13.4);
        TS_ASSERT((*loc1)[2] == 14.5);

        TS_ASSERT((*loc2)[0] == 15.6);
        TS_ASSERT((*loc2)[1] == 16.7);
        TS_ASSERT((*loc2)[2] == 17.8);

        TS_ASSERT((*loc3)[0] == 21.2);
        TS_ASSERT((*loc3)[1] == 22.3);
        TS_ASSERT((*loc3)[2] == 23.4);

        // Just for good measure, do some more copy construction and check again
        c_vector<double, 3> loc1_2 = *loc1;
        c_vector<double, 3> loc2_2 = *loc2;
        c_vector<double, 3> loc3_2 = *loc3;

        TS_ASSERT(loc1_2[0] == 12.3);
        TS_ASSERT(loc1_2[1] == 13.4);
        TS_ASSERT(loc1_2[2] == 14.5);

        TS_ASSERT(loc2_2[0] == 15.6);
        TS_ASSERT(loc2_2[1] == 16.7);
        TS_ASSERT(loc2_2[2] == 17.8);

        TS_ASSERT(loc3_2[0] == 21.2);
        TS_ASSERT(loc3_2[1] == 22.3);
        TS_ASSERT(loc3_2[2] == 23.4);

        // Tidy up memory
        loc1->~c_vector();
        loc2->~c_vector();
        loc3->~c_vector();
        operator delete[](raw_memory);
    }
};

#endif //_TESTNODE_HPP_
