/*

Copyright (c) 2005-2014, University of Oxford.
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
#ifndef TESTBOXCOLLECTION_HPP_
#define TESTBOXCOLLECTION_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "TetrahedralMesh.hpp"
#include "BoxCollection.hpp"
#include "TrianglesMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestBoxCollection : public CxxTest::TestSuite
{
public:
    void TestBox() throw (Exception)
    {
        c_vector<double, 2*2> box_size;
        box_size(0) = -0.1; // min x
        box_size(1) = 1.1; // max x
        box_size(2) = -0.1; // min y
        box_size(3) = 1.1; // max y

        Box<2> test_box(box_size);
        c_vector<double, 2*2> returned_min_max_values = test_box.rGetMinAndMaxValues();
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(returned_min_max_values(i), box_size(i));
        }

        c_vector<double, 2> node_location;
        node_location(0) = 0.5;
        node_location(1) = 0.5;

        Node<2> test_node(213, node_location);

        test_box.AddNode(&test_node);
        std::set< Node<2>* > nodes_contained_before = test_box.rGetNodesContained();

        TS_ASSERT_EQUALS(*(nodes_contained_before.begin()), &test_node);
        TS_ASSERT_EQUALS((*(nodes_contained_before.begin()))->GetIndex(), 213u);

        test_box.RemoveNode(&test_node);
        std::set< Node<2>* > nodes_contained_after = test_box.rGetNodesContained();
        TS_ASSERT(nodes_contained_after.empty());
    }


    void TestBoxGeneration1d() throw (Exception)
    {
        // Create a mesh
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(20);

        double cut_off_length = 5.0;
        c_vector<double, 2> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 20.15;

        BoxCollection<1> box_collection(cut_off_length, domain_size);

        box_collection.SetupAllLocalBoxes();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 5u);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            std::set< Node<1>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
            c_vector<double, 2> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

            for (std::set< Node<1>* >::iterator it_nodes_in_box = nodes_in_box.begin();
                 it_nodes_in_box != nodes_in_box.end();
                 it_nodes_in_box++)
            {
                Node<1>* current_node = *it_nodes_in_box;
                double x_position = current_node->rGetLocation()[0];

                double epsilon = 1e-12;

                TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
                TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
            }
        }

        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_1 = box_collection.GetLocalBoxes(1);
        std::set<unsigned> correct_answer_1;
        correct_answer_1.insert(0);
        correct_answer_1.insert(1);
        correct_answer_1.insert(2);
        TS_ASSERT_EQUALS(local_boxes_to_box_1, correct_answer_1);

        std::set<unsigned> local_boxes_to_box_4 = box_collection.GetLocalBoxes(4);
        std::set<unsigned> correct_answer_4;
        correct_answer_4.insert(3);
        correct_answer_4.insert(4);
        TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);

        c_vector<double,1> miles_away;
        miles_away(0) = 47323854;
        TS_ASSERT_THROWS_CONTAINS(box_collection.CalculateContainingBox(miles_away), "The point provided is outside all of the boxes");
    }


    // very simple test
    void TestAddElement() throw(Exception)
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.5, 1.0);

        double width = 0.4;
        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;

        BoxCollection<1> box_collection(width, domain_size);
        box_collection.rGetBox(0).AddElement(mesh.GetElement(0));
        TS_ASSERT_EQUALS(box_collection.rGetBox(0).rGetElementsContained().size(), 1u);
        TS_ASSERT_EQUALS(box_collection.rGetBox(1).rGetElementsContained().size(), 0u);
        TS_ASSERT_EQUALS(box_collection.rGetBox(2).rGetElementsContained().size(), 0u);
        TS_ASSERT_EQUALS(*(box_collection.rGetBox(0).rGetElementsContained().begin()), mesh.GetElement(0));
    }


    void TestSetupAllLocalBoxes2d() throw(Exception)
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;

        BoxCollection<2> box_collection(width, domain_size);

        assert(box_collection.GetNumBoxes()==12); // 4 * 3 boxes altogether

        box_collection.SetupAllLocalBoxes();

        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(4);
        correct_answer_0.insert(5);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.GetLocalBoxes(3);
        std::set<unsigned> correct_answer_3;
        correct_answer_3.insert(3);
        correct_answer_3.insert(2);
        correct_answer_3.insert(6);
        correct_answer_3.insert(7);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.GetLocalBoxes(5);
        std::set<unsigned> correct_answer_5;
        correct_answer_5.insert(0);
        correct_answer_5.insert(1);
        correct_answer_5.insert(2);
        correct_answer_5.insert(4);
        correct_answer_5.insert(5);
        correct_answer_5.insert(6);
        correct_answer_5.insert(8);
        correct_answer_5.insert(9);
        correct_answer_5.insert(10);
        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_10 = box_collection.GetLocalBoxes(10);
        std::set<unsigned> correct_answer_10;
        correct_answer_10.insert(5);
        correct_answer_10.insert(6);
        correct_answer_10.insert(7);
        correct_answer_10.insert(9);
        correct_answer_10.insert(10);
        correct_answer_10.insert(11);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);
    }


    void TestSetupAllLocalBoxes2dPeriodic() throw(Exception)
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;

        BoxCollection<2> box_collection(width, domain_size, true); // So periodic in X

        assert(box_collection.GetNumBoxes()==12); // 4 * 3 boxes altogether

        box_collection.SetupAllLocalBoxes();

        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(3);
        correct_answer_0.insert(4);
        correct_answer_0.insert(5);
        correct_answer_0.insert(7);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.GetLocalBoxes(3);
        std::set<unsigned> correct_answer_3;
        correct_answer_3.insert(0);
        correct_answer_3.insert(2);
        correct_answer_3.insert(3);
        correct_answer_3.insert(4);
        correct_answer_3.insert(6);
        correct_answer_3.insert(7);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.GetLocalBoxes(5);
        std::set<unsigned> correct_answer_5;
        correct_answer_5.insert(0);
        correct_answer_5.insert(1);
        correct_answer_5.insert(2);
        correct_answer_5.insert(4);
        correct_answer_5.insert(5);
        correct_answer_5.insert(6);
        correct_answer_5.insert(8);
        correct_answer_5.insert(9);
        correct_answer_5.insert(10);
        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_10 = box_collection.GetLocalBoxes(10);
        std::set<unsigned> correct_answer_10;
        correct_answer_10.insert(5);
        correct_answer_10.insert(6);
        correct_answer_10.insert(7);
        correct_answer_10.insert(9);
        correct_answer_10.insert(10);
        correct_answer_10.insert(11);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);

        std::set<unsigned> local_boxes_to_box_11 = box_collection.GetLocalBoxes(11);
        std::set<unsigned> correct_answer_11;
        correct_answer_11.insert(4);
        correct_answer_11.insert(6);
        correct_answer_11.insert(7);
        correct_answer_11.insert(8);
        correct_answer_11.insert(10);
        correct_answer_11.insert(11);
        TS_ASSERT_EQUALS(local_boxes_to_box_11, correct_answer_11);
    }


    void TestSetupAllLocalBoxes3d() throw(Exception)
    {
        double width = 1.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;
        domain_size(4) = 0;
        domain_size(5) = 2-0.01;

        BoxCollection<3> box_collection(width, domain_size);

        assert(box_collection.GetNumBoxes()==24); // 4 * 3 * 2 boxes altogether

        box_collection.SetupAllLocalBoxes();

        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(4);
        correct_answer_0.insert(5);
        correct_answer_0.insert(12);
        correct_answer_0.insert(13);
        correct_answer_0.insert(16);
        correct_answer_0.insert(17);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.GetLocalBoxes(3);
        std::set<unsigned> correct_answer_3;
        correct_answer_3.insert(3);
        correct_answer_3.insert(2);
        correct_answer_3.insert(6);
        correct_answer_3.insert(7);
        correct_answer_3.insert(14);
        correct_answer_3.insert(15);
        correct_answer_3.insert(18);
        correct_answer_3.insert(19);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.GetLocalBoxes(5);
        std::set<unsigned> correct_answer_5;
        correct_answer_5.insert(0);
        correct_answer_5.insert(1);
        correct_answer_5.insert(2);
        correct_answer_5.insert(4);
        correct_answer_5.insert(5);
        correct_answer_5.insert(6);
        correct_answer_5.insert(8);
        correct_answer_5.insert(9);
        correct_answer_5.insert(10);
        correct_answer_5.insert(12);
        correct_answer_5.insert(13);
        correct_answer_5.insert(14);
        correct_answer_5.insert(16);
        correct_answer_5.insert(17);
        correct_answer_5.insert(18);
        correct_answer_5.insert(20);
        correct_answer_5.insert(21);
        correct_answer_5.insert(22);

        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_19 = box_collection.GetLocalBoxes(19);
        std::set<unsigned> correct_answer_19;
        correct_answer_19.insert(2);
        correct_answer_19.insert(3);
        correct_answer_19.insert(6);
        correct_answer_19.insert(7);
        correct_answer_19.insert(10);
        correct_answer_19.insert(11);
        correct_answer_19.insert(14);
        correct_answer_19.insert(15);
        correct_answer_19.insert(18);
        correct_answer_19.insert(19);
        correct_answer_19.insert(22);
        correct_answer_19.insert(23);
        TS_ASSERT_EQUALS(local_boxes_to_box_19, correct_answer_19);

        std::set<unsigned> local_boxes_to_box_22 = box_collection.GetLocalBoxes(22);
        std::set<unsigned> correct_answer_22;
        correct_answer_22.insert(5);
        correct_answer_22.insert(6);
        correct_answer_22.insert(7);
        correct_answer_22.insert(9);
        correct_answer_22.insert(10);
        correct_answer_22.insert(11);
        correct_answer_22.insert(17);
        correct_answer_22.insert(18);
        correct_answer_22.insert(19);
        correct_answer_22.insert(21);
        correct_answer_22.insert(22);
        correct_answer_22.insert(23);
        TS_ASSERT_EQUALS(local_boxes_to_box_22, correct_answer_22);

    }

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    //
    //  Cancer/cell_based tests
    //  The following are tests written from this BoxCollection used to be
    //  cell_based/src/tissue/NodeBoxCollection and test the cell-based
    //  functionality
    //
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    void TestPairsReturned1d() throw (Exception)
    {
        std::vector< ChastePoint<1>* > points(5);
        points[0] = new ChastePoint<1>(0.2);
        points[1] = new ChastePoint<1>(0.7);
        points[2] = new ChastePoint<1>(1.3);
        points[3] = new ChastePoint<1>(2.1);
        points[4] = new ChastePoint<1>(6.9);

        std::vector<Node<1>* > nodes;
        for (unsigned i=0; i<points.size(); i++)
        {
            nodes.push_back(new Node<1>(i, *(points[i]), false));
        }

        double cut_off_length = 1.0;

        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 7.0;

        BoxCollection<1> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        std::vector< std::pair<Node<1>*, Node<1>* > > pairs_returned_vector;
        std::map<unsigned, std::set<unsigned> > neighbours_returned;
        box_collection.CalculateNodePairs(nodes,pairs_returned_vector, neighbours_returned);

        std::set< std::pair<Node<1>*, Node<1>* > > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            pairs_returned.insert(pairs_returned_vector[i]);
        }
        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
        neighbours_should_be[0].insert(1);
        neighbours_should_be[0].insert(2);
        neighbours_should_be[1].insert(0);
        neighbours_should_be[1].insert(2);
        neighbours_should_be[2].insert(0);
        neighbours_should_be[2].insert(1);
        neighbours_should_be[2].insert(3);
        neighbours_should_be[3].insert(2);
        neighbours_should_be[4] = std::set<unsigned>();

        TS_ASSERT_EQUALS(neighbours_should_be, neighbours_returned);

        std::set< std::pair<Node<1>*, Node<1>* > > pairs_should_be;
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[2]));
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));

        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

        // Check we empty boxes correctly
        box_collection.EmptyBoxes();
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            TS_ASSERT_EQUALS(box_collection.rGetBox(i).rGetNodesContained().size(), 0u);
        }

        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestBoxGeneration2d() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double cut_off_length = 0.2;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 1.15;
        domain_size(2) = -0.1;
        domain_size(3) = 1.15;

        BoxCollection<2> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 49u);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            std::set< Node<2>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
            c_vector<double, 2*2> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

            for (std::set< Node<2>* >::iterator it_nodes_in_box = nodes_in_box.begin();
                 it_nodes_in_box != nodes_in_box.end();
                 it_nodes_in_box++)
            {
                Node<2>* current_node = *it_nodes_in_box;
                double x_position = current_node->rGetLocation()[0];
                double y_position = current_node->rGetLocation()[1];

                double epsilon = 1e-12;

                TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
                TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
                TS_ASSERT_LESS_THAN(box_min_max_values(2)-epsilon, y_position);
                TS_ASSERT_LESS_THAN(y_position, box_min_max_values(3)+epsilon);
            }
        }

        // Have checked that all the local boxes are calculated correctly on a 5 by 6 grid - here we
        // hardcode a few checks on the 7 by 7 grid.
        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(7);
        correct_answer_0.insert(8);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_4 = box_collection.GetLocalBoxes(4);
        std::set<unsigned> correct_answer_4;
        correct_answer_4.insert(4);
        correct_answer_4.insert(5);
        correct_answer_4.insert(10);
        correct_answer_4.insert(11);
        correct_answer_4.insert(12);
        TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);

        std::set<unsigned> local_boxes_to_box_10 = box_collection.GetLocalBoxes(10);
        std::set<unsigned> correct_answer_10;
        correct_answer_10.insert(10);
        correct_answer_10.insert(11);
        correct_answer_10.insert(16);
        correct_answer_10.insert(17);
        correct_answer_10.insert(18);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);

        std::set<unsigned> local_boxes_to_box_48 = box_collection.GetLocalBoxes(48);
        std::set<unsigned> correct_answer_48;
        correct_answer_48.insert(48);
        TS_ASSERT_EQUALS(local_boxes_to_box_48, correct_answer_48);
    }

    /* This test insert to verify repeatability of BoxCollection floating point
     * calculations.  Note that failure of this test on a given architecture implies
     * the failure of node-based cell simulations
     */
    void TestLargeBoxCollection2d() throw (Exception)
    {
        double cut_off_length = 1e-3;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;
        domain_size(2) = 0.0;
        domain_size(3) = 1.0;

        BoxCollection<2> box_collection(cut_off_length, domain_size);
        TS_ASSERT_EQUALS(box_collection.mNumBoxesEachDirection[0], 1001u);
        TS_ASSERT_EQUALS(box_collection.mNumBoxesEachDirection[1], 1001u);

        c_vector<double, 2> probe;

        probe(0)=0.0; probe(1)=0.0;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 0u);

        probe(0)=1.0; probe(1)=0.0;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 1000u);

        probe(0)=0.0; probe(1)=1.0;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 1001000u);

        probe(0)=1.0; probe(1)=1.0;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 1002000u);

    }

    void TestPairsReturned2d() throw (Exception)
    {
        std::vector< ChastePoint<2>* > points(10);
        points[0] = new ChastePoint<2>(0.2, 3.7);
        points[1] = new ChastePoint<2>(0.5, 3.2);
        points[2] = new ChastePoint<2>(1.1, 1.99);
        points[3] = new ChastePoint<2>(1.3, 0.8);
        points[4] = new ChastePoint<2>(1.3, 0.3);
        points[5] = new ChastePoint<2>(2.2, 0.6);
        points[6] = new ChastePoint<2>(3.5, 0.2);
        points[7] = new ChastePoint<2>(2.6, 1.4);
        points[8] = new ChastePoint<2>(2.4, 1.5);
        points[9] = new ChastePoint<2>(3.3, 3.6);

        std::vector<Node<2>* > nodes;
        for (unsigned i=0; i<points.size(); i++)
        {
            nodes.push_back(new Node<2>(i, *(points[i]), false));
        }

        double cut_off_length = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 4.0-0.01;
        domain_size(2) = 0.0;
        domain_size(3) = 4.0-0.01; // so 4*4 boxes

        BoxCollection<2> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();


        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        std::vector< std::pair<Node<2>*, Node<2>* > > pairs_returned_vector;
        std::map<unsigned, std::set<unsigned> > neighbours_returned;

        box_collection.CalculateNodePairs(nodes,pairs_returned_vector, neighbours_returned);

        std::set< std::pair<Node<2>*, Node<2>* > > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            pairs_returned.insert(pairs_returned_vector[i]);
        }

        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
        neighbours_should_be[0].insert(1);
        neighbours_should_be[1].insert(0);
        neighbours_should_be[2].insert(3);
        neighbours_should_be[2].insert(4);
        neighbours_should_be[2].insert(5);
        neighbours_should_be[2].insert(7);
        neighbours_should_be[2].insert(8);
        neighbours_should_be[3].insert(2);
        neighbours_should_be[3].insert(4);
        neighbours_should_be[3].insert(5);
        neighbours_should_be[3].insert(7);
        neighbours_should_be[3].insert(8);
        neighbours_should_be[4].insert(2);
        neighbours_should_be[4].insert(3);
        neighbours_should_be[4].insert(5);
        neighbours_should_be[4].insert(7);
        neighbours_should_be[4].insert(8);
        neighbours_should_be[5].insert(2);
        neighbours_should_be[5].insert(3);
        neighbours_should_be[5].insert(4);
        neighbours_should_be[5].insert(6);
        neighbours_should_be[5].insert(7);
        neighbours_should_be[5].insert(8);
        neighbours_should_be[6].insert(5);
        neighbours_should_be[6].insert(7);
        neighbours_should_be[6].insert(8);
        neighbours_should_be[7].insert(2);
        neighbours_should_be[7].insert(3);
        neighbours_should_be[7].insert(4);
        neighbours_should_be[7].insert(5);
        neighbours_should_be[7].insert(6);
        neighbours_should_be[7].insert(8);
        neighbours_should_be[8].insert(2);
        neighbours_should_be[8].insert(3);
        neighbours_should_be[8].insert(4);
        neighbours_should_be[8].insert(5);
        neighbours_should_be[8].insert(6);
        neighbours_should_be[8].insert(7);
        neighbours_should_be[9] = std::set<unsigned>();

        TS_ASSERT_EQUALS(neighbours_should_be, neighbours_returned);

        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[6]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[7],nodes[8]));

        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestPairsReturned3d() throw (Exception)
    {
        // 3D cube of nodes, set up so that there is one node in each of the 3x3x3 boxes.
        std::vector<Node<3>* > nodes;
        for (unsigned k=0; k<3; k++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned i=0; i<3; i++)
                {
                    nodes.push_back(new Node<3>(i + 3*j + 9*k, false, 0.75 + 1.5*i, 0.75 + 1.5*j , 0.75 + 1.5*k));
                }
            }
        }

        double cut_off_length = 1.5;

        c_vector<double, 2*3> domain_size;  // 3x3x3 boxes
        domain_size(0) = 0.0;
        domain_size(1) = 4.4;
        domain_size(2) = 0.0;
        domain_size(3) = 4.4;
        domain_size(4) = 0.0;
        domain_size(5) = 4.4;

        BoxCollection<3> box_collection(cut_off_length, domain_size);
        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        // Make sure there is exactly one node in each box.
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            TS_ASSERT_EQUALS(box_collection.rGetBox(i).rGetNodesContained().size(), 1u);
        }

        // Calculate which pairs of node should be pairs
        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            for (unsigned j=0; j<nodes.size(); j++)
            {
                if((i < j) && norm_2(nodes[i]->rGetLocation() - nodes[j]->rGetLocation()) < 2.6)    // sqrt ( 1.5^2 + 1.5^2 + 1.5^2) rounded up.
                {
                    neighbours_should_be[i].insert(j);
                    neighbours_should_be[j].insert(i);
                }
            }
        }

        std::vector< std::pair<Node<3>*, Node<3>* > > pairs_returned_vector;
        std::map<unsigned, std::set<unsigned> > neighbours_returned;

        box_collection.CalculateNodePairs(nodes,pairs_returned_vector, neighbours_returned);

        std::set< std::pair<Node<3>*, Node<3>* > > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            pairs_returned.insert(pairs_returned_vector[i]);
        }

        // Check that the correct pairs of node 13 (central node) are in the pairs
        std::vector<unsigned> pairs_of_13;
        pairs_of_13.push_back(5);
        pairs_of_13.push_back(6);
        pairs_of_13.push_back(7);
        pairs_of_13.push_back(8);
        pairs_of_13.push_back(14);
        pairs_of_13.push_back(15);
        pairs_of_13.push_back(16);
        pairs_of_13.push_back(17);
        pairs_of_13.push_back(22);
        pairs_of_13.push_back(23);
        pairs_of_13.push_back(24);
        pairs_of_13.push_back(25);
        pairs_of_13.push_back(26);

        for (unsigned i=0; i<pairs_of_13.size(); i++)
        {
            std::pair<Node<3>*, Node<3>* > pair(nodes[13], nodes[pairs_of_13[i]]);
            TS_ASSERT(pairs_returned.find(pair) != pairs_returned.end());
        }

        // And check that others are not pairs
        std::vector<unsigned> not_pairs_of_13;
        not_pairs_of_13.push_back(0);
        not_pairs_of_13.push_back(1);
        not_pairs_of_13.push_back(2);
        not_pairs_of_13.push_back(3);
        not_pairs_of_13.push_back(4);
        not_pairs_of_13.push_back(9);
        not_pairs_of_13.push_back(10);
        not_pairs_of_13.push_back(11);
        not_pairs_of_13.push_back(13);
        not_pairs_of_13.push_back(18);
        not_pairs_of_13.push_back(19);
        not_pairs_of_13.push_back(20);
        not_pairs_of_13.push_back(21);

        for (unsigned i=0; i<not_pairs_of_13.size(); i++)
        {
            std::pair<Node<3>*, Node<3>* > pair(nodes[13], nodes[not_pairs_of_13[i]]);
            TS_ASSERT(pairs_returned.find(pair) == pairs_returned.end());
        }

        // Check the neighbour lists
        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(neighbours_should_be[i], neighbours_returned[i]);
        }

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestPairsReturned2dPeriodic() throw (Exception)
    {
        std::vector< ChastePoint<2>* > points(10);
        points[0] = new ChastePoint<2>(0.2, 3.7);
        points[1] = new ChastePoint<2>(0.5, 3.2);
        points[2] = new ChastePoint<2>(1.1, 1.99);
        points[3] = new ChastePoint<2>(1.3, 0.8);
        points[4] = new ChastePoint<2>(1.3, 0.3);
        points[5] = new ChastePoint<2>(2.2, 0.6);
        points[6] = new ChastePoint<2>(3.5, 0.2);
        points[7] = new ChastePoint<2>(2.6, 1.4);
        points[8] = new ChastePoint<2>(2.4, 1.5);
        points[9] = new ChastePoint<2>(3.3, 3.6);

        std::vector<Node<2>* > nodes;
        for (unsigned i=0; i<points.size(); i++)
        {
            nodes.push_back(new Node<2>(i, *(points[i]), false));
        }

        double cut_off_length = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 4.0-0.01;
        domain_size(2) = 0.0;
        domain_size(3) = 4.0-0.01;// so 4*4 boxes

        BoxCollection<2> box_collection(cut_off_length, domain_size, true); // Periodic in X

        box_collection.SetupLocalBoxesHalfOnly();


        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        std::vector< std::pair<Node<2>*, Node<2>* > > pairs_returned_vector;
        std::map<unsigned, std::set<unsigned> > neighbours_returned;

        box_collection.CalculateNodePairs(nodes, pairs_returned_vector, neighbours_returned);

        std::set< std::pair<Node<2>*, Node<2>* > > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            pairs_returned.insert(pairs_returned_vector[i]);
        }

        std::map<unsigned, std::set<unsigned> > neighbours_should_be;

        neighbours_should_be[0].insert(1);
        neighbours_should_be[0].insert(9);

        neighbours_should_be[1].insert(0);
        neighbours_should_be[1].insert(9);

        neighbours_should_be[2].insert(3);
        neighbours_should_be[2].insert(4);
        neighbours_should_be[2].insert(5);
        neighbours_should_be[2].insert(7);
        neighbours_should_be[2].insert(8);

        neighbours_should_be[3].insert(2);
        neighbours_should_be[3].insert(4);
        neighbours_should_be[3].insert(5);
        neighbours_should_be[3].insert(7);
        neighbours_should_be[3].insert(8);

        neighbours_should_be[4].insert(2);
        neighbours_should_be[4].insert(3);
        neighbours_should_be[4].insert(5);
        neighbours_should_be[4].insert(7);
        neighbours_should_be[4].insert(8);

        neighbours_should_be[5].insert(2);
        neighbours_should_be[5].insert(3);
        neighbours_should_be[5].insert(4);
        neighbours_should_be[5].insert(6);
        neighbours_should_be[5].insert(7);
        neighbours_should_be[5].insert(8);

        neighbours_should_be[6].insert(5);
        neighbours_should_be[6].insert(7);
        neighbours_should_be[6].insert(8);

        neighbours_should_be[7].insert(2);
        neighbours_should_be[7].insert(3);
        neighbours_should_be[7].insert(4);
        neighbours_should_be[7].insert(5);
        neighbours_should_be[7].insert(6);
        neighbours_should_be[7].insert(8);

        neighbours_should_be[8].insert(2);
        neighbours_should_be[8].insert(3);
        neighbours_should_be[8].insert(4);
        neighbours_should_be[8].insert(5);
        neighbours_should_be[8].insert(6);
        neighbours_should_be[8].insert(7);

        neighbours_should_be[9].insert(0);
        neighbours_should_be[9].insert(1);

        TS_ASSERT_EQUALS(neighbours_should_be, neighbours_returned);

        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));

        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[8]));

        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[2]));

        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[2]));

        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[6]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[2]));

        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[8]));

        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[7],nodes[8]));

        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[9],nodes[0]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[9],nodes[1]));

        TS_ASSERT_EQUALS(pairs_should_be.size(), pairs_returned.size());
        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }


    void TestBoxGeneration3d() throw (Exception)
    {
        // Create a mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(4,5,6);

        double cut_off_length = 2.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 4.15;
        domain_size(2) = -0.1;
        domain_size(3) = 5.15;
        domain_size(4) = -0.1;
        domain_size(5) = 6.15;

        BoxCollection<3> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        // Check the GetDomainSize method is working.
        const c_vector<double, 2*3> returned_domain_size = box_collection.rGetDomainSize();
        for (unsigned i=0; i < 2*3; i++)
        {
            TS_ASSERT_DELTA(returned_domain_size[i], domain_size[i], 1e-10);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 36u);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            std::set< Node<3>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
            c_vector<double, 2*3> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

            for (std::set< Node<3>* >::iterator it_nodes_in_box = nodes_in_box.begin();
                 it_nodes_in_box != nodes_in_box.end();
                 it_nodes_in_box++)
            {
                Node<3>* current_node = *it_nodes_in_box;
                double x_position = current_node->rGetLocation()[0];
                double y_position = current_node->rGetLocation()[1];
                double z_position = current_node->rGetLocation()[2];

                double epsilon = 1e-12;

                TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
                TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
                TS_ASSERT_LESS_THAN(box_min_max_values(2)-epsilon, y_position);
                TS_ASSERT_LESS_THAN(y_position, box_min_max_values(3)+epsilon);
                TS_ASSERT_LESS_THAN(box_min_max_values(4)-epsilon, z_position);
                TS_ASSERT_LESS_THAN(z_position, box_min_max_values(5)+epsilon);
            }
        }

        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(3);
        correct_answer_0.insert(4);
        correct_answer_0.insert(9);
        correct_answer_0.insert(10);
        correct_answer_0.insert(12);
        correct_answer_0.insert(13);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_13 = box_collection.GetLocalBoxes(13);
        std::set<unsigned> correct_answer_13;
        correct_answer_13.insert(5);
        correct_answer_13.insert(6);
        correct_answer_13.insert(7);
        correct_answer_13.insert(8);
        correct_answer_13.insert(13);
        correct_answer_13.insert(14);
        correct_answer_13.insert(15);
        correct_answer_13.insert(16);
        correct_answer_13.insert(17);
        correct_answer_13.insert(22);
        correct_answer_13.insert(23);
        correct_answer_13.insert(24);
        correct_answer_13.insert(25);
        correct_answer_13.insert(26);
        TS_ASSERT_EQUALS(local_boxes_to_box_13, correct_answer_13);

        std::set<unsigned> local_boxes_to_box_34 = box_collection.GetLocalBoxes(34);
        std::set<unsigned> correct_answer_34;
        correct_answer_34.insert(26);
        correct_answer_34.insert(34);
        correct_answer_34.insert(35);
        TS_ASSERT_EQUALS(local_boxes_to_box_34, correct_answer_34);

        std::set<unsigned> local_boxes_to_box_35 = box_collection.GetLocalBoxes(35);
        std::set<unsigned> correct_answer_35;
        correct_answer_35.insert(35);
        TS_ASSERT_EQUALS(local_boxes_to_box_35, correct_answer_35);
    }
};

#endif /*TESTBOXCOLLECTION_HPP_*/
