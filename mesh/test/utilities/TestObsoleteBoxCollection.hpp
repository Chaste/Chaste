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

#ifndef TESTOBSOLETEBOXCOLLECTION_HPP_
#define TESTOBSOLETEBOXCOLLECTION_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "TetrahedralMesh.hpp"
#include "ObsoleteBoxCollection.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestObsoleteBoxCollection : public CxxTest::TestSuite
{
public:

    void TestBox()
    {
        Box<2> test_box;

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

    void TestIndexing()
    {
        // 1D
        {
            c_vector<double, 2> domain_size;
            domain_size(0) = 0.0; // min x
            domain_size(1) = 1.0; // max x

            ObsoleteBoxCollection<1> box_collection(0.123, domain_size);

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 9u);

            for (unsigned i = 0 ; i < box_collection.GetNumBoxes() ; i++)
            {
                c_vector<int, 1u> grid_indices = box_collection.GetGridIndices(i);
                TS_ASSERT_EQUALS(grid_indices(0), (int)i);
                TS_ASSERT_EQUALS(box_collection.GetLinearIndex(grid_indices), i);
            }
        }

        // 2D
        {
            c_vector<double, 4> domain_size;
            domain_size(0) = 0.0; // min x
            domain_size(1) = 1.0; // max x
            domain_size(2) = 0.0; // min y
            domain_size(3) = 1.0; // max y

            ObsoleteBoxCollection<2> box_collection(0.123, domain_size);
            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 81u);
        }

        // 3D
        {
            c_vector<double, 6> domain_size;
            domain_size(0) = 0.0; // min x
            domain_size(1) = 1.0; // max x
            domain_size(2) = 0.0; // min y
            domain_size(3) = 1.0; // max y
            domain_size(4) = 0.0; // min z
            domain_size(5) = 1.0; // max z

            ObsoleteBoxCollection<3> box_collection(0.123, domain_size);
            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 729u);

            c_vector<double, 6> get_domain_size = box_collection.rGetDomainSize();
            for (unsigned i=0; i<6; i++)
            {
                TS_ASSERT_DELTA(get_domain_size(i), domain_size(i), 1e-15);
            }
        }
    }

    void TestBoxGeneration1d()
    {
        // Create a mesh
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(20);

        double cut_off_length = 5.0;
        c_vector<double, 2> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 20.15;

        ObsoleteBoxCollection<1> box_collection(cut_off_length, domain_size);
        box_collection.SetupAllLocalBoxes();

        // Coverage of IsOwned()
        TS_ASSERT_EQUALS(box_collection.IsOwned(NULL), true);
        TS_ASSERT_EQUALS(box_collection.IsBoxOwned(0), true);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 5u);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            std::set< Node<1>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
            c_vector<double, 2> box_min_max_values;
            box_min_max_values(0) = i * cut_off_length - 0.1;
            box_min_max_values(1) = (i+1) * cut_off_length - 0.1;

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

        std::set<unsigned> local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_1 = box_collection.rGetLocalBoxes(1);
        std::set<unsigned> correct_answer_1;
        correct_answer_1.insert(0);
        correct_answer_1.insert(1);
        correct_answer_1.insert(2);
        TS_ASSERT_EQUALS(local_boxes_to_box_1, correct_answer_1);

        std::set<unsigned> local_boxes_to_box_4 = box_collection.rGetLocalBoxes(4);
        std::set<unsigned> correct_answer_4;
        correct_answer_4.insert(3);
        correct_answer_4.insert(4);
        TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);

        c_vector<double,1> miles_away;
        miles_away(0) = 47323854;
        TS_ASSERT_THROWS_CONTAINS(box_collection.CalculateContainingBox(miles_away), "Location in dimension 0 is");
    }

    void TestAddElement()
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.5, 1.0);

        double width = 0.4;
        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;

        ObsoleteBoxCollection<1> box_collection(width, domain_size);
        box_collection.rGetBox(0).AddElement(mesh.GetElement(0));

        TS_ASSERT_EQUALS(box_collection.rGetBox(0).rGetElementsContained().size(), 1u);
        TS_ASSERT_EQUALS(box_collection.rGetBox(1).rGetElementsContained().size(), 0u);
        TS_ASSERT_EQUALS(box_collection.rGetBox(2).rGetElementsContained().size(), 0u);
        TS_ASSERT_EQUALS(*(box_collection.rGetBox(0).rGetElementsContained().begin()), mesh.GetElement(0));
    }

    void TestSetupAllLocalBoxes2d()
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;

        ObsoleteBoxCollection<2> box_collection(width, domain_size);
        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 12u); // 4 * 3 boxes altogether

        box_collection.SetupAllLocalBoxes();

        std::set<unsigned> local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);
        unsigned correct_0[]= {0, 1, 4, 5};
        std::set<unsigned> correct_answer_0(correct_0, correct_0 + 4);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
        unsigned correct_3[]= {3, 2, 6, 7};
        std::set<unsigned> correct_answer_3(correct_3, correct_3 + 4);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.rGetLocalBoxes(5);
        unsigned correct_5[] = {0, 1, 2, 4, 5, 6, 8, 9, 10};
        std::set<unsigned> correct_answer_5(correct_5, correct_5 + 9);
        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_10 = box_collection.rGetLocalBoxes(10);
        unsigned correct_10[] = {5, 6, 7, 9, 10, 11};
        std::set<unsigned> correct_answer_10(correct_10, correct_10 + 6);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);
    }

    void TestSetupAllLocalBoxes2dPeriodic()
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;

        ObsoleteBoxCollection<2> box_collection(width, domain_size, true); // So periodic in X
        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 12u); // 4 * 3 boxes altogether

        box_collection.SetupAllLocalBoxes();

        std::set<unsigned> local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);
        unsigned correct_0[] = {0, 1, 3, 4, 5, 7};
        std::set<unsigned> correct_answer_0(correct_0, correct_0 + 6);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
        unsigned correct_3[] = {0, 2, 3, 4, 6, 7};
        std::set<unsigned> correct_answer_3(correct_3, correct_3 + 6);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.rGetLocalBoxes(5);
        unsigned correct_5[] = {0, 1, 2, 4, 5, 6, 8, 9, 10};
        std::set<unsigned> correct_answer_5(correct_5, correct_5 + 9);
        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_10 = box_collection.rGetLocalBoxes(10);
        unsigned correct_10[] = {4, 5, 6, 7, 8, 9, 10, 11};
        std::set<unsigned> correct_answer_10(correct_10, correct_10 + 8);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);

        std::set<unsigned> local_boxes_to_box_11 = box_collection.rGetLocalBoxes(11);
        unsigned correct_11[] = {4, 6, 7, 8, 10, 11};
        std::set<unsigned> correct_answer_11(correct_11, correct_11 + 6);
        TS_ASSERT_EQUALS(local_boxes_to_box_11, correct_answer_11);
    }

    void TestConvertBetweenLinearAndGridIndices()
    {
        // 1D
        {
            c_vector<double, 2 * 1> domain_size;
            domain_size(0) = 0.0;
            domain_size(1) = 0.6;

            // This interaction distance will force 6 boxes one with nearly no overlap
            double interaction_distance = 0.1001;

            ObsoleteBoxCollection<1> box_collection(interaction_distance, domain_size);

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 6u);

            // Test GetLinearIndex(GetGridIndices()) returns the same value
            for (unsigned box_idx = 0 ; box_idx < box_collection.GetNumBoxes() ; box_idx++)
            {
                TS_ASSERT_EQUALS(box_collection.GetLinearIndex(box_collection.GetGridIndices(box_idx)), box_idx);
            }
        }

        // 2D
        {
            c_vector<double, 2 * 2> domain_size;
            domain_size(0) = 0.0;
            domain_size(1) = 0.6;
            domain_size(2) = 0.0;
            domain_size(3) = 0.6;

            // This interaction distance will force 6 boxes in each dim, one with nearly no overlap
            double interaction_distance = 0.1001;

            ObsoleteBoxCollection<2> box_collection(interaction_distance, domain_size);

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 36u);

            // Test GetLinearIndex(GetGridIndices()) returns the same value
            for (unsigned box_idx = 0 ; box_idx < box_collection.GetNumBoxes() ; box_idx++)
            {
                TS_ASSERT_EQUALS(box_collection.GetLinearIndex(box_collection.GetGridIndices(box_idx)), box_idx);
            }
        }

        // 3D
        {
            c_vector<double, 2 * 3> domain_size;
            domain_size(0) = 0.0;
            domain_size(1) = 0.4;
            domain_size(2) = 0.0;
            domain_size(3) = 0.3;
            domain_size(4) = 0.0;
            domain_size(5) = 0.2;

            // This interaction distance will force 6 boxes in each dim, one with nearly no overlap
            double interaction_distance = 0.1001;

            ObsoleteBoxCollection<3> box_collection(interaction_distance, domain_size);

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 24u);

            // Test GetLinearIndex(GetGridIndices()) returns the same value
            for (unsigned box_idx = 0 ; box_idx < box_collection.GetNumBoxes() ; box_idx++)
            {
                TS_ASSERT_EQUALS(box_collection.GetLinearIndex(box_collection.GetGridIndices(box_idx)), box_idx);
            }
        }
    }

    void TestIsBoxInDomain()
    {
        c_vector<double, 2 * 3> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 0.4;
        domain_size(2) = 0.0;
        domain_size(3) = 0.3;
        domain_size(4) = 0.0;
        domain_size(5) = 0.2;

        double interaction_distance = 0.1001;

        ObsoleteBoxCollection<3> box_collection(interaction_distance, domain_size);

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 24u);

        // Test that some boxes are indeed in the domain
        TS_ASSERT(box_collection.IsBoxInDomain(box_collection.GetGridIndices(19)));
        TS_ASSERT(box_collection.IsBoxInDomain(box_collection.GetGridIndices(11)));

        // Test that the method correctly states a 'box' is not in the domain
        c_vector<int, 3> test_point_0;
        test_point_0(0) = 0;
        test_point_0(1) = 1;
        test_point_0(2) = -1;

        TS_ASSERT_EQUALS(box_collection.IsBoxInDomain(test_point_0), false);
        TS_ASSERT_EQUALS(box_collection.IsBoxInDomain(box_collection.GetGridIndices(19) + test_point_0), true);
    }

    void TestSetupAllLocalBoxes3d()
    {
        double width = 1.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;
        domain_size(4) = 0;
        domain_size(5) = 2-0.01;

        ObsoleteBoxCollection<3> box_collection(width, domain_size);

        assert(box_collection.GetNumBoxes()==24); // 4 * 3 * 2 boxes altogether

        box_collection.SetupAllLocalBoxes();

        std::set<unsigned> local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);
        unsigned correct_0[] = {0, 1, 4, 5, 12, 13, 16, 17};
        std::set<unsigned> correct_answer_0(correct_0, correct_0 + 8);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
        unsigned correct_3[] = {3, 2, 6, 7, 14, 15, 18, 19};
        std::set<unsigned> correct_answer_3(correct_3, correct_3 + 8);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.rGetLocalBoxes(5);
        unsigned correct_5[] = {0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 20, 21, 22};
        std::set<unsigned> correct_answer_5(correct_5, correct_5 + 18);
        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_19 = box_collection.rGetLocalBoxes(19);
        unsigned correct_19[] = {2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23};
        std::set<unsigned> correct_answer_19(correct_19, correct_19 + 12);
        TS_ASSERT_EQUALS(local_boxes_to_box_19, correct_answer_19);

        std::set<unsigned> local_boxes_to_box_22 = box_collection.rGetLocalBoxes(22);
        unsigned correct_22[] = {5, 6, 7, 9, 10, 11, 17, 18, 19, 21, 22, 23};
        std::set<unsigned> correct_answer_22(correct_22, correct_22 + 12);
        TS_ASSERT_EQUALS(local_boxes_to_box_22, correct_answer_22);
    }

    void TestPairsReturned1d()
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

        ObsoleteBoxCollection<1> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        std::vector< std::pair<Node<1>*, Node<1>* > > pairs_returned_vector;

        box_collection.CalculateNodePairs(nodes,pairs_returned_vector);

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
        for (unsigned i=0; i<nodes.size(); i++){
            std::vector<unsigned> expected(neighbours_should_be[i].begin(), neighbours_should_be[i].end());
            TS_ASSERT_EQUALS(nodes[i]->rGetNeighbours(), expected);
        }

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

        // Avoid memory leaks
        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestBoxGeneration2d()
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

        ObsoleteBoxCollection<2> box_collection(cut_off_length, domain_size);

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

            c_vector<double, 2*2> box_min_max_values;
            c_vector<unsigned, 2> indices = box_collection.GetGridIndices(i);
            box_min_max_values(0) = indices(0)*cut_off_length - 0.1;
            box_min_max_values(1) = (indices(0)+1)*cut_off_length - 0.1;
            box_min_max_values(2) = indices(1)*cut_off_length - 0.1;
            box_min_max_values(3) = (indices(1)+1)*cut_off_length - 0.1;


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
    }

    /*
     * This test verifies repeatability of ObsoleteBoxCollection floating point
     * calculations. Failure of this test on a given architecture implies
     * failure of node-based cell simulations.
     */
    void TestLargeObsoleteBoxCollection2d()
    {
        double cut_off_length = 1e-3;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;
        domain_size(2) = 0.0;
        domain_size(3) = 1.0;

        ObsoleteBoxCollection<2> box_collection(cut_off_length, domain_size);
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

    void TestPairsReturned2d()
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

        ObsoleteBoxCollection<2> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        std::vector< std::pair<Node<2>*, Node<2>* > > pairs_returned_vector;

        box_collection.CalculateNodePairs(nodes,pairs_returned_vector);

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
        for (unsigned i=0; i<nodes.size(); i++){
            std::vector<unsigned> expected(neighbours_should_be[i].begin(), neighbours_should_be[i].end());
            TS_ASSERT_EQUALS(nodes[i]->rGetNeighbours(), expected);
        }

        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[6]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[7],nodes[6]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[7],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[8],nodes[6]));

        TS_ASSERT_EQUALS(pairs_should_be.size(), pairs_returned.size());
        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

        // Avoid memory leaks
        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestPairsReturned3d()
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

        ObsoleteBoxCollection<3> box_collection(cut_off_length, domain_size);
        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        // Make sure there is exactly one node in each box
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            TS_ASSERT_EQUALS(box_collection.rGetBox(i).rGetNodesContained().size(), 1u);
        }

        // Calculate which pairs of nodes should be pairs
        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            for (unsigned j=0; j<nodes.size(); j++)
            {
                if ((i < j) && norm_2(nodes[i]->rGetLocation() - nodes[j]->rGetLocation()) < 2.6)    // sqrt ( 1.5^2 + 1.5^2 + 1.5^2) rounded up.
                {
                    neighbours_should_be[i].insert(j);
                    neighbours_should_be[j].insert(i);
                }
            }
        }

        std::vector< std::pair<Node<3>*, Node<3>*> > pairs_returned_vector;

        box_collection.CalculateNodePairs(nodes, pairs_returned_vector);

        std::set< std::pair<unsigned, unsigned> > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            std::pair<Node<3>*, Node<3>* > this_pair = pairs_returned_vector[i];
            unsigned this_pair_first = (this_pair.first)->GetIndex();
            unsigned this_pair_second = (this_pair.second)->GetIndex();
            std::pair<unsigned, unsigned> this_pair_indices(this_pair_first, this_pair_second);

            pairs_returned.insert(this_pair_indices);
        }

        // Check that the correct pairs of node 13 (central node) are in the pairs
        std::vector<unsigned> pairs_of_13;
        pairs_of_13.push_back(2);
        pairs_of_13.push_back(5);
        pairs_of_13.push_back(7);
        pairs_of_13.push_back(8);
        pairs_of_13.push_back(11);
        pairs_of_13.push_back(14);
        pairs_of_13.push_back(16);
        pairs_of_13.push_back(17);
        pairs_of_13.push_back(20);
        pairs_of_13.push_back(22);
        pairs_of_13.push_back(23);
        pairs_of_13.push_back(25);
        pairs_of_13.push_back(26);

        for (unsigned i=0; i<pairs_of_13.size(); i++)
        {
            std::pair<unsigned, unsigned> pair(nodes[13]->GetIndex(), nodes[pairs_of_13[i]]->GetIndex());
            TS_ASSERT(pairs_returned.find(pair) != pairs_returned.end());
        }

        // And check that others are not pairs
        std::vector<unsigned> not_pairs_of_13;
        not_pairs_of_13.push_back(0);
        not_pairs_of_13.push_back(1);
        not_pairs_of_13.push_back(3);
        not_pairs_of_13.push_back(4);
        not_pairs_of_13.push_back(6);
        not_pairs_of_13.push_back(9);
        not_pairs_of_13.push_back(10);
        not_pairs_of_13.push_back(12);
        not_pairs_of_13.push_back(13);
        not_pairs_of_13.push_back(15);
        not_pairs_of_13.push_back(18);
        not_pairs_of_13.push_back(19);
        not_pairs_of_13.push_back(21);
        not_pairs_of_13.push_back(24);

        for (unsigned i=0; i<not_pairs_of_13.size(); i++)
        {
            std::pair<unsigned, unsigned> pair(nodes[13]->GetIndex(), nodes[not_pairs_of_13[i]]->GetIndex());
            TS_ASSERT(pairs_returned.find(pair) == pairs_returned.end());
        }

        // Check the neighbour lists
        for (unsigned i=0; i<nodes.size(); i++){
            std::vector<unsigned> expected(neighbours_should_be[i].begin(), neighbours_should_be[i].end());
            TS_ASSERT_EQUALS(nodes[i]->rGetNeighbours(), expected);
        }

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestLocalBoxesHalfOnly2D()
    {
        // Define parameters for the box collection
        c_vector<double, 2 * 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;
        domain_size(2) = 0.0;
        domain_size(3) = 1.0;

        // This interaction distance will force 10 boxes in each dim, one with nearly no overlap
        double interaction_distance = 0.1001;

        // Create box collection without periodicity
        {
            ObsoleteBoxCollection<2> box_collection(interaction_distance, domain_size, false, false);
            box_collection.SetupLocalBoxesHalfOnly();

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 100u);

            std::set<unsigned> calculated_neighbours;

            // Test box 0
            calculated_neighbours = box_collection.rGetLocalBoxes(0);
            unsigned correct_0[] = {0, 10, 11, 1};
            std::set<unsigned> correct_neighbours_0(correct_0, correct_0 + 4);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_0);

            // Test representative box
            calculated_neighbours = box_collection.rGetLocalBoxes(54);
            unsigned correct_54[] = {54, 64, 65, 55, 45};
            std::set<unsigned> correct_neighbours_54(correct_54, correct_54 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_54);

            // Test top-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(92);
            unsigned correct_92[] = {92, 93, 83};
            std::set<unsigned> correct_neighbours_92(correct_92, correct_92 + 3);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_92);

            // Test penultimate-top box
            calculated_neighbours = box_collection.rGetLocalBoxes(81);
            unsigned correct_81[] = {81, 91, 92, 82, 72};
            std::set<unsigned> correct_neighbours_81(correct_81, correct_81 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_81);

            // Test right-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(49);
            unsigned correct_49[] = {49, 59};
            std::set<unsigned> correct_neighbours_49(correct_49, correct_49 + 2);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_49);

            // Test penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(58);
            unsigned correct_58[] = {58, 68, 69, 59, 49};
            std::set<unsigned> correct_neighbours_58(correct_58, correct_58 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_58);

            // Test top-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(99);
            unsigned correct_99[] = {99};
            std::set<unsigned> correct_neighbours_99(correct_99, correct_99 + 1);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_99);

            // Test penultimate-top penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(88);
            unsigned correct_88[] = {88, 98, 99, 89, 79};
            std::set<unsigned> correct_neighbours_88(correct_88, correct_88 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_88);
        }

        // Create box collection with periodicity in x
        {
            ObsoleteBoxCollection<2> box_collection(interaction_distance, domain_size, true, false);
            box_collection.SetupLocalBoxesHalfOnly();

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 100u);

            std::set<unsigned> calculated_neighbours;

            // Test box 0
            calculated_neighbours = box_collection.rGetLocalBoxes(0);
            unsigned correct_0[] = {0, 10, 11, 1};
            std::set<unsigned> correct_neighbours_0(correct_0, correct_0 + 4);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_0);

            // Test representative box
            calculated_neighbours = box_collection.rGetLocalBoxes(54);
            unsigned correct_54[] = {54, 64, 65, 55, 45};
            std::set<unsigned> correct_neighbours_54(correct_54, correct_54 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_54);

            // Test top-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(92);
            unsigned correct_92[] = {92, 93, 83};
            std::set<unsigned> correct_neighbours_92(correct_92, correct_92 + 3);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_92);

            // Test penultimate-top box
            calculated_neighbours = box_collection.rGetLocalBoxes(81);
            unsigned correct_81[] = {81, 91, 92, 82, 72};
            std::set<unsigned> correct_neighbours_81(correct_81, correct_81 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_81);

            // Test right-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(49);
            unsigned correct_49[] = {49, 59, 50, 40, 30};
            std::set<unsigned> correct_neighbours_49(correct_49, correct_49 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_49);

            // Test penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(58);
            unsigned correct_58[] = {58, 68, 69, 59, 49, 60, 50, 40};
            std::set<unsigned> correct_neighbours_58(correct_58, correct_58 + 8);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_58);

            // Test top-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(99);
            unsigned correct_99[] = {99, 90, 80};
            std::set<unsigned> correct_neighbours_99(correct_99, correct_99 + 3);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_99);

            // Test penultimate-top penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(88);
            unsigned correct_88[] = {88, 98, 99, 89, 79, 90, 80, 70};
            std::set<unsigned> correct_neighbours_88(correct_88, correct_88 + 8);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_88);
        }

        // Create box collection with periodicity in y
        {
            ObsoleteBoxCollection<2> box_collection(interaction_distance, domain_size, false, true);
            box_collection.SetupLocalBoxesHalfOnly();

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 100u);

            std::set<unsigned> calculated_neighbours;

            // Test box 0
            calculated_neighbours = box_collection.rGetLocalBoxes(0);
            unsigned correct_0[] = {0, 10, 11, 1, 91};
            std::set<unsigned> correct_neighbours_0(correct_0, correct_0 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_0);

            // Test representative box
            calculated_neighbours = box_collection.rGetLocalBoxes(54);
            unsigned correct_54[] = {54, 64, 65, 55, 45};
            std::set<unsigned> correct_neighbours_54(correct_54, correct_54 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_54);

            // Test top-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(92);
            unsigned correct_92[] = {92, 2, 3, 93, 83};
            std::set<unsigned> correct_neighbours_92(correct_92, correct_92 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_92);

            // Test penultimate-top box
            calculated_neighbours = box_collection.rGetLocalBoxes(81);
            unsigned correct_81[] = {81, 91, 92, 82, 72, 1, 2};
            std::set<unsigned> correct_neighbours_81(correct_81, correct_81 + 7);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_81);

            // Test right-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(49);
            unsigned correct_49[] = {49, 59};
            std::set<unsigned> correct_neighbours_49(correct_49, correct_49 + 2);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_49);

            // Test penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(58);
            unsigned correct_58[] = {58, 68, 69, 59, 49};
            std::set<unsigned> correct_neighbours_58(correct_58, correct_58 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_58);

            // Test top-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(99);
            unsigned correct_99[] = {99, 9};
            std::set<unsigned> correct_neighbours_99(correct_99, correct_99 + 2);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_99);

            // Test penultimate-top penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(88);
            unsigned correct_88[] = {88, 98, 99, 89, 79, 8, 9};
            std::set<unsigned> correct_neighbours_88(correct_88, correct_88 + 7);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_88);
        }

        // Create box collection with periodicity in x and y
        {
            ObsoleteBoxCollection<2> box_collection(interaction_distance, domain_size, true, true);
            box_collection.SetupLocalBoxesHalfOnly();

            TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 100u);

            std::set<unsigned> calculated_neighbours;

            // Test box 0
            calculated_neighbours = box_collection.rGetLocalBoxes(0);
            unsigned correct_0[] = {0, 10, 11, 1, 91};
            std::set<unsigned> correct_neighbours_0(correct_0, correct_0 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_0);

            // Test representative box
            calculated_neighbours = box_collection.rGetLocalBoxes(54);
            unsigned correct_54[] = {54, 64, 65, 55, 45};
            std::set<unsigned> correct_neighbours_54(correct_54, correct_54 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_54);

            // Test top-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(92);
            unsigned correct_92[] = {92, 2, 3, 93, 83};
            std::set<unsigned> correct_neighbours_92(correct_92, correct_92 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_92);

            // Test penultimate-top box
            calculated_neighbours = box_collection.rGetLocalBoxes(81);
            unsigned correct_81[] = {81, 91, 92, 82, 72, 1, 2};
            std::set<unsigned> correct_neighbours_81(correct_81, correct_81 + 7);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_81);

            // Test right-edge box
            calculated_neighbours = box_collection.rGetLocalBoxes(49);
            unsigned correct_49[] = {49, 59, 50, 40, 30};
            std::set<unsigned> correct_neighbours_49(correct_49, correct_49 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_49);

            // Test penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(58);
            unsigned correct_58[] = {58, 68, 69, 59, 49, 60, 50, 40};
            std::set<unsigned> correct_neighbours_58(correct_58, correct_58 + 8);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_58);

            // Test top-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(99);
            unsigned correct_99[] = {99, 9, 0, 90, 80};
            std::set<unsigned> correct_neighbours_99(correct_99, correct_99 + 5);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_99);

            // Test penultimate-top penultimate-right box
            calculated_neighbours = box_collection.rGetLocalBoxes(88);
            unsigned correct_88[] = {88, 98, 99, 89, 79, 8, 9, 90, 80, 70, 0};
            std::set<unsigned> correct_neighbours_88(correct_88, correct_88 + 11);
            TS_ASSERT_EQUALS(calculated_neighbours, correct_neighbours_88);
        }
    }

    void TestNodesPairs2DWithPeriodicity()
    {
        // Set up a box collection
        c_vector<double, 2 * 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;
        domain_size(2) = 0.0;
        domain_size(3) = 1.0;

        double delta = 1e-10;
        double box_size = 0.1 + delta; // this will force 4 boxes in each dim, one with nearly no overlap

        ObsoleteBoxCollection<2> box_collection(box_size, domain_size, true, true);
        box_collection.SetupLocalBoxesHalfOnly();

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 100u);

        // For each test point in the domain, we place another point within the interaction distance
        unsigned num_offsets_to_test = 16;
        std::vector<c_vector<double, 2> > offsets_to_test(num_offsets_to_test, zero_vector<double>(2));
        for (unsigned offset = 0 ; offset < num_offsets_to_test ; offset++)
        {
            double theta = 2.0 * M_PI * double(offset) / double(num_offsets_to_test);

            offsets_to_test[offset](0) = cos(theta);
            offsets_to_test[offset](1) = sin(theta);
        }

        /*
         * Systematically test points in the domain at twice the box resolution, and for each
         * point test each of the offsets calculated above, placed at the interaction distance.
         * In each case, the nodes should be considered neighbours.
         *
         * This checks 10 * 10 * 16 = 12800 node-pairs
         */
        for (unsigned pos_x = 0 ; pos_x < 20 ; pos_x++)
        {
            for (unsigned pos_y = 0 ; pos_y < 20 ; pos_y++)
            {
                // Random first location
                c_vector<double, 2> node_a_location;
                node_a_location(0) = (0.25 + 0.5 * double(pos_x)) * box_size;
                node_a_location(1) = (0.25 + 0.5 * double(pos_y)) * box_size;

                for (unsigned offset = 0 ; offset < num_offsets_to_test ; offset++)
                {
                    // Offset a second position within the interaction distance
                    c_vector<double, 2> node_b_location = node_a_location + (box_size - delta) * offsets_to_test[offset];

                    // Account for periodicity
                    node_b_location[0] = fmod(node_b_location[0] + 1.0, 1.0);
                    node_b_location[1] = fmod(node_b_location[1] + 1.0, 1.0);

                    std::vector<Node<2>* > nodes;
                    nodes.push_back(new Node<2>(0, node_a_location, false));
                    nodes.push_back(new Node<2>(1, node_b_location, false));

                    std::vector< std::pair<Node<2>*, Node<2>* > > pairs_returned_vector;

                    box_collection.CalculateNodePairs(nodes, pairs_returned_vector);

                    TS_ASSERT(pairs_returned_vector.size() == 1);

                    // clean up memory
                    delete nodes[0];
                    delete nodes[1];
                }
            }
        }
    }

    void TestNodesPairs3DWithPeriodicity()
    {
        // Set up a box collection
        c_vector<double, 2 * 3> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 0.4;
        domain_size(2) = 0.0;
        domain_size(3) = 0.4;
        domain_size(4) = 0.0;
        domain_size(5) = 0.4;

        double delta = 1e-10;
        double box_size = 0.1 + delta; // this will force 10 boxes in each dim, one with nearly no overlap

        ObsoleteBoxCollection<3> box_collection(box_size, domain_size, true, true, true);
        box_collection.SetupLocalBoxesHalfOnly();

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 64u);

        // For each test point in the domain, we place another point within the interaction distance
        unsigned num_offsets_theta = 8;
        unsigned num_offsets_phi = 8;
        std::vector<c_vector<double, 3> > offsets_to_test(num_offsets_theta * num_offsets_phi, zero_vector<double>(3));

        for (unsigned offset_phi = 0 ; offset_phi < num_offsets_phi ; offset_phi++)
        {
            double phi = M_PI * double(offset_phi) / double(num_offsets_phi);

            for (unsigned offset_theta = 0 ; offset_theta < num_offsets_theta ; offset_theta++)
            {
                double theta = 2.0 * M_PI * double(offset_theta) / double(num_offsets_theta);

                unsigned index = offset_phi * num_offsets_theta + offset_theta;

                offsets_to_test[index](0) = cos(theta);
                offsets_to_test[index](1) = cos(phi) * sin(theta);
                offsets_to_test[index](2) = -sin(phi) * sin(theta);
            }
        }

        /*
         * Systematically test points in the domain at twice the box resolution, and for each
         * point test each of the offsets calculated above, placed at the interaction distance.
         * In each case, the nodes should be considered neighbours.
         *
         * This checks 8 * 8 * 8 * 16 * 8 = 65536 node-pairs
         */
        for (unsigned pos_x = 0 ; pos_x < 6 ; pos_x++)
        {
            for (unsigned pos_y = 0 ; pos_y < 6 ; pos_y++)
            {
                for (unsigned pos_z = 0 ; pos_z < 6 ; pos_z++)
                {
                    // First location
                    c_vector<double, 3> node_a_location;
                    node_a_location(0) = (0.25 + 0.5 * double(pos_x)) * box_size;
                    node_a_location(1) = (0.25 + 0.5 * double(pos_y)) * box_size;
                    node_a_location(2) = (0.25 + 0.5 * double(pos_z)) * box_size;

                    for (unsigned offset = 0 ; offset < offsets_to_test.size() ; offset++)
                    {
                        // Offset a second position within the interaction distance
                        c_vector<double, 3> node_b_location = node_a_location + (box_size - delta) * offsets_to_test[offset];

                        // Account for periodicity
                        node_b_location[0] = fmod(node_b_location[0] + 1.0, 1.0);
                        node_b_location[1] = fmod(node_b_location[1] + 1.0, 1.0);
                        node_b_location[2] = fmod(node_b_location[2] + 1.0, 1.0);

                        std::vector<Node<3>* > nodes;
                        nodes.push_back(new Node<3>(0, node_a_location, false));
                        nodes.push_back(new Node<3>(1, node_b_location, false));

                        std::vector< std::pair<Node<3>*, Node<3>* > > pairs_returned_vector;

//                        box_collection.CalculateNodePairs(nodes, pairs_returned_vector, neighbours_returned);

                        ///\todo Add a test here #2725

                        // Free memory
                        delete nodes[0];
                        delete nodes[1];
                    }
                }
            }
        }
    }
};

#endif /*TESTOBSOLETEBOXCOLLECTION_HPP_*/
