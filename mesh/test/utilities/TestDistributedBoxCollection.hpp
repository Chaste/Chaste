/*

Copyright (c) 2005-2023, University of Oxford.
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
#ifndef TESTDISTRIBUTEDBOXCOLLECTION_HPP_
#define TESTDISTRIBUTEDBOXCOLLECTION_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "TetrahedralMesh.hpp"
#include "DistributedBoxCollection.hpp"
#include "TrianglesMeshReader.hpp"
#include "ArchiveOpener.hpp"
#include "Warnings.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDistributedBoxCollection : public CxxTest::TestSuite
{
private:

    template<unsigned DIM>
    void DoUpdateHaloBoxes(unsigned numProcs)
    {
        if (3 < numProcs)
        {
            TS_TRACE("This test is only designed for 3 or fewer processes.");
            return;
        }

        // Construct a 3(x3(x3)) DistributedBoxCollection
        c_vector<double, 2*DIM> domain_size;
        for (unsigned i=0; i<DIM; i++)
        {
            domain_size[2*i] = 0.0;
            domain_size[2*i+1] = 3.0;
        }

        double box_width = 1.0;

        DistributedBoxCollection<DIM> box_collection(box_width, domain_size);


        // Put a node in each local box.
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            if (box_collection.IsBoxOwned(i))
            {
                // Create a new node and add it to the box.
                c_vector<double, DIM> node_location = zero_vector<double>(DIM);
                box_collection.rGetBox(i).AddNode(new Node<DIM>(i, node_location));
            }
            if (box_collection.IsHaloBox(i))
            {
                TS_ASSERT_THROWS_NOTHING(box_collection.rGetHaloBox(i));
            }
        }

        box_collection.UpdateHaloBoxes();

        if (!PetscTools::AmTopMost())
        {
            TS_ASSERT_EQUALS(box_collection.rGetHaloNodesRight().size(), pow(3.0, (double)DIM-1));
        }
        if (!PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(box_collection.rGetHaloNodesLeft().size(), pow(3.0, (double)DIM-1));
        }

        // Tidy up.
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            if (box_collection.IsBoxOwned(i))
            {
                std::set<Node<DIM>* > nodes = box_collection.rGetBox(i).rGetNodesContained();
                for (typename std::set<Node<DIM>* >::iterator node_iter = nodes.begin();
                     node_iter != nodes.end();
                     ++node_iter)
                {
                    delete (*node_iter);
                }
            }
        }
    }

    template<unsigned DIM>
    void DoSetupHaloBoxes(unsigned numProcs)
    {
        if (3 < numProcs)
        {
            TS_TRACE("This test is only designed for 3 or fewer processes.");
            return;
        }

        // Construct a 3(x3(x3)) DistributedBoxCollection
        c_vector<double, 2*DIM> domain_size;
        for (unsigned i=0; i<DIM; i++)
        {
            domain_size[2*i] = 0.0;
            domain_size[2*i+1] = 3.0;
        }

        double box_width = 1.0;

        DistributedBoxCollection<DIM> box_collection(box_width, domain_size);

        // Work out how many rows of boxes I have.
        std::vector<unsigned> stacks_vector;
        for (unsigned i=0;i<3;i++)
        {
            stacks_vector.push_back(i);
        }

        DistributedVectorFactory factory(3);

        Vec petsc_vec = PetscTools::CreateVec(stacks_vector.size());

        // Add data from mNumBoxesEachDirection(DIM-1)
        double* p_ret;
        VecGetArray(petsc_vec, &p_ret);
        int lo, hi;
        VecGetOwnershipRange(petsc_vec, &lo, &hi);

        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_ret[local_index] = stacks_vector[global_index];
        }
        VecRestoreArray(petsc_vec, &p_ret);

        DistributedVector* p_dist_vector = new DistributedVector(petsc_vec, &factory);

        unsigned min_row_index = (unsigned)(*p_dist_vector)[lo];
        unsigned max_row_index = (unsigned)(*p_dist_vector)[hi-1];

        // Work out how many halos there should be. # boundary processes * 3^(DIM-1)
        unsigned num_boundary_processes = 2-(unsigned)(PetscTools::AmMaster()) - (unsigned)(PetscTools::AmTopMost());
        unsigned correct_num_halos = num_boundary_processes*(unsigned)pow(3.0,(double)DIM-1);

        // Work out what the halo boxes should be.
        std::vector<unsigned> halos_should_be_right;
        std::vector<unsigned> halos_should_be_left;

        if (!PetscTools::AmTopMost())
        {
            // Only halos to the right. Code below calculates a layer of boxes.
            std::vector<unsigned> dimension_counters(DIM-1, 0);
            for (unsigned i=0; i<(unsigned)pow(3.0,(double)DIM-1); /*Increment later*/)
            {
                c_vector<unsigned, DIM> box_coords;
                for (int d=0; d< (int)DIM-1; d++)
                {
                    box_coords[d] = dimension_counters[d];
                }
                box_coords[DIM-1] = max_row_index;

                // Increment counters.
                i++;
                for (int var = 0; var < (int)DIM-1; var++)
                {
                  if (i%(unsigned)pow(3.0, (double)var) == 0)
                  {
                      dimension_counters[var] = ((dimension_counters[var]+1)%3);
                  }
                }
                unsigned global_index = box_collection.CalculateGlobalIndex(box_coords);
                halos_should_be_right.push_back(global_index);
            }
        }
        if (!PetscTools::AmMaster())
        {
            // Only halos to the right. Code below calculates a layer of boxes.
            std::vector<unsigned> dimension_counters(DIM-1, 0);
            for (unsigned i=0; i<(unsigned)pow(3.0,(double)DIM-1); /*Increment later*/)
            {
                c_vector<unsigned, DIM> box_coords;
                for (int d=0; d< (int)DIM-1; d++)
                {
                    box_coords[d] = dimension_counters[d];
                }
                box_coords[DIM-1] = min_row_index;

                // Increment counters.
                i++;
                for (int var = 0; var < (int)DIM-1; var++)
                {
                  if (i%(unsigned)pow(3.0, (double)var) == 0)
                  {
                      dimension_counters[var] = ((dimension_counters[var]+1)%3);
                  }
                }
                unsigned global_index = box_collection.CalculateGlobalIndex(box_coords);
                halos_should_be_left.push_back(global_index);
            }
        }

        TS_ASSERT_EQUALS(halos_should_be_right, box_collection.mHalosRight);
        TS_ASSERT_EQUALS(halos_should_be_left, box_collection.mHalosLeft);
        TS_ASSERT_EQUALS(box_collection.mHaloBoxes.size(),correct_num_halos);

        // Tidy up
        delete p_dist_vector;
        PetscTools::Destroy(petsc_vec);
    }

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


    void TestBoxGeneration1d()
    {
        // Create a mesh
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(20);

        double cut_off_length = 5.0;
        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 25.0;

        DistributedBoxCollection<1> box_collection(cut_off_length, domain_size);

        box_collection.SetupAllLocalBoxes();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            if (box_collection.IsBoxOwned(box_index))
            {
                box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
            }
        }

        TS_ASSERT_EQUALS(0u, box_collection.CalculateGridIndices(0)[0]);

        // The expected number of boxes is 5.
        // However, if there are more processes than required then the top-most processes will be
        // granted a box each anyway and the domain size will be swollen.
        unsigned expected_number_of_boxes = std::max(5u, PetscTools::GetNumProcs());
        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), expected_number_of_boxes);

        // Make sure the default number of rows are set for each process.
        std::vector<unsigned> rows_vector;
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            rows_vector.push_back(i);
        }
        Vec petsc_vec = PetscTools::CreateVec(rows_vector.size());
        int lo, hi;
        VecGetOwnershipRange(petsc_vec, &lo, &hi);
        PetscTools::Destroy(petsc_vec);

        unsigned local_rows = (unsigned)(hi-lo);

        TS_ASSERT_EQUALS(box_collection.GetNumLocalBoxes(), local_rows);

        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);
            std::set<unsigned> correct_answer_0;
            correct_answer_0.insert(0);
            correct_answer_0.insert(1); // Halo above
            TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
        }

        if (box_collection.IsBoxOwned(1))
        {
            std::set<unsigned>& local_boxes_to_box_1 = box_collection.rGetLocalBoxes(1);
            std::set<unsigned> correct_answer_1;
            correct_answer_1.insert(0); // Halo below
            correct_answer_1.insert(1);
            correct_answer_1.insert(2);// Halo above
            TS_ASSERT_EQUALS(local_boxes_to_box_1, correct_answer_1);
        }

        if (box_collection.IsBoxOwned(4))
        {
            std::set<unsigned>& local_boxes_to_box_4 = box_collection.rGetLocalBoxes(4);
            std::set<unsigned> correct_answer_4;
            correct_answer_4.insert(3); // Halo below
            correct_answer_4.insert(4);
            if (PetscTools::GetNumProcs() > 5u)
            {
                // There's a process (spinning) which requires an extra halo box
                correct_answer_4.insert(5);
            }
            TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);
        }

        c_vector<double,1> miles_away;
        miles_away(0) = 47323854;
        TS_ASSERT_THROWS_CONTAINS(box_collection.CalculateContainingBox(miles_away), "The point provided is outside all of the boxes");

        // Test whether we can correctly identify interior boxes
        if (PetscTools::IsSequential())
        {
            // In serial everything is an interior box.
            for (unsigned i=0; i<5u; i++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(i), true);
            }
        }
        else
        {
            TS_ASSERT_EQUALS(box_collection.IsInteriorBox(lo), false);
            for (int i=lo+1; i<hi-2; i++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(i), true);
            }
            TS_ASSERT_EQUALS(box_collection.IsInteriorBox(hi-1), false);
        }
    }


    // very simple test
    void TestAddElement()
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.5, 1.0);

        double width = 0.4;
        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.2;

        DistributedBoxCollection<1> box_collection(width, domain_size);

        if (box_collection.IsBoxOwned(0))
        {
            box_collection.rGetBox(0).AddElement(mesh.GetElement(0));
            TS_ASSERT_EQUALS(box_collection.rGetBox(0).rGetElementsContained().size(), 1u);
            TS_ASSERT_EQUALS(*(box_collection.rGetBox(0).rGetElementsContained().begin()), mesh.GetElement(0));
        }
        if (box_collection.IsBoxOwned(1))
        {
            TS_ASSERT_EQUALS(box_collection.rGetBox(1).rGetElementsContained().size(), 0u);
        }
        if (box_collection.IsBoxOwned(2))
        {
            TS_ASSERT_EQUALS(box_collection.rGetBox(2).rGetElementsContained().size(), 0u);
        }
    }

    void TestSetupAllLocalBoxes2d()
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 4.0;
        domain_size(2) = 0.0;
        domain_size(3) = 3.0;

        Warnings::Instance()->QuietDestroy();
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        DistributedBoxCollection<2> box_collection(width, domain_size);

        // We also want to consider periodic domains
        DistributedBoxCollection<2> box_collection_pdc_X(width, domain_size,true);
        DistributedBoxCollection<2> box_collection_pdc_Y(width, domain_size,false,true);
        DistributedBoxCollection<2> box_collection_pdc_XY(width, domain_size,true,true);

        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),
                             "There are more processes than convenient for the domain/mesh/box size.  The domain size has been swollen.");
            Warnings::Instance()->QuietDestroy();
        }
        // Number of slices is 3, unless there are more than 3 processes.
        // Hence the expected number of boxes is 12, but will grow with the number of processes.
        unsigned expected_number_boxes = 4u * std::max(3u, PetscTools::GetNumProcs());
        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), expected_number_boxes); // 4 * 3 boxes altogether, normally

        box_collection.SetupAllLocalBoxes();
        box_collection_pdc_X.SetupAllLocalBoxes();
        box_collection_pdc_Y.SetupAllLocalBoxes();
        box_collection_pdc_XY.SetupAllLocalBoxes();

        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);

            std::set<unsigned> correct_answer_0;
            correct_answer_0.insert(0);
            correct_answer_0.insert(1);
            correct_answer_0.insert(4); // Halo above
            correct_answer_0.insert(5); // Halo above
            TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

            // Periodic answers
            std::set<unsigned> correct_answer_0_pdc_X = correct_answer_0;
            correct_answer_0_pdc_X.insert(3);
            correct_answer_0_pdc_X.insert(7);
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(0), correct_answer_0_pdc_X);

            std::set<unsigned> correct_answer_0_pdc_Y = correct_answer_0;
            correct_answer_0_pdc_Y.insert(8);
            correct_answer_0_pdc_Y.insert(9);
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(0), correct_answer_0_pdc_Y);

            std::set<unsigned> correct_answer_0_pdc_XY = correct_answer_0_pdc_X;
            correct_answer_0_pdc_XY.insert(correct_answer_0_pdc_Y.begin(),correct_answer_0_pdc_Y.end());
            correct_answer_0_pdc_XY.insert(11);
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(0), correct_answer_0_pdc_XY);
        }
        if (box_collection.IsBoxOwned(3))
        {
            std::set<unsigned>& local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
            std::set<unsigned> correct_answer_3;
            correct_answer_3.insert(3);
            correct_answer_3.insert(2);
            correct_answer_3.insert(6);
            correct_answer_3.insert(7);
            TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

            // Periodic answers
            std::set<unsigned> correct_answer_3_pdc_X = correct_answer_3;
            correct_answer_3_pdc_X.insert(0);
            correct_answer_3_pdc_X.insert(4);
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(3), correct_answer_3_pdc_X);

            std::set<unsigned> correct_answer_3_pdc_Y = correct_answer_3;
            correct_answer_3_pdc_Y.insert(10);
            correct_answer_3_pdc_Y.insert(11);
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(3), correct_answer_3_pdc_Y);

            std::set<unsigned> correct_answer_3_pdc_XY = correct_answer_3_pdc_X;
            correct_answer_3_pdc_XY.insert(correct_answer_3_pdc_Y.begin(),correct_answer_3_pdc_Y.end());
            correct_answer_3_pdc_XY.insert(8);
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(3), correct_answer_3_pdc_XY);
        }
        if (box_collection.IsBoxOwned(5))
        {
            std::set<unsigned>& local_boxes_to_box_5 = box_collection.rGetLocalBoxes(5);
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

            // Periodic answers should be the same
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(5), correct_answer_5);
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(5), correct_answer_5);
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(5), correct_answer_5);
        }
        if (box_collection.IsBoxOwned(10))
        {
            std::set<unsigned>& local_boxes_to_box_10 = box_collection.rGetLocalBoxes(10);
            std::set<unsigned> correct_answer_10;
            correct_answer_10.insert(5); // Halo below 9
            correct_answer_10.insert(6); // Halo below 10
            correct_answer_10.insert(7); // Halo below 11
            correct_answer_10.insert(9);
            correct_answer_10.insert(10);
            correct_answer_10.insert(11);
            if (PetscTools::GetNumProcs() > 3u)
            {
                // There's a process (spinning) which requires an extra halo slice (12, 13, 14, 15)
                correct_answer_10.insert(13); // Halo above 9
                correct_answer_10.insert(14); // Halo above 10
                correct_answer_10.insert(15); // Halo above 11
            }
            TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);

            // Periodic answers
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(10), correct_answer_10);

            correct_answer_10.insert(1);
            correct_answer_10.insert(2);
            correct_answer_10.insert(3);
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(10), correct_answer_10);
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(10), correct_answer_10);
        }
    }


    void TestSetupAllLocalBoxes2dPeriodic()
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 4.0;
        domain_size(2) = 0.0;
        domain_size(3) = 3.0;

        Warnings::Instance()->QuietDestroy();
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        DistributedBoxCollection<2> box_collection(width, domain_size);

        // We also want to consider periodic domains
        DistributedBoxCollection<2> box_collection_pdc_X(width, domain_size,true);
        DistributedBoxCollection<2> box_collection_pdc_Y(width, domain_size,false,true);
        DistributedBoxCollection<2> box_collection_pdc_XY(width, domain_size,true,true);

        // Check periodic getter method
        TS_ASSERT_EQUALS(box_collection_pdc_X.GetIsPeriodicInX(),true);
        TS_ASSERT_EQUALS(box_collection_pdc_X.GetIsPeriodicInY(),false);
        TS_ASSERT_EQUALS(box_collection_pdc_X.GetIsPeriodicInZ(),false);

        TS_ASSERT_EQUALS(box_collection_pdc_Y.GetIsPeriodicInX(),false);
        TS_ASSERT_EQUALS(box_collection_pdc_Y.GetIsPeriodicInY(),true);
        TS_ASSERT_EQUALS(box_collection_pdc_Y.GetIsPeriodicInZ(),false);

        TS_ASSERT_EQUALS(box_collection_pdc_XY.GetIsPeriodicInX(),true);
        TS_ASSERT_EQUALS(box_collection_pdc_XY.GetIsPeriodicInY(),true);
        TS_ASSERT_EQUALS(box_collection_pdc_XY.GetIsPeriodicInZ(),false);


        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),
                             "There are more processes than convenient for the domain/mesh/box size.  The domain size has been swollen.");
            Warnings::Instance()->QuietDestroy();
        }
        // Number of slices is 3, unless there are more than 3 processes.
        // Hence the expected number of boxes is 12, but will grow with the number of processes.
        unsigned expected_number_boxes = 4u * std::max(3u, PetscTools::GetNumProcs());
        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), expected_number_boxes); // 4 * 3 boxes altogether, normally

        box_collection.SetupLocalBoxesHalfOnly();
        box_collection_pdc_X.SetupLocalBoxesHalfOnly();
        box_collection_pdc_Y.SetupLocalBoxesHalfOnly();
        box_collection_pdc_XY.SetupLocalBoxesHalfOnly();

        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);

            std::set<unsigned> correct_answer_0;
            correct_answer_0.insert(0);
            correct_answer_0.insert(1);
            correct_answer_0.insert(4); // Halo above
            correct_answer_0.insert(5); // Halo above
            TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

            // Periodic answers
            std::set<unsigned> correct_answer_0_x = correct_answer_0;
            correct_answer_0_x.insert(7);
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(0), correct_answer_0_x);
            if ( PetscTools::GetNumProcs() > 1 )
            {
                correct_answer_0.insert(8);
                correct_answer_0.insert(9);
            }
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(0), correct_answer_0);
            correct_answer_0.insert(7);
            if ( PetscTools::GetNumProcs() > 1 )
            {
                correct_answer_0.insert(11);
            }
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(0), correct_answer_0);
        }
        if (box_collection.IsBoxOwned(3))
        {
            std::set<unsigned>& local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
            std::set<unsigned> correct_answer_3;
            correct_answer_3.insert(3);
            correct_answer_3.insert(6);
            correct_answer_3.insert(7);
            TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

            // Periodic answers
            std::set<unsigned> correct_answer_3_y = correct_answer_3;
            if ( PetscTools::GetNumProcs() > 1 )
            {
                correct_answer_3_y.insert(10);
                correct_answer_3_y.insert(11);
            }
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(3), correct_answer_3_y);
            correct_answer_3.insert(0);
            correct_answer_3.insert(4);
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(3), correct_answer_3);
            if ( PetscTools::GetNumProcs() > 1 )
            {
                correct_answer_3.insert(8);
                correct_answer_3.insert(10);
                correct_answer_3.insert(11);
            }
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(3), correct_answer_3);
        }
        if (box_collection.IsBoxOwned(5))
        {
            std::set<unsigned>& local_boxes_to_box_5 = box_collection.rGetLocalBoxes(5);
            std::set<unsigned> correct_answer_5;
            correct_answer_5.insert(5);
            correct_answer_5.insert(6);
            correct_answer_5.insert(8);
            correct_answer_5.insert(9);
            correct_answer_5.insert(10);
            if ( PetscTools::GetMyRank() > 0 )
            {
                correct_answer_5.insert(0);
                correct_answer_5.insert(1);
                correct_answer_5.insert(2);
            }
            TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

            // Periodic answers should be the same
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(5), correct_answer_5);
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(5), correct_answer_5);
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(5), correct_answer_5);
        }
        if (box_collection.IsBoxOwned(10))
        {
            std::set<unsigned>& local_boxes_to_box_10 = box_collection.rGetLocalBoxes(10);
            std::set<unsigned> correct_answer_10;
            correct_answer_10.insert(10);
            correct_answer_10.insert(11);
            if ( PetscTools::GetMyRank() > 0 )
            {
                correct_answer_10.insert(5);
                correct_answer_10.insert(6);
                correct_answer_10.insert(7);
            }
            if (PetscTools::GetNumProcs() > 3u)
            {
                // There's a process (spinning) which requires an extra halo slice (12, 13, 14, 15)
                correct_answer_10.insert(13); // Halo above 9
                correct_answer_10.insert(14); // Halo above 10
                correct_answer_10.insert(15); // Halo above 11
            }
            TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);

            // Periodic answers
            TS_ASSERT_EQUALS(box_collection_pdc_X.rGetLocalBoxes(10), correct_answer_10);
            correct_answer_10.insert(1);
            correct_answer_10.insert(2);
            correct_answer_10.insert(3);
            TS_ASSERT_EQUALS(box_collection_pdc_Y.rGetLocalBoxes(10), correct_answer_10);
            TS_ASSERT_EQUALS(box_collection_pdc_XY.rGetLocalBoxes(10), correct_answer_10);
        }
    }


    void TestSetupAllLocalBoxes3d()
    {
        if (PetscTools::GetNumProcs() > 2)
        {
            TS_TRACE("This test is only designed for up to 2 processes");
            return;
        }
        double width = 1.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4.0;
        domain_size(2) = 0;
        domain_size(3) = 3.0;
        domain_size(4) = 0;
        domain_size(5) = 5.0;

        DistributedBoxCollection<3> box_collection(width, domain_size);

        // Also want to look at periodic in z box collections
        DistributedBoxCollection<3> box_collection_pdc_Z(width, domain_size,false,false,true);
        DistributedBoxCollection<3> box_collection_pdc_XZ(width, domain_size,true,false,true);
        DistributedBoxCollection<3> box_collection_pdc_XYZ(width, domain_size,true,true,true);

        // Check periodic getter method
        TS_ASSERT_EQUALS(box_collection_pdc_Z.GetIsPeriodicInX(),false);
        TS_ASSERT_EQUALS(box_collection_pdc_Z.GetIsPeriodicInY(),false);
        TS_ASSERT_EQUALS(box_collection_pdc_Z.GetIsPeriodicInZ(),true);

        TS_ASSERT_EQUALS(box_collection_pdc_XZ.GetIsPeriodicInX(),true);
        TS_ASSERT_EQUALS(box_collection_pdc_XZ.GetIsPeriodicInY(),false);
        TS_ASSERT_EQUALS(box_collection_pdc_XZ.GetIsPeriodicInZ(),true);

        TS_ASSERT_EQUALS(box_collection_pdc_XYZ.GetIsPeriodicInX(),true);
        TS_ASSERT_EQUALS(box_collection_pdc_XYZ.GetIsPeriodicInY(),true);
        TS_ASSERT_EQUALS(box_collection_pdc_XYZ.GetIsPeriodicInZ(),true);



        assert(box_collection.GetNumBoxes()==60); // 4 * 3 * 5 boxes altogether

        box_collection.SetupAllLocalBoxes();
        box_collection_pdc_Z.SetupAllLocalBoxes();
        box_collection_pdc_XZ.SetupAllLocalBoxes();
        box_collection_pdc_XYZ.SetupAllLocalBoxes();

        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);

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

            // Now look at the periodic
            correct_answer_0.insert(48);
            correct_answer_0.insert(49);
            correct_answer_0.insert(52);
            correct_answer_0.insert(53);
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(0), correct_answer_0);

            correct_answer_0.insert(3);
            correct_answer_0.insert(7);
            correct_answer_0.insert(15);
            correct_answer_0.insert(19);
            correct_answer_0.insert(51);
            correct_answer_0.insert(55);
            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(0), correct_answer_0);

            correct_answer_0.insert(8);
            correct_answer_0.insert(9);
            correct_answer_0.insert(11);
            correct_answer_0.insert(20);
            correct_answer_0.insert(21);
            correct_answer_0.insert(23);
            correct_answer_0.insert(56);
            correct_answer_0.insert(57);
            correct_answer_0.insert(59);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(0), correct_answer_0);
        }
        if (box_collection.IsBoxOwned(3))
        {
            std::set<unsigned>& local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
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
        }
        if (box_collection.IsBoxOwned(5))
        {
            std::set<unsigned>& local_boxes_to_box_5 = box_collection.rGetLocalBoxes(5);
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

            // Now look at the periodic
            correct_answer_5.insert(48);
            correct_answer_5.insert(49);
            correct_answer_5.insert(50);
            correct_answer_5.insert(52);
            correct_answer_5.insert(53);
            correct_answer_5.insert(54);
            correct_answer_5.insert(56);
            correct_answer_5.insert(57);
            correct_answer_5.insert(58);
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(5), correct_answer_5);

            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(5), correct_answer_5);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(5), correct_answer_5);
        }
        if (box_collection.IsBoxOwned(19))
        {
            std::set<unsigned>& local_boxes_to_box_19 = box_collection.rGetLocalBoxes(19);
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
            correct_answer_19.insert(26);
            correct_answer_19.insert(27);
            correct_answer_19.insert(30);
            correct_answer_19.insert(31);
            correct_answer_19.insert(34);
            correct_answer_19.insert(35);
            TS_ASSERT_EQUALS(local_boxes_to_box_19, correct_answer_19);

            // Now look at the periodic
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(19), correct_answer_19);
            correct_answer_19.insert(0);
            correct_answer_19.insert(4);
            correct_answer_19.insert(8);
            correct_answer_19.insert(12);
            correct_answer_19.insert(16);
            correct_answer_19.insert(20);
            correct_answer_19.insert(24);
            correct_answer_19.insert(28);
            correct_answer_19.insert(32);
            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(19), correct_answer_19);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(19), correct_answer_19);
        }
        if (box_collection.IsBoxOwned(58))
        {
            std::set<unsigned>& local_boxes_to_box_58 = box_collection.rGetLocalBoxes(58);
            std::set<unsigned> correct_answer_58;
            correct_answer_58.insert(41);
            correct_answer_58.insert(42);
            correct_answer_58.insert(43);
            correct_answer_58.insert(45);
            correct_answer_58.insert(46);
            correct_answer_58.insert(47);
            correct_answer_58.insert(53);
            correct_answer_58.insert(54);
            correct_answer_58.insert(55);
            correct_answer_58.insert(57);
            correct_answer_58.insert(58);
            correct_answer_58.insert(59);
            TS_ASSERT_EQUALS(local_boxes_to_box_58, correct_answer_58);

            // Now look at the periodic
            correct_answer_58.insert(5);
            correct_answer_58.insert(6);
            correct_answer_58.insert(7);
            correct_answer_58.insert(9);
            correct_answer_58.insert(10);
            correct_answer_58.insert(11);
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(58), correct_answer_58);
            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(58), correct_answer_58);
            correct_answer_58.insert(1);
            correct_answer_58.insert(2);
            correct_answer_58.insert(3);
            correct_answer_58.insert(37);
            correct_answer_58.insert(38);
            correct_answer_58.insert(39);
            correct_answer_58.insert(45);
            correct_answer_58.insert(46);
            correct_answer_58.insert(47);
            correct_answer_58.insert(49);
            correct_answer_58.insert(50);
            correct_answer_58.insert(51);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(58), correct_answer_58);
        }

    }

    void TestSetupLocalBoxesHalfOnly3d()
    {
        if (PetscTools::GetNumProcs() > 3)
        {
            /*
             * Note will run for more but the numebr of processors changes the local boxes as
             * boxes not on the same processor need to be included.
             *
             * With 1 processor you only have boxes above and to the left.
             * With 2 processoes you have three rows on 0 and 2 row on 1.
             * With 3 processors you have two rows on 0 two rows on 1
             * and the top row on processor.
             */
            TS_TRACE("This test is only designed for up to 3 processes");
            return;
        }
        double width = 1.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4.0;
        domain_size(2) = 0;
        domain_size(3) = 3.0;
        domain_size(4) = 0;
        domain_size(5) = 5.0;

        DistributedBoxCollection<3> box_collection(width, domain_size);

        // Also want to look at periodic in z box collections
        DistributedBoxCollection<3> box_collection_pdc_Z(width, domain_size,false,false,true);
        DistributedBoxCollection<3> box_collection_pdc_XZ(width, domain_size,true,false,true);
        DistributedBoxCollection<3> box_collection_pdc_XYZ(width, domain_size,true,true,true);

        assert(box_collection.GetNumBoxes()==60); // 4 * 3 * 5 boxes altogether

        box_collection.SetupLocalBoxesHalfOnly();
        box_collection_pdc_Z.SetupLocalBoxesHalfOnly();
        box_collection_pdc_XZ.SetupLocalBoxesHalfOnly();
        box_collection_pdc_XYZ.SetupLocalBoxesHalfOnly();

        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);

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

            // Now look at the periodic
            if (PetscTools::GetNumProcs() > 1)
            {
                // If more than one processor then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_0.insert(48);
                correct_answer_0.insert(49);
                correct_answer_0.insert(52);
                correct_answer_0.insert(53);
            }
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(0), correct_answer_0);

            if (PetscTools::GetNumProcs() > 1)
            {
                // If more than one processor then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_0.insert(51);
                correct_answer_0.insert(55);
            }
            correct_answer_0.insert(7);
            correct_answer_0.insert(15);
            correct_answer_0.insert(19);
            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(0), correct_answer_0);

            if (PetscTools::GetNumProcs() > 1)
            {
                // If more than one processor then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_0.insert(56);
                correct_answer_0.insert(57);
                correct_answer_0.insert(59);
            }
            correct_answer_0.insert(20);
            correct_answer_0.insert(21);
            correct_answer_0.insert(23);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(0), correct_answer_0);
        }
        if (box_collection.IsBoxOwned(3))
        {
            std::set<unsigned>& local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
            std::set<unsigned> correct_answer_3;
            correct_answer_3.insert(3);
            correct_answer_3.insert(6);
            correct_answer_3.insert(7);
            correct_answer_3.insert(14);
            correct_answer_3.insert(15);
            correct_answer_3.insert(18);
            correct_answer_3.insert(19);
            TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

            // Now look at the periodic

            if (PetscTools::GetNumProcs() > 1)
            {
                // If more than one processor then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_3.insert(50);
                correct_answer_3.insert(51);
                correct_answer_3.insert(54);
                correct_answer_3.insert(55);
            }
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(3), correct_answer_3);


            if (PetscTools::GetNumProcs() > 1)
            {
                // If more than one processor then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_3.insert(48);
                correct_answer_3.insert(52);
            }
            correct_answer_3.insert(0);
            correct_answer_3.insert(4);
            correct_answer_3.insert(12);
            correct_answer_3.insert(16);
            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(3), correct_answer_3);


            if (PetscTools::GetNumProcs() > 1)
            {
                // If more than one processor then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_3.insert(58);
                correct_answer_3.insert(59);
                correct_answer_3.insert(56);
            }
            correct_answer_3.insert(20);
            correct_answer_3.insert(22);
            correct_answer_3.insert(23);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(3), correct_answer_3);
        }
        if (box_collection.IsBoxOwned(5))
        {
            std::set<unsigned>& local_boxes_to_box_5 = box_collection.rGetLocalBoxes(5);
            std::set<unsigned> correct_answer_5;
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

            // Now look at the periodic
            if (PetscTools::GetNumProcs() > 1)
            {
                // If more than one processor then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_5.insert(48);
                correct_answer_5.insert(49);
                correct_answer_5.insert(50);
                correct_answer_5.insert(52);
                correct_answer_5.insert(53);
                correct_answer_5.insert(54);
                correct_answer_5.insert(56);
                correct_answer_5.insert(57);
                correct_answer_5.insert(58);
            }
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(5), correct_answer_5);
            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(5), correct_answer_5);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(5), correct_answer_5);
        }
        if (box_collection.IsBoxOwned(19))
        {
            std::set<unsigned>& local_boxes_to_box_19 = box_collection.rGetLocalBoxes(19);
            std::set<unsigned> correct_answer_19;
            correct_answer_19.insert(19);
            correct_answer_19.insert(22);
            correct_answer_19.insert(23);
            correct_answer_19.insert(26);
            correct_answer_19.insert(27);
            correct_answer_19.insert(30);
            correct_answer_19.insert(31);
            correct_answer_19.insert(34);
            correct_answer_19.insert(35);
            TS_ASSERT_EQUALS(local_boxes_to_box_19, correct_answer_19);

            // Now look at the periodic
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(19), correct_answer_19);
            correct_answer_19.insert(16);
            correct_answer_19.insert(20);
            correct_answer_19.insert(24);
            correct_answer_19.insert(28);
            correct_answer_19.insert(32);
            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(19), correct_answer_19);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(19), correct_answer_19);
        }
        if (box_collection.IsBoxOwned(58))
        {
            std::set<unsigned>& local_boxes_to_box_58 = box_collection.rGetLocalBoxes(58);
            std::set<unsigned> correct_answer_58;

            if (PetscTools::GetNumProcs() > 2 )
            {
                // If more than two processors then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_58.insert(41);
                correct_answer_58.insert(42);
                correct_answer_58.insert(43);
                correct_answer_58.insert(45);
                correct_answer_58.insert(46);
                correct_answer_58.insert(47);
            }
            correct_answer_58.insert(58);
            correct_answer_58.insert(59);
            TS_ASSERT_EQUALS(local_boxes_to_box_58, correct_answer_58);

            // Now look at the periodic
            correct_answer_58.insert(5);
            correct_answer_58.insert(6);
            correct_answer_58.insert(7);
            correct_answer_58.insert(9);
            correct_answer_58.insert(10);
            correct_answer_58.insert(11);
            TS_ASSERT_EQUALS(box_collection_pdc_Z.rGetLocalBoxes(58), correct_answer_58);

            TS_ASSERT_EQUALS(box_collection_pdc_XZ.rGetLocalBoxes(58), correct_answer_58);

            if (PetscTools::GetNumProcs() > 2 )
            {
                // If more than two processors then you need to include the
                // lower and left boxes that are on different processors
                correct_answer_58.insert(37);
                correct_answer_58.insert(38);
                correct_answer_58.insert(39);
            }
            correct_answer_58.insert(1);
            correct_answer_58.insert(2);
            correct_answer_58.insert(3);
            correct_answer_58.insert(49);
            correct_answer_58.insert(50);
            correct_answer_58.insert(51);
            TS_ASSERT_EQUALS(box_collection_pdc_XYZ.rGetLocalBoxes(58), correct_answer_58);
        }

    }

    void TestNodesPairs2DWithPeriodicity()
    {
        // Re-implemented from TestObsoleteBoxCollection

        // Set up a box collection
        c_vector<double, 2 * 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;
        domain_size(2) = 0.0;
        domain_size(3) = 1.0;

        double delta = 1e-10;
        double box_size = 0.1 + delta; // this will force 4 boxes in each dim, one with nearly no overlap

        DistributedBoxCollection<2> box_collection(box_size, domain_size, true, true);
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
         */
        for (unsigned pos_x = 0 ; pos_x < 20 ; pos_x++)
        {
            for (unsigned pos_y = 0 ; pos_y < 20 ; pos_y++)
            {
                // Random first location
                c_vector<double, 2> node_a_location = zero_vector<double>(2);
                node_a_location(0) = (0.25 + 0.5 * double(pos_x)) * box_size;
                node_a_location(1) = (0.25 + 0.5 * double(pos_y)) * box_size;

                // Check it is owned by the processor
                if (box_collection.IsBoxOwned( box_collection.CalculateContainingBox(node_a_location) ))
                {
                    for (unsigned offset = 0 ; offset < num_offsets_to_test ; offset++)
                    {
                        // Offset a second position within the interaction distance
                        c_vector<double, 2> node_b_location;
                        node_b_location = node_a_location + (box_size - delta) * offsets_to_test[offset];

                        // Account for periodicity
                        node_b_location[0] = fmod(node_b_location[0] + 1.0, 1.0);
                        node_b_location[1] = fmod(node_b_location[1] + 1.0, 1.0);

                        // Create nodes
                        std::vector<Node<2>* > nodes;
                        nodes.push_back(new Node<2>(0, node_a_location, false));
                        nodes.push_back(new Node<2>(1, node_b_location, false));

                        // Get associated boxes
                        std::vector< Box<2>* > boxes;
                        boxes.push_back(&box_collection.rGetBox(box_collection.CalculateContainingBox(node_a_location)));
                        boxes.push_back(&box_collection.rGetBox(box_collection.CalculateContainingBox(node_b_location)));

                        // Add nodes to boxes
                        boxes[0]->AddNode(nodes[0]);
                        boxes[1]->AddNode(nodes[1]);

                        std::vector< std::pair<Node<2>*, Node<2>* > > pairs_returned_vector;

                        box_collection.CalculateNodePairs(nodes, pairs_returned_vector);
                        if ( pairs_returned_vector.size() != 1 )
                        {
                            printf("For node locations (%f,%f) and (%f,%f) on process %i, no pair found. Will be in boxes %i and %i\n",
                                        node_a_location[0],node_a_location[1],node_b_location[0],node_b_location[1], (int)PetscTools::GetMyRank(),
                                        box_collection.CalculateContainingBox(node_a_location), box_collection.CalculateContainingBox(node_b_location));
                        }
                        TS_ASSERT(pairs_returned_vector.size() == 1);

                        // Remove nodes from box
                        boxes[0]->RemoveNode(nodes[0]);
                        boxes[1]->RemoveNode(nodes[1]);

                        // clean up memory
                        delete nodes[0];
                        delete nodes[1];
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    //
    //  Parallel functionality tests.
    //
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    void TestSetupHaloBoxes1d2d3d()
    {
        unsigned num_procs = PetscTools::GetNumProcs();
        DoSetupHaloBoxes<1>(num_procs);
        DoSetupHaloBoxes<2>(num_procs);
        DoSetupHaloBoxes<3>(num_procs);
    }

    void TestUpdateHaloBoxes1d2d3d()
    {
        unsigned num_procs = PetscTools::GetNumProcs();
        DoUpdateHaloBoxes<1>(num_procs);
        DoUpdateHaloBoxes<2>(num_procs);
        DoUpdateHaloBoxes<3>(num_procs);
    }

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    //
    //  Cell-based tests
    //  The following are tests written from this DistributedBoxCollection used to be
    //  cell_based/src/tissue/NodeDistributedBoxCollection and test the cell-based
    //  functionality
    //
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
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

        DistributedBoxCollection<1> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();
        TS_ASSERT_THROWS_THIS(box_collection.SetupLocalBoxesHalfOnly(), "Local Boxes Are Already Set");

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            if (box_collection.IsBoxOwned(box_index))
            {
                box_collection.IsOwned(nodes[i]);
                TS_ASSERT_EQUALS(box_collection.GetProcessOwningNode(nodes[i]), PetscTools::GetMyRank());

                box_collection.rGetBox(box_index).AddNode(nodes[i]);
            }
        }

        std::vector< std::pair<Node<1>*, Node<1>* > > pairs_returned_vector;
        box_collection.CalculateNodePairs(nodes,pairs_returned_vector);

        std::set< std::pair<Node<1>*, Node<1>* > > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            pairs_returned.insert(pairs_returned_vector[i]);
        }

        switch (PetscTools::GetNumProcs())
        {
            case 1:
            {
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

                for (unsigned i=0; i<nodes.size(); i++)
                {
                    std::vector<unsigned> expected(neighbours_should_be[i].begin(), neighbours_should_be[i].end());
                    TS_ASSERT_EQUALS(nodes[i]->rGetNeighbours(), expected);
                }

                std::set< std::pair<Node<1>*, Node<1>* > > pairs_should_be;
                pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
                pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[2]));
                pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
                pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));

                TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

                break;
            }
            case  2:
            {
                std::map<unsigned, std::set<unsigned> > neighbours_should_be;
                std::set< std::pair<Node<1>*, Node<1>* > > pairs_should_be;
                switch (PetscTools::GetMyRank())
                {
                    case 0:
                    {
                        neighbours_should_be[0].insert(1);
                        neighbours_should_be[0].insert(2);
                        neighbours_should_be[1].insert(0);
                        neighbours_should_be[1].insert(2);
                        neighbours_should_be[2].insert(0);
                        neighbours_should_be[2].insert(1);
                        neighbours_should_be[2].insert(3);
                        neighbours_should_be[3].insert(2);

                        unsigned nodes_on_this_process[4] = {0,1,2,3};

                        for (unsigned i=0; i<4; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[2]));
                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));

                        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);
                        break;

                    }
                    case 1:
                    {
                        neighbours_should_be[4] = std::set<unsigned>();

                        unsigned nodes_on_this_process[1] = {4};

                        for (unsigned i=0; i<1; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);
                        break;
                    }
                    default:
                    {
                        NEVER_REACHED;
                        break;
                    }
                }
                break;
            }

            case  3:
            {
                std::map<unsigned, std::set<unsigned> > neighbours_should_be;
                std::set< std::pair<Node<1>*, Node<1>* > > pairs_should_be;
                switch (PetscTools::GetMyRank())
                {
                    case 0:
                    {
                        neighbours_should_be[0].insert(1);
                        neighbours_should_be[0].insert(2);
                        neighbours_should_be[1].insert(0);
                        neighbours_should_be[1].insert(2);
                        neighbours_should_be[2].insert(0);
                        neighbours_should_be[2].insert(1);
                        neighbours_should_be[2].insert(3);
                        neighbours_should_be[3].insert(2);

                        unsigned nodes_on_this_process[4] = {0,1,2,3};

                        for (unsigned i=0; i<4; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[2]));
                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
                        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));

                        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);
                        break;

                    }
                    case 1:
                    {
                        // no nodes on this process

                        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);
                        break;
                    }
                    case 2:
                    {
                        neighbours_should_be[4] = std::set<unsigned>();

                        unsigned nodes_on_this_process[1] = {4};

                        for (unsigned i=0; i<1; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);
                        break;
                    }
                    default:
                    {
                        NEVER_REACHED;
                        break;
                    }
                }
                break;
            }

            default:
            {
                TS_TRACE("This test works with fewer than 4 processes.");
                return;
            }
        }

        // Tidy up.
        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestBoxGeneration2d()
    {
        if (PetscTools::GetNumProcs() > 3)
        {
            TS_TRACE("TestBoxGeneration2d only designed for up to 3 processes.");
            return;
        }
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double cut_off_length = 0.2;

        c_vector<double, 2*2> domain_size;
        domain_size[0] = -0.1;
        domain_size(1) = 1.09;
        domain_size(2) = -0.1;
        domain_size(3) = 1.09;

        DistributedBoxCollection<2> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            if (box_collection.IsBoxOwned(box_index))
            {
                box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
            }
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 36u);

        // Have checked that all the local boxes are calculated correctly on a 5 by 6 grid - here we
        // hardcode a few checks on the 6 by 6 grid.
        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);
            std::set<unsigned> correct_answer_0;
            correct_answer_0.insert(0);
            correct_answer_0.insert(1);
            correct_answer_0.insert(6);
            correct_answer_0.insert(7);
            TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
        }
        if (box_collection.IsBoxOwned(4))
        {
            std::set<unsigned>& local_boxes_to_box_4 = box_collection.rGetLocalBoxes(4);
            std::set<unsigned> correct_answer_4;
            correct_answer_4.insert(4);
            correct_answer_4.insert(5);
            correct_answer_4.insert(9);
            correct_answer_4.insert(10);
            correct_answer_4.insert(11);
            TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);
        }
        if (box_collection.IsBoxOwned(10))
        {
            std::set<unsigned>& local_boxes_to_box_10 = box_collection.rGetLocalBoxes(10);
            std::set<unsigned> correct_answer_10;
            correct_answer_10.insert(10);
            correct_answer_10.insert(11);
            correct_answer_10.insert(15);
            correct_answer_10.insert(16);
            correct_answer_10.insert(17);
            TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);
        }
        if (box_collection.IsBoxOwned(35))
        {
            std::set<unsigned>& local_boxes_to_box_35 = box_collection.rGetLocalBoxes(35);
            std::set<unsigned> correct_answer_35;
            correct_answer_35.insert(35);
            TS_ASSERT_EQUALS(local_boxes_to_box_35, correct_answer_35);
        }

        // Test whether we can correctly identify interior boxes
        if (PetscTools::IsSequential())
        {
            // In serial everything is an interior box.
            for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(i), true);
            }
        }
        else
        {
            int lo = box_collection.mpDistributedBoxStackFactory->GetLow();
            int hi = box_collection.mpDistributedBoxStackFactory->GetHigh();
            int num_face = box_collection.mNumBoxesInAFace;

            int counter;
            for (counter = lo; counter < lo + num_face; counter++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(counter), false);
            }
            for ( /*Carry on from last loop */;counter < hi-num_face + 1; counter++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(counter), true);
            }
            for ( /*Carry on from last loop */;counter < hi; counter++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(counter), false);
            }
        }
    }

    /* This test insert to verify repeatability of BoxCollection floating point
     * calculations.  Note that failure of this test on a given architecture implies
     * the failure of node-based cell simulations
     */
    void TestLargeBoxCollection2d()
    {
        double cut_off_length = 1e-3;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;
        domain_size(2) = 0.0;
        domain_size(3) = 1.0;

        DistributedBoxCollection<2> box_collection(cut_off_length, domain_size);
        TS_ASSERT_EQUALS(box_collection.mNumBoxesEachDirection[0], 1000u);
        TS_ASSERT_EQUALS(box_collection.mNumBoxesEachDirection[1], 1000u);

        c_vector<double, 2> probe;

        probe(0)=0.0; probe(1)=0.0;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 0u);

        probe(0)=1.0-1e-3; probe(1)=0.0;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 999u);

        probe(0)=0.0; probe(1)=1.0-1e-3;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 999000u);

        probe(0)=1.0-1e-3; probe(1)=1.0-1e-3;
        TS_ASSERT_EQUALS(box_collection.CalculateContainingBox(probe), 999999u);

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
        domain_size(1) = 4.0;
        domain_size(2) = 0.0;
        domain_size(3) = 4.0; // so 4*4 boxes

        DistributedBoxCollection<2> box_collection(cut_off_length, domain_size);
        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            if (box_collection.IsBoxOwned(box_index))
            {
                TS_ASSERT(box_collection.IsOwned(nodes[i]));
                TS_ASSERT_EQUALS(box_collection.GetProcessOwningNode(nodes[i]), PetscTools::GetMyRank());

                box_collection.rGetBox(box_index).AddNode(nodes[i]);
            }
        }

        std::vector< std::pair<Node<2>*, Node<2>* > > pairs_returned_vector;

        box_collection.CalculateNodePairs(nodes,pairs_returned_vector);

        // Put these into a set as they are easier to compare as they are ordered
        std::set< std::pair<Node<2>*, Node<2>* > > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            pairs_returned.insert(pairs_returned_vector[i]);
        }

        switch (PetscTools::GetNumProcs())
        {
            case 1:
            {
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

                for (unsigned i=0; i<nodes.size(); i++)
                {
                    std::vector<unsigned> expected(neighbours_should_be[i].begin(), neighbours_should_be[i].end());
                    TS_ASSERT_EQUALS(nodes[i]->rGetNeighbours(), expected);
                }

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
                break;
            }
            case 2:
            {
                switch (PetscTools::GetMyRank())
                {
                    case 0:
                    {
                        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
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

                        unsigned nodes_on_this_process[7] = {2,3,4,5,6,7,8};

                        for (unsigned i=0; i<7; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
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
                        break;
                    }
                    case 1:
                    {
                        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
                        neighbours_should_be[0].insert(1);
                        neighbours_should_be[1].insert(0);

                        neighbours_should_be[9] = std::set<unsigned>();

                        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
                        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));


                        unsigned nodes_on_this_process[2] = {0,1};

                        for (unsigned i=0; i<2; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        break;
                    }
                    default:
                    {
                        NEVER_REACHED;
                        break;
                    }
                }
                break;
            }
            case 3:
            {
                switch (PetscTools::GetMyRank())
                {
                    case 0:
                    {
                        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
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

                        unsigned nodes_on_this_process[7] = {2,3,4,5,6,7,8};

                        for (unsigned i=0; i<7; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
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
                        break;
                    }
                    case 1:
                    {
                        // No pairs on this strip.
                        break;
                    }
                    case 2:
                    {
                        std::map<unsigned, std::set<unsigned> > neighbours_should_be;
                        neighbours_should_be[0].insert(1);
                        neighbours_should_be[1].insert(0);

                        neighbours_should_be[9] = std::set<unsigned>();
                        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
                        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));

                        unsigned nodes_on_this_process[2] = {0,1};

                        for (unsigned i=0; i<2; i++)
                        {
                            unsigned node_index = nodes_on_this_process[i];
                            std::vector<unsigned> expected(neighbours_should_be[node_index].begin(),
                                                           neighbours_should_be[node_index].end());
                            TS_ASSERT_EQUALS(nodes[node_index]->rGetNeighbours(), expected);
                        }

                        break;
                    }
                    default:
                    {
                        NEVER_REACHED;
                        break;
                    }
                }
                break;
            }
            default:
            {
                // 4 or more procs
                TS_TRACE("This test is not designed for more than 3 processes");
                break;
            }
        }

        // Check we empty boxes correctly
        box_collection.EmptyBoxes();

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            if (box_collection.IsBoxOwned(i) )
            {
                TS_ASSERT_EQUALS(box_collection.rGetBox(i).rGetNodesContained().size(), 0u);
            }
        }

        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestPairsReturned3d()
    {
        if (PetscTools::GetNumProcs() > 3)
        {
            TS_TRACE("TestPairsReturned3d only designed for up to 3 processes");
            return;
        }
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
        domain_size(1) = 4.5;
        domain_size(2) = 0.0;
        domain_size(3) = 4.5;
        domain_size(4) = 0.0;
        domain_size(5) = 4.5;

        DistributedBoxCollection<3> box_collection(cut_off_length, domain_size);
        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            if (box_collection.IsBoxOwned(box_index))
            {
                // Check the box collection knows which nodes it should own.
                TS_ASSERT(box_collection.IsOwned(nodes[i]));
                TS_ASSERT_EQUALS(box_collection.GetProcessOwningNode(nodes[i]), PetscTools::GetMyRank());

                box_collection.rGetBox(box_index).AddNode(nodes[i]);
            }
            // Add as a halo if appropriate.
            if (box_collection.IsHaloBox(box_index))
            {
                box_collection.rGetHaloBox(box_index).AddNode(nodes[i]);
            }
        }

        // Make sure there is exactly one node in each box.
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            if (box_collection.IsBoxOwned(i))
            {
                TS_ASSERT_EQUALS(box_collection.rGetBox(i).rGetNodesContained().size(), 1u);
            }
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

        std::vector< std::pair<Node<3>*, Node<3>* > > pairs_returned_vector;

        box_collection.CalculateNodePairs(nodes, pairs_returned_vector);

        std::set< std::pair<Node<3>*, Node<3>* > > pairs_returned;
        for (unsigned i=0; i<pairs_returned_vector.size(); i++)
        {
            pairs_returned.insert(pairs_returned_vector[i]);
        }

        // Check that the correct pairs of node 13 (central node) are in the pairs
        std::vector<unsigned> pairs_of_13;
        pairs_of_13.push_back(14);
        pairs_of_13.push_back(15);
        pairs_of_13.push_back(16);
        pairs_of_13.push_back(17);
        pairs_of_13.push_back(18);
        pairs_of_13.push_back(19);
        pairs_of_13.push_back(20);
        pairs_of_13.push_back(21);
        pairs_of_13.push_back(22);
        pairs_of_13.push_back(23);
        pairs_of_13.push_back(24);
        pairs_of_13.push_back(25);
        pairs_of_13.push_back(26);

        // Add some extra nodes in parallel
        if (PetscTools::IsParallel())
        {
            pairs_of_13.push_back(18);
            pairs_of_13.push_back(19);
            pairs_of_13.push_back(20);
            pairs_of_13.push_back(21);

            if (PetscTools::GetNumProcs() > 2)
            {
                pairs_of_13.push_back(0);
                pairs_of_13.push_back(1);
                pairs_of_13.push_back(2);
                pairs_of_13.push_back(3);
                pairs_of_13.push_back(4);
            }
        }
        // Node 13 should be in box 13
        if (box_collection.IsBoxOwned(13))
        {
            // And check that others are not pairs
            std::vector<unsigned> not_pairs_of_13;
            if (PetscTools::GetNumProcs() < 3)
            {
                not_pairs_of_13.push_back(0);
                not_pairs_of_13.push_back(1);
                not_pairs_of_13.push_back(2);
                not_pairs_of_13.push_back(3);
                not_pairs_of_13.push_back(4);
            }

            not_pairs_of_13.push_back(9);
            not_pairs_of_13.push_back(10);
            not_pairs_of_13.push_back(11);
            not_pairs_of_13.push_back(12);
            not_pairs_of_13.push_back(13);

            for (unsigned i=0; i<not_pairs_of_13.size(); i++)
            {
                std::pair<Node<3>*, Node<3>* > pair(nodes[13], nodes[not_pairs_of_13[i]]);
                TS_ASSERT(pairs_returned.find(pair) == pairs_returned.end());
            }
        }

        // Check the neighbour lists
        for (unsigned i=0; i<nodes.size(); i++)
        {
            if (box_collection.IsBoxOwned(i))
            {
                std::vector<unsigned> expected(neighbours_should_be[i].begin(), neighbours_should_be[i].end());
                TS_ASSERT_EQUALS(nodes[i]->rGetNeighbours(), expected);
            }
        }

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestSplitNeighbourCalculation()
    {
        std::vector<Node<2>* > nodes;
        for (unsigned j=0; j<3; j++)
        {
            for (unsigned i=0; i<3; i++)
            {
                nodes.push_back(new Node<2>(i + 3*j, false, 0.75 + 1.5*i, 0.75 + 1.5*j));
            }
        }


        double cut_off_length = 1.5;

        c_vector<double, 4> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 4.5;
        domain_size(2) = 0.0;
        domain_size(3) = 4.5;

        DistributedBoxCollection<2> box_collection(cut_off_length, domain_size);
        box_collection.SetupLocalBoxesHalfOnly();


        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            if (box_collection.IsBoxOwned(box_index))
            {
                // Check the box collection knows which nodes it should own.
                TS_ASSERT(box_collection.IsOwned(nodes[i]));
                TS_ASSERT_EQUALS(box_collection.GetProcessOwningNode(nodes[i]), PetscTools::GetMyRank());

                box_collection.rGetBox(box_index).AddNode(nodes[i]);
            }
            // Add as a halo if appropriate.
            if (box_collection.IsHaloBox(box_index))
            {
                box_collection.rGetHaloBox(box_index).AddNode(nodes[i]);
            }
        }

        std::vector< std::pair<Node<2>*, Node<2>* > > pairs_returned_vector;

        box_collection.CalculateInteriorNodePairs(nodes,pairs_returned_vector);

        /* On 2 processes only neighbours of the base layer of boxes should have been calculated,
         * however there are no interior boxes on any processor as the bottom row is counted as a boundary! see #2677
         */

        // Here nodes 0-5 are on process 0 and 6-8 are on process 1
        if (PetscTools::GetNumProcs() == 2)
        {
            if (PetscTools::AmMaster())
            {
                std::vector<unsigned> neighbours_of_0 = nodes[0]->rGetNeighbours();
                std::set<unsigned> neighbours_should_be;
                neighbours_should_be.insert(1);
                neighbours_should_be.insert(3);
                neighbours_should_be.insert(4);
                std::vector<unsigned> expected(neighbours_should_be.begin(), neighbours_should_be.end());
                // This assert wasn't included and this test did nothing need to look at why it doensn't pass
                //TS_ASSERT_EQUALS(neighbours_of_0, expected);
            }

            // Check all nodes are empty!!!
            if (PetscTools::GetMyRank() == 0)
            {
                for (unsigned i=0; i<6; i++)
                {
                   TS_ASSERT_EQUALS((nodes[i]->rGetNeighbours()).size(), 0u);
                }
            }
            if (PetscTools::GetMyRank() == 1)
            {
                for (unsigned i=6; i<9; i++)
                {
                    TS_ASSERT_EQUALS((nodes[i]->rGetNeighbours()).size(), 0u);
                }
            }
        }

        /* On 3 processes nothing should have been calculated as no boxes are interior*/
        if (PetscTools::GetNumProcs() == 3)
        {
            // Note can only check on nodes owned by the process as NodeAttributed isn't set up on all nodes.
            if (PetscTools::GetMyRank() == 0)
            {
                for (unsigned i=0; i<3; i++)
                {
                   TS_ASSERT_EQUALS((nodes[i]->rGetNeighbours()).size(), 0u);
                }
            }
            if (PetscTools::GetMyRank() == 1)
            {
                for (unsigned i=3; i<6; i++)
                {
                    TS_ASSERT_EQUALS((nodes[i]->rGetNeighbours()).size(), 0u);
                }
            }
            if (PetscTools::GetMyRank() == 2)
            {
                for (unsigned i=6; i<9; i++)
                {
                    TS_ASSERT_EQUALS((nodes[i]->rGetNeighbours()).size(), 0u);
                }
            }

        }

        box_collection.CalculateBoundaryNodePairs(nodes,pairs_returned_vector);

        if (PetscTools::GetNumProcs() == 2)
        {
            if (PetscTools::AmMaster())
            {
                std::vector<unsigned> neighbours_of_4 = nodes[4]->rGetNeighbours();
                std::set<unsigned> neighbours_should_be;
                neighbours_should_be.insert(0);
                neighbours_should_be.insert(1);
                neighbours_should_be.insert(2);
                neighbours_should_be.insert(3);
                neighbours_should_be.insert(5);
                neighbours_should_be.insert(6);
                neighbours_should_be.insert(7);
                neighbours_should_be.insert(8);

                std::vector<unsigned> expected(neighbours_should_be.begin(), neighbours_should_be.end());
                TS_ASSERT_EQUALS(neighbours_of_4, expected);
            }
        }
        /* Test the all node neighbours have been calculated on 3 processes */
        if (PetscTools::GetNumProcs() == 3)
        {
            if (PetscTools::GetMyRank() == 1)
            {
                std::vector<unsigned> neighbours_of_4 = nodes[4]->rGetNeighbours();
                std::set<unsigned> neighbours_should_be;
                neighbours_should_be.insert(0);
                neighbours_should_be.insert(1);
                neighbours_should_be.insert(2);
                neighbours_should_be.insert(3);
                neighbours_should_be.insert(5);
                neighbours_should_be.insert(6);
                neighbours_should_be.insert(7);
                neighbours_should_be.insert(8);

                std::vector<unsigned> expected(neighbours_should_be.begin(), neighbours_should_be.end());
                TS_ASSERT_EQUALS(neighbours_of_4, expected);
            }
        }
        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestPairsReturned2dPeriodic()
    {
        EXIT_IF_PARALLEL;

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
        domain_size(1) = 4.0;
        domain_size(2) = 0.0;
        domain_size(3) = 4.0;// so 4*4 boxes

        DistributedBoxCollection<2> box_collection(cut_off_length, domain_size, true); // Periodic in X

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

        for (unsigned i=0; i<nodes.size(); i++)
        {
            std::vector<unsigned> expected(neighbours_should_be[i].begin(), neighbours_should_be[i].end());
            TS_ASSERT_EQUALS(nodes[i]->rGetNeighbours(), expected);
        }

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

        // Loop over and check each pair exists
        TS_ASSERT_EQUALS(pairs_should_be.size(), pairs_returned.size());
        for ( std::set< std::pair<Node<2>*, Node<2>* > >::iterator it = pairs_should_be.begin(); it!= pairs_should_be.end(); ++it )
        {
            std::set< std::pair<Node<2>*, Node<2>* > >::iterator pair_location = pairs_returned.find( *it );
            if ( pair_location == pairs_returned.end() )
            {
                // Need to check the pair isn't added as the opposite pairing
                pair_location = pairs_returned.find( std::make_pair( (*it).second, (*it).first ));
            }
            TS_ASSERT_DIFFERS( pair_location, pairs_returned.end() );
        }

        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestBoxGeneration3d()
{
        // Create a mesh
        TetrahedralMesh<3,3> mesh;

        double cut_off_length = 2.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 5.9;
        domain_size(2) = -0.1;
        domain_size(3) = 5.9;
        domain_size(4) = -0.1;
        domain_size(5) = 7.9;

        DistributedBoxCollection<3> box_collection(cut_off_length, domain_size);

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 9u * std::max(4u, PetscTools::GetNumProcs()));
        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            if (box_collection.IsBoxOwned(box_index))
            {
                box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
            }
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 9u * std::max(4u, PetscTools::GetNumProcs()));

        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);
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
        }
        if (PetscTools::GetNumProcs() <= 3u)
        {
            // If there are more processors, then there are more slices and hence more boxes.
            // We won't test for that here - there is a 2-d equivalent test.
            if (box_collection.IsBoxOwned(13))
            {
                std::set<unsigned>& local_boxes_to_box_13 = box_collection.rGetLocalBoxes(13);
                std::set<unsigned> correct_answer_13;
                correct_answer_13.insert(13);
                correct_answer_13.insert(14);
                correct_answer_13.insert(15);
                correct_answer_13.insert(16);
                correct_answer_13.insert(17);
                correct_answer_13.insert(18);
                correct_answer_13.insert(19);
                correct_answer_13.insert(20);
                correct_answer_13.insert(21);
                correct_answer_13.insert(22);
                correct_answer_13.insert(23);
                correct_answer_13.insert(24);
                correct_answer_13.insert(25);
                correct_answer_13.insert(26);
                if (PetscTools::GetNumProcs() > 1)
                {
                    correct_answer_13.insert(18);
                    correct_answer_13.insert(19);
                    correct_answer_13.insert(20);
                    correct_answer_13.insert(21);
                }
                TS_ASSERT_EQUALS(local_boxes_to_box_13, correct_answer_13);
            }

            if (box_collection.IsBoxOwned(34))
            {
                std::set<unsigned>& local_boxes_to_box_34 = box_collection.rGetLocalBoxes(34);
                std::set<unsigned> correct_answer_34;
                correct_answer_34.insert(34);
                correct_answer_34.insert(35);
                if (PetscTools::GetNumProcs() == 3)
                {
                    correct_answer_34.insert(21);
                    correct_answer_34.insert(22);
                    correct_answer_34.insert(23);
                    correct_answer_34.insert(24);
                    correct_answer_34.insert(25);
                    correct_answer_34.insert(26);
                }
                TS_ASSERT_EQUALS(local_boxes_to_box_34, correct_answer_34);
            }

            if (box_collection.IsBoxOwned(35))
            {
                std::set<unsigned>& local_boxes_to_box_35 = box_collection.rGetLocalBoxes(35);
                std::set<unsigned> correct_answer_35;
                correct_answer_35.insert(35);
                if (PetscTools::GetNumProcs() == 3)
                {
                    correct_answer_35.insert(22);
                    correct_answer_35.insert(23);
                    correct_answer_35.insert(25);
                    correct_answer_35.insert(26);
                }
                TS_ASSERT_EQUALS(local_boxes_to_box_35, correct_answer_35);
            }
        }
        // Test whether we can correctly identify interior boxes
        if (PetscTools::IsSequential())
        {
            // In serial everything is an interior box.
            for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(i), true);
            }
        }
        else
        {
            int lo = box_collection.mpDistributedBoxStackFactory->GetLow();
            int hi = box_collection.mpDistributedBoxStackFactory->GetHigh();
            int num_face = box_collection.mNumBoxesInAFace;

            int counter;
            for (counter = lo; counter < lo + num_face; counter++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(counter), false);
            }
            for ( /*Carry on from last loop */;counter < hi-num_face + 1; counter++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(counter), true);
            }
            for ( /*Carry on from last loop */;counter < hi; counter++)
            {
                TS_ASSERT_EQUALS(box_collection.IsInteriorBox(counter), false);
            }
        }
    }

    void TestArchivingDistributedBoxCollection()
    {
         FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
         std::string archive_file = "box_collection.arch";
         unsigned num_boxes = 0;

         double cut_off_length = 1.6;
         c_vector<double, 6> domain_size;
         for (unsigned i=0; i<3; i++)
         {
             domain_size[2*i] = 0.0;
             domain_size[2*i+1] = 4.8;
         }

         {
             DistributedBoxCollection<3>* p_box_collection = new DistributedBoxCollection<3>(cut_off_length, domain_size);
             p_box_collection->SetupLocalBoxesHalfOnly();
             num_boxes = p_box_collection->GetNumBoxes();
             // 4.8/1.6 = 3, so expect 3 times 3x3slices
             TS_ASSERT_EQUALS(num_boxes, 9u * std::max(3u, PetscTools::GetNumProcs()));

             {
                 // Create an output archive
                 ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
                 boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

                 // Make a const pointer
                 DistributedBoxCollection<3>* const p_const_box_collection = p_box_collection;
                 (*p_arch) << p_const_box_collection;
             }

             // Tidy up
             delete p_box_collection;
         }

         // Note if the domain has been swollen:
         if (PetscTools::GetNumProcs() > 3u)
         {
             domain_size[2*2+1] = 1.6 * PetscTools::GetNumProcs();
         }

         {
             DistributedBoxCollection<3>* p_box_collection;

             // Restore the cell population
             ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
             boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

             (*p_arch) >> p_box_collection;
             TS_ASSERT_EQUALS(num_boxes, p_box_collection->GetNumBoxes());

             for (unsigned i=0; i<3; i++)
             {
                 TS_ASSERT_DELTA(0.0, p_box_collection->rGetDomainSize()[2*i], 1e-4);
                 TS_ASSERT_DELTA(domain_size[2*i+1], p_box_collection->rGetDomainSize()[2*i+1], 1e-4);
             }

             delete p_box_collection;
         }
    }

    void TestLoadBalanceFunction()
    {
        // This test is designed for 3 process environment. Tests that an equal spread of load results
        // in an equal spread of the domain size between processes.
        if (PetscTools::GetNumProcs() == 3)
        {
            double cut_off_length = 1.0;

            c_vector<double, 2> domain_size;
            domain_size(0) = 0.0;
            domain_size(1) = 9.0;

            DistributedBoxCollection<1> box_collection(cut_off_length, domain_size);

            std::vector<int> local_loads;

            // First two processes have two rows of boxes, each with a load of 10.
            if (PetscTools::GetMyRank() < 2)
            {
                local_loads.push_back(10);
                local_loads.push_back(10);
            }
            // The top process has five rows of boxes, each with a load of 10.
            if (PetscTools::GetMyRank() == 2)
            {
                local_loads.push_back(10);
                local_loads.push_back(10);
                local_loads.push_back(10);
                local_loads.push_back(10);
                local_loads.push_back(10);
            }

            /*
             *  We expect the following number of rows to be allocated as we iterate the load-balance function, which re-assigns rows one at a time:
             *
             *  Iter    p0  p1  p2
             *  0       2   2   5   // Initial condition
             *  1       2   3   4
             *  2       3   3   3   // In 'equilibrium'
             */

            int local_rows = box_collection.LoadBalance(local_loads);
            if (PetscTools::AmMaster())
            {
                TS_ASSERT_EQUALS(local_rows, 2);
            }
            else if (PetscTools::AmTopMost())
            {
                TS_ASSERT_EQUALS(local_rows, 4);
            }
            else
            {
                TS_ASSERT_EQUALS(local_rows, 3);
            }

            // Update the load vector.
            local_loads.clear();
            local_loads.resize(local_rows, 10);

            local_rows = box_collection.LoadBalance(local_loads);
            TS_ASSERT_EQUALS(local_rows, 3);
        }
    }

    /**
     * Because of the nature of load-bal algorithm, if a domain size becomes smaller than 2 rows we
     * can encounter problems of domain shrinking to zero size at the next load balance. This test
     * makes sure that this cannot happen.
     */
    void TestLoadBalanceMaintainsMinimumLocalRegion()
    {
        if (PetscTools::GetNumProcs() == 3)
        {
            double cut_off_length = 1.0;

            c_vector<double, 2> domain_size;
            domain_size(0) = 0.0;
            domain_size(1) = 9.0;

            DistributedBoxCollection<1> box_collection(cut_off_length, domain_size);

            // Set up loads as:     1   1    |   1     100    |    1   1
            std::vector<int> local_loads;

            local_loads.resize(2,1);

            if (PetscTools::GetMyRank() == 1)
            {
                local_loads[1] = 100;
            }

            int new_rows = box_collection.LoadBalance(local_loads);
            TS_ASSERT(new_rows > 1);
        }
    }

    void TestGetDistributionOfNodes()
    {
        double cut_off_length = 1.0;

        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 9.0;

        std::vector<Node<1>* > nodes;

        DistributedBoxCollection<1> box_collection(cut_off_length, domain_size);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            if (box_collection.IsBoxOwned(i))
            {
                // Add the same number of nodes as the box index.
                for (unsigned k=0; k<i; k++)
                {
                    nodes.push_back(new Node<1>(i, false));
                    box_collection.rGetBox(i).AddNode(nodes[k]);
                }
            }
        }

        std::vector<int> local_distribution = box_collection.CalculateNumberOfNodesInEachStrip();
        int counter = 0;
        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            if (box_collection.IsBoxOwned(i))
            {
                TS_ASSERT_EQUALS(local_distribution[counter], (int)i);
                counter++;
            }
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTDISTRIBUTEDBOXCOLLECTION_HPP_*/
