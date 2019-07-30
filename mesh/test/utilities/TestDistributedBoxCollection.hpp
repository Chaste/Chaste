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

        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);

            std::set<unsigned> correct_answer_0;
            correct_answer_0.insert(0);
            correct_answer_0.insert(1);
            correct_answer_0.insert(4); // Halo above
            correct_answer_0.insert(5); // Halo above
            TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
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
        }
    }


    void TestSetupAllLocalBoxes2dPeriodic()
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4.0;
        domain_size(2) = 0;
        domain_size(3) = 3.0;

        Warnings::Instance()->QuietDestroy();
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        DistributedBoxCollection<2> box_collection(width, domain_size, true); // So periodic in X

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
        if (box_collection.IsBoxOwned(0))
        {
            std::set<unsigned>& local_boxes_to_box_0 = box_collection.rGetLocalBoxes(0);

            std::set<unsigned> correct_answer_0;
            correct_answer_0.insert(0);
            correct_answer_0.insert(1);
            correct_answer_0.insert(3);
            correct_answer_0.insert(4);
            correct_answer_0.insert(5);
            correct_answer_0.insert(7);
            TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);
        }
        if (box_collection.IsBoxOwned(3))
        {
            std::set<unsigned>& local_boxes_to_box_3 = box_collection.rGetLocalBoxes(3);
            std::set<unsigned> correct_answer_3;
            correct_answer_3.insert(0);
            correct_answer_3.insert(2);
            correct_answer_3.insert(3);
            correct_answer_3.insert(4);
            correct_answer_3.insert(6);
            correct_answer_3.insert(7);
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
            TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);
        }
        if (box_collection.IsBoxOwned(10))
        {
            std::set<unsigned>& local_boxes_to_box_10 = box_collection.rGetLocalBoxes(10);
            std::set<unsigned> correct_answer_10;
            correct_answer_10.insert(5);
            correct_answer_10.insert(6);
            correct_answer_10.insert(7);
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
        }
        if (box_collection.IsBoxOwned(11))
        {
            std::set<unsigned>& local_boxes_to_box_11 = box_collection.rGetLocalBoxes(11);
            std::set<unsigned> correct_answer_11;
            correct_answer_11.insert(4);
            correct_answer_11.insert(6);
            correct_answer_11.insert(7);
            correct_answer_11.insert(8);
            correct_answer_11.insert(10);
            correct_answer_11.insert(11);
            if (PetscTools::GetNumProcs() > 3u)
            {
                // There's a process (spinning) which requires an extra halo slice (12, 13, 14, 15)
                correct_answer_11.insert(14); // Halo above 10
                correct_answer_11.insert(15); // Halo above 11
                correct_answer_11.insert(12); // Halo above 8 (Periodic in x)
            }
           TS_ASSERT_EQUALS(local_boxes_to_box_11, correct_answer_11);
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
        domain_size(5) = 2.0;

        DistributedBoxCollection<3> box_collection(width, domain_size);

        assert(box_collection.GetNumBoxes()==24); // 4 * 3 * 2 boxes altogether

        box_collection.SetupAllLocalBoxes();
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
            TS_ASSERT_EQUALS(local_boxes_to_box_19, correct_answer_19);
        }
        if (box_collection.IsBoxOwned(22))
        {
            std::set<unsigned>& local_boxes_to_box_22 = box_collection.rGetLocalBoxes(22);
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
        domain_size(0) = -0.1;
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

            if (PetscTools::GetNumProcs() < 2)
            {
                not_pairs_of_13.push_back(18);
                not_pairs_of_13.push_back(19);
                not_pairs_of_13.push_back(20);
                not_pairs_of_13.push_back(21);
            }

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

        TS_ASSERT_EQUALS(pairs_should_be.size(), pairs_returned.size());
        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

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
        mesh.ConstructCuboid(4,5,6);

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
                correct_answer_34.insert(26);
                correct_answer_34.insert(34);
                correct_answer_34.insert(35);
                if (PetscTools::GetNumProcs() == 3)
                {
                    correct_answer_34.insert(21);
                    correct_answer_34.insert(22);
                    correct_answer_34.insert(23);
                    correct_answer_34.insert(24);
                    correct_answer_34.insert(25);
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
