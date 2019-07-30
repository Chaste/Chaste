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

#ifndef TESTT2SWAPCELLKILLER_HPP_
#define TESTT2SWAPCELLKILLER_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "T2SwapCellKiller.hpp"
#include "OffLatticeSimulation.hpp"
#include "FileComparison.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Warnings.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "FakePetscSetup.hpp"

/**
 * This class contains tests for the cell killer that performs T2 swaps.
 */
class TestT2SwapCellKiller : public AbstractCellBasedTestSuite
{
public:

    /**
     * This tests the the killer method CheckAndLabelCellsForApoptosisOrDeath()
     * method for performing T2 swaps (element removal).
     */
    void TestKillerForT2Swap()
    {
        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.1, 0.05));
        nodes.push_back(new Node<2>(4, false, 0.9, 0.05));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.475));

        // Make one triangular and three trapezium elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {3, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_2[4] = {2, 0, 3, 5};
        unsigned node_indices_elem_3[4] = {0, 1, 4, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        vertex_mesh.SetT2Threshold(0.01);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // The population should have 4 cells
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u)

        // Give the Population to the cell killer
        T2SwapCellKiller<2> cell_killer(&cell_population);

        // Perform swaps
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // We should not have had any T2 swaps yet
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // We move the inner vertices a bit inwards
        c_vector<double, 2>& new_location_0 = vertex_elements[0]->GetNode(0)->rGetModifiableLocation();
        new_location_0(0) = 0.499;
        new_location_0(1) = 0.249;

        c_vector<double, 2>& new_location_1 = vertex_elements[0]->GetNode(1)->rGetModifiableLocation();
        new_location_1(0) = 0.501;
        new_location_1(1) = 0.249;

        c_vector<double, 2>& new_location_2 = vertex_elements[0]->GetNode(2)->rGetModifiableLocation();
        new_location_2(0) = 0.5;
        new_location_2(1) = 0.251;

        // T2 swaps should now be able to happen

        // Perform swaps
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        // Test that the merged node has the correct location
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 0.4999, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.2496, 1e-3);

        // Test tracking of T2 swaps
        std::vector< c_vector<double, 2> > t2_locations = cell_population.GetLocationsOfT2Swaps();
        std::vector< unsigned > t2_cell_ids = cell_population.GetCellIdsOfT2Swaps();
        TS_ASSERT_EQUALS(t2_locations.size(), 1u);
        TS_ASSERT_EQUALS(t2_cell_ids.size(), 1u);
        TS_ASSERT_EQUALS(t2_cell_ids[0], 0u);
        TS_ASSERT_DELTA(t2_locations[0][0], 0.4999, 1e-3);
        TS_ASSERT_DELTA(t2_locations[0][1], 0.2496, 1e-3);

        // Test T2 swap Location clearing
        cell_population.ClearLocationsAndCellIdsOfT2Swaps();
        t2_locations = cell_population.GetLocationsOfT2Swaps();
        t2_cell_ids = cell_population.GetCellIdsOfT2Swaps();
        TS_ASSERT_EQUALS(t2_locations.size(), 0u);
        TS_ASSERT_EQUALS(t2_cell_ids.size(), 0u);

        // Test that each element contains the correct number of nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);

        // Test that each element contains the correct nodes following the rearrangement
        // (note that the node indices are reset since element 0 was deleted from the mesh)
        unsigned node_indices_element_0[3] = {1, 2, 6};
        unsigned node_indices_element_1[3] = {2, 0, 6};
        unsigned node_indices_element_2[3] = {0, 1, 6};
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
        }

        // Test that the cell corresponding to element 0 has been
        // marked as deleted and the others haven't
        TS_ASSERT(cell_population.GetCellUsingLocationIndex(0)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(1)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(2)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(3)->IsDead());

        // We shouldn't have overseen any cells
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(),4u);
    }

    void TestKillerForT2SwapInSimulation()
    {
        /**
         * This is performs a single T2 swap in a simulation and tests that the cells and vertex elements are
         * deleted correctly.
         */
        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.1, 0.05));
        nodes.push_back(new Node<2>(4, false, 0.9, 0.05));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.475));

        // Make one triangular and three trapezium elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {3, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_2[4] = {2, 0, 3, 5};
        unsigned node_indices_elem_3[4] = {0, 1, 4, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        vertex_mesh.SetT2Threshold(0.01);

        // Prevent T1 swaps from occurring by setting the threshold distance between vertices to be very small
        vertex_mesh.SetCellRearrangementThreshold(0.00001);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // make a simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestT2SwapCellKillerInSimulation");
        simulator.SetEndTime(0.003);

        // Perform swaps
        simulator.Solve();

        // We should not have had any T2 swaps yet
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Move the inner vertices inwards
        c_vector<double, 2>& new_location_0 = vertex_elements[0]->GetNode(0)->rGetModifiableLocation();
        new_location_0(0) = 0.499;
        new_location_0(1) = 0.249;

        c_vector<double, 2>& new_location_1 = vertex_elements[0]->GetNode(1)->rGetModifiableLocation();
        new_location_1(0) = 0.501;
        new_location_1(1) = 0.249;

        c_vector<double, 2>& new_location_2 = vertex_elements[0]->GetNode(2)->rGetModifiableLocation();
        new_location_2(0) = 0.5;
        new_location_2(1) = 0.251;

        // T2 swaps should now happen
        simulator.SetEndTime(0.005);

        // Perform swaps
        simulator.Solve();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        // Test that the merged node has the correct location
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[0], 0.4999, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[1], 0.2496, 1e-3);

        // Test that each element contains the correct number of nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);

        // Test that each element contains the correct nodes following the rearrangement
        // (note that the node indices are reset since element 0 was deleted from the mesh)
        unsigned node_indices_element_0[3] = {1, 2, 3};
        unsigned node_indices_element_1[3] = {2, 0, 3};
        unsigned node_indices_element_2[3] = {0, 1, 3};
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
        }

        // Check that we have the right number of cells and that none of the cells are marked as dead
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 3u)
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(0)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(1)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(2)->IsDead());

        // We also do not have any undeleted cells
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(),3u);
    }

    void TestKillerForMultipleT2Swaps()
    {
        /**
         * Create a mesh comprising ten nodes contained in six elements, two of which are small triangles,
         * as shown below. We will test that the CheckAndLabelCellsForApoptosisOrDeath() method works correctly in the case where multiple
         * T2 swaps are required. After remeshing, the two triangular elements should be removed from the mesh.
         *
         *   ________          _________
         *  |\      /|        |\       /|
         *  | |\__/| | -----> | \_____/ |
         *  | |/  \| |        | /     \ |
         *  |/______\|        |/_______\|
         *
         */
        std::vector<Node<2>*> nodes;
        // Create the nodes
        // The boolean is true for boundary nodes
        nodes.push_back(new Node<2>(0, true,  0.0, 0.0));
        nodes.push_back(new Node<2>(1, true,  1.0, 0.0));
        nodes.push_back(new Node<2>(2, true,  1.0, 1.0));
        nodes.push_back(new Node<2>(3, true,  0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.3, 0.45));
        nodes.push_back(new Node<2>(5, false, 0.4, 0.5));
        nodes.push_back(new Node<2>(6, false, 0.3, 0.55));
        nodes.push_back(new Node<2>(7, false, 0.6, 0.5));
        nodes.push_back(new Node<2>(8, false, 0.7, 0.45));
        nodes.push_back(new Node<2>(9, false, 0.7, 0.55));

        // Construct elements from the nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3, nodes_elem_4, nodes_elem_5;
        unsigned node_indices_elem_0[6] = {0, 1, 8, 7, 5, 4};
        unsigned node_indices_elem_1[6] = {2, 3, 6, 5, 7, 9};
        unsigned node_indices_elem_2[4] = {1, 2, 9, 8};
        unsigned node_indices_elem_3[4] = {0, 4, 6, 3};
        unsigned node_indices_elem_4[3] = {4, 5, 6};
        unsigned node_indices_elem_5[3] = {7, 8, 9};

        for (unsigned i=0; i<6; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            if (i < 4)
            {
                nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
                nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
            }
            if (i < 3)
            {
                nodes_elem_4.push_back(nodes[node_indices_elem_4[i]]);
                nodes_elem_5.push_back(nodes[node_indices_elem_5[i]]);
            }
        }

        // create a mesh
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(4, nodes_elem_4));
        vertex_elements.push_back(new VertexElement<2,2>(5, nodes_elem_5));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Set the remeshing threshold parameter values so that T2 swaps do occur, but T1 swaps do not
        vertex_mesh.SetT2Threshold(0.01);
        vertex_mesh.SetCellRearrangementThreshold(0.00001);

        // Test that the numbers of nodes and elements are correct
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // the population should have 6 cells
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 6u)

        // ... and give it to the killer
        T2SwapCellKiller<2> cell_killer(&cell_population);

        // Perform swaps
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Test that after remeshing we have the correct number of nodes and elements
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test that after remeshing, each element contains the correct number of nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);

        // Test that each merged node is correctly located at the centroid of the element it has replaced
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 1.0/3.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-3);

        // The new node is either added to the end of the nodes vector or it replaces the deleted node with the
        // largest index among the deleted nodes. Example: above we had 9 nodes and added a 10th, and then deleted
        // nodes 4, 5, and 6. So now we wrote into node 6 before we deleted node 7, 8, and 9.
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 2.0/3.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.5, 1e-3);

        // Test tracking of T2 swaps
        std::vector< c_vector<double, 2> > t2_locations = cell_population.GetLocationsOfT2Swaps();
        std::vector< unsigned > t2_cell_ids = cell_population.GetCellIdsOfT2Swaps();
        TS_ASSERT_EQUALS(t2_locations.size(), 2u);
        TS_ASSERT_EQUALS(t2_cell_ids.size(), 2u);
        TS_ASSERT_EQUALS(t2_cell_ids[0], 4u);
        TS_ASSERT_EQUALS(t2_cell_ids[1], 5u);
        TS_ASSERT_DELTA(t2_locations[0][0], 1.0/3.0, 1e-3);
        TS_ASSERT_DELTA(t2_locations[0][1], 0.5, 1e-3);
        TS_ASSERT_DELTA(t2_locations[1][0], 2.0/3.0, 1e-3);
        TS_ASSERT_DELTA(t2_locations[1][1], 0.5, 1e-3);

        // Test T2 swap Location clearing
        cell_population.ClearLocationsAndCellIdsOfT2Swaps();
        t2_locations = cell_population.GetLocationsOfT2Swaps();
        t2_cell_ids = cell_population.GetCellIdsOfT2Swaps();
        TS_ASSERT_EQUALS(t2_locations.size(), 0u);
        TS_ASSERT_EQUALS(t2_cell_ids.size(), 0u);

        // Test that after remeshing, each element contains the correct nodes
        unsigned node_indices_element_0[4] = {0, 1, 6, 10};
        unsigned node_indices_element_1[4] = {2, 3, 10, 6};
        unsigned node_indices_element_2[3] = {1, 2, 6};
        unsigned node_indices_element_3[3] = {0, 10, 3};
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
        }
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), node_indices_element_0[3]);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), node_indices_element_1[3]);

        // Test that the cell corresponding to element 4 and 5 has been
        // marked as deleted and the others haven't
        TS_ASSERT(cell_population.GetCellUsingLocationIndex(4)->IsDead());
        TS_ASSERT(cell_population.GetCellUsingLocationIndex(5)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(0)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(1)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(2)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(3)->IsDead());

        // None of the cells should be deleted yet
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 6u);
    }

    void TestKillerForMultipleT2SwapInSimulation()
    {
        /**
         * We conduct the same test as before, but now within an OffLatticeSimulation. We make sure that
         * all killed cells and the corresponding vertex elements get correctly deleted.
         *
         *   ________          _________
         *  |\      /|        |\       /|
         *  | |\__/| | -----> | \_____/ |
         *  | |/  \| |        | /     \ |
         *  |/______\|        |/_______\|
         *
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true,  0.0, 0.0));
        nodes.push_back(new Node<2>(1, true,  1.0, 0.0));
        nodes.push_back(new Node<2>(2, true,  1.0, 1.0));
        nodes.push_back(new Node<2>(3, true,  0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.3, 0.45));
        nodes.push_back(new Node<2>(5, false, 0.4, 0.5));
        nodes.push_back(new Node<2>(6, false, 0.3, 0.55));
        nodes.push_back(new Node<2>(7, false, 0.6, 0.5));
        nodes.push_back(new Node<2>(8, false, 0.7, 0.45));
        nodes.push_back(new Node<2>(9, false, 0.7, 0.55));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3, nodes_elem_4, nodes_elem_5;
        unsigned node_indices_elem_0[6] = {0, 1, 8, 7, 5, 4};
        unsigned node_indices_elem_1[6] = {2, 3, 6, 5, 7, 9};
        unsigned node_indices_elem_2[4] = {1, 2, 9, 8};
        unsigned node_indices_elem_3[4] = {0, 4, 6, 3};
        unsigned node_indices_elem_4[3] = {4, 5, 6};
        unsigned node_indices_elem_5[3] = {7, 8, 9};

        for (unsigned i=0; i<6; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            if (i < 4)
            {
                nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
                nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
            }
            if (i < 3)
            {
                nodes_elem_4.push_back(nodes[node_indices_elem_4[i]]);
                nodes_elem_5.push_back(nodes[node_indices_elem_5[i]]);
            }
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(4, nodes_elem_4));
        vertex_elements.push_back(new VertexElement<2,2>(5, nodes_elem_5));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Set the remeshing threshold parameter values so that T2 swaps do occur, but T1 swaps do not
        vertex_mesh.SetT2Threshold(0.01);
        vertex_mesh.SetCellRearrangementThreshold(0.00001);

        // Test that the numbers of nodes and elements are correct before and after remeshing
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // the population should have 6 cells
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 6u)

        // make a simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestT2SwapCellKillerMultipleT2SwapsInSimulation");
        simulator.SetEndTime(0.003);

        // Perform swaps
        simulator.Solve();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test that after remeshing, each element contains the correct number of nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);

        // Test that each merged node is correctly located at the centroid of the element it has replaced
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 1.0/3.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // The node that we wrote the latest should have a smaller index than the one that we wrote first.
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 2.0/3.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-3);

        // Test that after remeshing, each element contains the correct nodes
        unsigned node_indices_element_0[4] = {0, 1, 4, 5};
        unsigned node_indices_element_1[4] = {2, 3, 5, 4};
        unsigned node_indices_element_2[3] = {1, 2, 4};
        unsigned node_indices_element_3[3] = {0, 5, 3};
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
        }
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), node_indices_element_0[3]);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), node_indices_element_1[3]);

        // Test that the cell corresponding to element 4 and 5 have been
        // deleted and none of the cells is marked as dead

        TS_ASSERT_EQUALS(cell_population.rGetCells().size(),4u);

        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(0)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(1)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(2)->IsDead());
        TS_ASSERT(!cell_population.GetCellUsingLocationIndex(3)->IsDead());
    }

    /**
     * This tests that T1 swaps rearrange to form a triangular element for a T2 swap
     */
    void TestPrepareForT2Swap()
    {
        /*
         * Create a mesh comprising eight nodes contained in four trapezium elements and a central
         * square element, as shown below. We will test that a T1 swap correctly results in a triangle
         * element, which is then correctly removed by a T2 swap.
         *  _______
         * |\  2  /|
         * | \___/ |
         * | |   | |
         * |3| 4 |1|
         * | |___| |
         * | / 0 \ |
         * |/_____\|
         */
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, -1.0, -1.0));
        nodes.push_back(new Node<2>(1, false,  1.0, -1.0));
        nodes.push_back(new Node<2>(2, false,  1.0,  1.0));
        nodes.push_back(new Node<2>(3, false, -1.0,  1.0));
        nodes.push_back(new Node<2>(4, false, -0.1, -0.1));
        nodes.push_back(new Node<2>(5, false,  0.1, -0.1));
        nodes.push_back(new Node<2>(6, false,  0.1,  0.1));
        nodes.push_back(new Node<2>(7, false, -0.1,  0.1));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3, nodes_elem_4;
        unsigned node_indices_elem_0[4] = {0, 1, 5, 4};
        unsigned node_indices_elem_1[4] = {1, 2, 6, 5};
        unsigned node_indices_elem_2[4] = {2, 3, 7, 6};
        unsigned node_indices_elem_3[4] = {0, 4, 7, 3};
        unsigned node_indices_elem_4[4] = {4, 5, 6, 7};
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
            nodes_elem_4.push_back(nodes[node_indices_elem_4[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(4, nodes_elem_4));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Set the threshold distance between nodes for T1 swaps to be longer than the inner edges of the mesh
        vertex_mesh.SetCellRearrangementThreshold(0.25);

        // Set the threshold element area to be small enough to prevent T2 swaps
        vertex_mesh.SetT2Threshold(0.001);

        /*
         * Calling ReMesh() once should result in a T1 swap occurring for nodes 4 and 5,
         * so that the mesh should now be as shown below.
         *  ______
         * |\ 2  /|
         * | \__/ |
         * |  \/  |
         * | 3 | 1|
         * |  /\  |
         * | / 0\ |
         * |/____\|
         */
        vertex_mesh.ReMesh();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], -0.2875, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.0875, 1e-3);

        // Test that each element contains the correct number of nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 5u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[3] = {0, 1, 4};
        unsigned node_indices_element_1[5] = {1, 2, 6, 5, 4};
        unsigned node_indices_element_2[4] = {2, 3, 7, 6};
        unsigned node_indices_element_3[5] = {0, 4, 5, 7, 3};
        unsigned node_indices_element_4[3] = {5, 6, 7};
        for (unsigned i=0; i<5; i++)
        {
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(i), node_indices_element_4[i]);
            }
            if (i < 4)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
        }

        // Reset the threshold element area to allow T2 swaps
        vertex_mesh.SetT2Threshold(0.1);

        /*
         * Calling ReMesh() once more should result in a T2 swap occurring for element 4,
         * so that the mesh should now be as shown below.
         *  ______
         * |\ 2  /|
         * | \  / |
         * |  \/  |
         * |3  | 1|
         * |  /\  |
         * | / 0\ |
         * |/____\|
         */

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // ... and give it to the killer
        T2SwapCellKiller<2> cell_killer(&cell_population);

        // Perform swaps
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test that the merged node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(8)->rGetLocation()[1], 0.2875/3.0, 1e-3);

        // Test that each element contains the correct number of nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 4u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned new_node_indices_element_0[3] = {0, 1, 4};
        unsigned new_node_indices_element_1[4] = {1, 2, 8, 4};
        unsigned new_node_indices_element_2[3] = {2, 3, 8};
        unsigned new_node_indices_element_3[4] = {0, 4, 8, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), new_node_indices_element_0[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), new_node_indices_element_2[i]);
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), new_node_indices_element_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), new_node_indices_element_3[i]);
        }
    }

    void TestT2SwapCellKillerException()
    {
        // Create a cell population whose type should not be used with a T2SwapCellKiller
        HoneycombMeshGenerator generator(4, 4, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> population(*p_mesh, cells);

        // Test that the correct exception is thrown if we try to use a T2SwapCellKiller with the population
        TS_ASSERT_THROWS_THIS(MAKE_PTR_ARGS(T2SwapCellKiller<2>, p_killer, (&population)),
            "A T2SwapCellKiller should only be used together with a VertexBasedCellPopulation.");
    }

    void TestArchivingOfT2SwapCellKiller()
    {
        // Set up singleton classes
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "T2swap_cell_killer.arch";
        ArchiveLocationInfo::SetMeshFilename("vertex_mesh_T2swap");

        {
            HoneycombVertexMeshGenerator generator(4,4);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            T2SwapCellKiller<2> cell_killer(&cell_population);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Serialize via pointer
            T2SwapCellKiller<2>* const p_cell_killer = &cell_killer;
            (*p_arch) << p_cell_killer;
        }

        {
            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            T2SwapCellKiller<2>* p_cell_killer;

            // Restore from the archive
            (*p_arch) >> p_cell_killer;

            TS_ASSERT(p_cell_killer != NULL);
            TS_ASSERT(p_cell_killer->GetCellPopulation() != NULL);

            // Tidy up
            delete p_cell_killer->mpCellPopulation;
            delete p_cell_killer;
        }
    }
};

#endif /*TESTT2SWAPCELLKILLER_HPP_*/
