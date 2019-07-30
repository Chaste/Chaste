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

#ifndef TESTMESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_
#define TESTMESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "CellId.hpp"
#include "FileComparison.hpp"
#include "ApoptoticCellProperty.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellPopulationAreaWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"


class TestMeshBasedCellPopulationWithGhostNodes : public AbstractCellBasedTestSuite
{
public:
    /*
     * Here we set up a test with 5 nodes, make a cell for each. We then set cell
     * 0 to be associated with node 1 instead of node 0, and Validate() throws an
     * exception. We then set node 0 to be a ghost node, and Validate() passes.
     */
    void TestValidateMeshBasedCellPopulationWithGhostNodes()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes()-1);

        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<cells.size(); i++)
        {
            cell_location_indices.push_back(i);
        }

        // Fails as the cell population constructor is not given the location indices
        // corresponding to real cells, so cannot work out which nodes are
        // ghost nodes
        std::vector<CellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(MeshBasedCellPopulationWithGhostNodes<2> dodgy_cell_population(mesh, cells_copy),
                "Node 4 does not appear to be a ghost node or have a cell associated with it");

        // Passes as the cell population constructor automatically works out which
        // cells are ghost nodes using the mesh and cell_location_indices
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(mesh, cells, cell_location_indices);

        // Here we set the ghost nodes to what they already are
        std::set<unsigned> ghost_node_indices;
        ghost_node_indices.insert(mesh.GetNumNodes()-1u);
        cell_population.SetGhostNodes(ghost_node_indices);

        // So validate passes at the moment
        cell_population.Validate();

        // Test GetCellUsingLocationIndex()

        TS_ASSERT_THROWS_NOTHING(cell_population.GetCellUsingLocationIndex(0)); // real cell
        TS_ASSERT_THROWS_THIS(cell_population.GetCellUsingLocationIndex(mesh.GetNumNodes()-1u),"Location index input argument does not correspond to a Cell"); // ghost node

        // Now we label a real cell's node as a ghost
        ghost_node_indices.insert(1u);

        // Validate detects this inconsistency
        TS_ASSERT_THROWS_THIS(cell_population.SetGhostNodes(ghost_node_indices),"Node 1 is labelled as a ghost node and has a cell attached");
    }

    // Test with ghost nodes, checking that the Iterator doesn't loop over ghost nodes
    void TestMeshBasedCellPopulationWithGhostNodesSetup()
    {
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices_set;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            node_indices.insert(p_mesh->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            location_indices_set.insert(location_indices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices_set.begin(), location_indices_set.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        std::vector<bool> is_ghost_node(p_mesh->GetNumNodes(), false);
        for (std::set<unsigned>::iterator it=ghost_node_indices.begin();
             it!=ghost_node_indices.end();
             it++)
        {
            is_ghost_node[*it] = true;
        }

        // Note in the following that some compilers (clang on recent OSX) can't cope with operator== overloading in CxxTest TS_ASSERT_EQUALS
        TS_ASSERT(cell_population.rGetGhostNodes() == is_ghost_node);
        //TS_ASSERT_EQUALS(cell_population.rGetGhostNodes(), is_ghost_node); // might not compile on recent clang

        // Test the GetGhostNodeIndices method
        std::set<unsigned> ghost_node_indices2 = cell_population.GetGhostNodeIndices();
        TS_ASSERT_EQUALS(ghost_node_indices, ghost_node_indices2);

        // Check the iterator doesn't loop over ghost nodes
        unsigned counter = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            TS_ASSERT( !(is_ghost_node[node_index]) );
            //TS_ASSERT_EQUALS(is_ghost_node[node_index], false); // might not compile on recent clang
            counter++;
        }

        TS_ASSERT_EQUALS(counter, cell_population.GetNumRealCells());

        // Check counter = num_nodes - num_ghost_nodes
        TS_ASSERT_EQUALS(counter + ghost_node_indices.size(), p_mesh->GetNumNodes());

        TS_ASSERT_EQUALS(cell_population.rGetGhostNodes().size(), p_mesh->GetNumNodes());

        // Test GetNeighbouringLocationIndices() method
        std::set<unsigned> expected_neighbours_of_cell_0;
        expected_neighbours_of_cell_0.insert(23);
        expected_neighbours_of_cell_0.insert(32);

        std::set<unsigned> neighbours_of_cell_0 = cell_population.GetNeighbouringLocationIndices(*(cell_population.Begin()));
        TS_ASSERT(neighbours_of_cell_0 == expected_neighbours_of_cell_0);
    }

    void TestCellPopulationIteratorWithNoCells()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create vector of cell location indices
        std::vector<unsigned> cell_location_indices;
        cell_location_indices.push_back(80);

        // Create a single cell
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());
        cells[0]->StartApoptosis();

        // Create a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(mesh, cells, cell_location_indices);

        // Iterate over cell population and check there is a single cell
        unsigned counter = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            counter++;
        }
        TS_ASSERT_EQUALS(counter, 1u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().empty(), false);

        // Increment simulation time and update cell population
        p_simulation_time->IncrementTimeOneStep();

        unsigned num_cells_removed = cell_population.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_cells_removed, 1u);

        cell_population.Update();

        // Iterate over cell population and check there are now no cells
        counter = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            counter++;
        }
        TS_ASSERT_EQUALS(counter, 0u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().empty(), true);
    }

    void TestAreaBasedVisocityOnAPeriodicMesh()
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 2;

        double scale_factor = 1.2;
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, scale_factor);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellPropertyRegistry::Instance()->Clear();
        RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

        // Set up cells
        std::vector<CellPtr> cells;
        cells.clear();
        unsigned num_cells = location_indices.empty() ? p_mesh->GetNumNodes() : location_indices.size();
        cells.reserve(num_cells);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            unsigned generation;
            double y = 0.0;

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }

            FixedG1GenerationalCellCycleModel* p_cell_cycle_model = new FixedG1GenerationalCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            double birth_time = -p_random_num_gen->ranf();
            if (y <= 0.3)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                generation = 0;
                birth_time *= typical_stem_cycle_time; // hours
            }
            else if (y < 2.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 1;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 3.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 2;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 4.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 3;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else
            {
                if (p_cell_cycle_model->CanCellTerminallyDifferentiate())
                {
                    p_cell->SetCellProliferativeType(p_diff_type);
                }
                else
                {
                    p_cell->SetCellProliferativeType(p_transit_type);
                }
                generation = 4;
                birth_time *= typical_transit_cycle_time; // hours
            }

            p_cell_cycle_model->SetGeneration(generation);
            p_cell->SetBirthTime(birth_time);

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                cells.push_back(p_cell);
            }
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        GeneralisedLinearSpringForce<2> linear_force;

        // It seems quite difficult to test this on a periodic mesh,
        // so just check the areas of all the cells are correct.

        cell_population.CreateVoronoiTessellation();

        TS_ASSERT_EQUALS(cell_population.GetVoronoiTessellation()->GetNumElements(), p_mesh->GetNumNodes());

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            double area = cell_population.GetVolumeOfVoronoiElement(node_index);
            TS_ASSERT_DELTA(area, sqrt(3.0)*scale_factor*scale_factor/2, 1e-6);
        }
    }

    void TestRemoveDeadCellsAndReMeshWithGhostNodes()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create vector of cell location indices
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=10; i<mesh.GetNumNodes(); i++)
        {
            if (i!=80)
            {
                cell_location_indices.push_back(i);
            }
        }

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());
        cells[27]->StartApoptosis();

        // Create a cell population, with some random ghost nodes
        MeshBasedCellPopulationWithGhostNodes<2> cell_population_with_ghost_nodes(mesh, cells, cell_location_indices);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);

        // Num real cells should be num_nodes (81) - num_ghosts (11) = 70
        TS_ASSERT_EQUALS(cell_population_with_ghost_nodes.GetNumRealCells(), 70u);
        TS_ASSERT_EQUALS(cell_population_with_ghost_nodes.GetNumAllCells(), 70u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed_with_ghost_nodes = cell_population_with_ghost_nodes.RemoveDeadCells();

        TS_ASSERT_EQUALS(num_removed_with_ghost_nodes, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
        TS_ASSERT_DIFFERS(cell_population_with_ghost_nodes.rGetCells().size(), cells.size()); // CellPopulation now copies cells
        TS_ASSERT_DIFFERS(cell_population_with_ghost_nodes.GetNumAllCells(), cells.size()); // CellPopulation now copies cells

        // Num real cells should be num_nodes (81) - num_ghosts (11) - 1 deleted node = 69
        TS_ASSERT_EQUALS(cell_population_with_ghost_nodes.GetNumRealCells(), 69u);
        TS_ASSERT_EQUALS(cell_population_with_ghost_nodes.rGetGhostNodes().size(), mesh.GetNumAllNodes());

        cell_population_with_ghost_nodes.Update();

        // For coverage
        NodeMap map(mesh.GetNumAllNodes());
        map.ResetToIdentity();
        cell_population_with_ghost_nodes.UpdateGhostNodesAfterReMesh(map);

        // Num real cells should be new_num_nodes (80) - num_ghosts (11)
        TS_ASSERT_EQUALS(cell_population_with_ghost_nodes.GetNumRealCells(), 69u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumAllNodes());
        TS_ASSERT_EQUALS(cell_population_with_ghost_nodes.rGetGhostNodes().size(), mesh.GetNumNodes());

        // Nodes 0-9 should not been renumbered so are still ghost nodes.
        // the ghost node at node 80 is now at 79 as node 27 was deleted..
        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            // True (ie should be a ghost) if i<10 or i==79, else false
            TS_ASSERT_EQUALS(cell_population_with_ghost_nodes.IsGhostNode(i), ((i<10)||(i==79)));
        }

        // Finally, check the cells node indices have updated

        // We expect the cell node indices to be {10,11,...,79}
        std::set<unsigned> expected_node_indices;
        for (unsigned i=0; i<cell_population_with_ghost_nodes.GetNumRealCells(); i++)
        {
            expected_node_indices.insert(i+10);
        }

        // Get actual cell node indices
        std::set<unsigned> node_indices_with_ghost_nodes;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population_with_ghost_nodes.Begin();
             cell_iter != cell_population_with_ghost_nodes.End();
             ++cell_iter)
        {
            // Record node index corresponding to cell
            unsigned node_index_with_ghost_nodes = cell_population_with_ghost_nodes.GetLocationIndexUsingCell(*cell_iter);
            node_indices_with_ghost_nodes.insert(node_index_with_ghost_nodes);
        }

        TS_ASSERT_EQUALS(node_indices_with_ghost_nodes, expected_node_indices);
    }

    void TestAddAndRemoveAndAddWithOutUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create vector of cell location indices
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=10; i<mesh.GetNumNodes(); i++)
        {
            if (i!=80)
            {
                cell_location_indices.push_back(i);
            }
        }

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());
        cells[27]->StartApoptosis();

        // Create a cell population, with some random ghost nodes
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(mesh, cells, cell_location_indices);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 70u);

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_stem_type);
        p_new_cell->SetBirthTime(0);

        c_vector<double,2> new_location;
        new_location[0] = 0.3433453454443;
        new_location[0] = 0.3435346344234;
        cell_population.AddCell(p_new_cell, cell_population.rGetCells().front()); // random choice of parent

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 71u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = cell_population.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_removed, 1u);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 70u);

        FixedG1GenerationalCellCycleModel* p_model2 = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell2(new Cell(p_state, p_model2));
        p_new_cell2->SetCellProliferativeType(p_stem_type);
        p_new_cell2->SetBirthTime(0);

        c_vector<double,2> new_location2;
        new_location2[0] = 0.6433453454443;
        new_location2[0] = 0.6435346344234;
        cell_population.AddCell(p_new_cell2, cell_population.rGetCells().front()); // random choice of parent

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 71u);
    }

    void TestSpringIterator2d()
    {
        // Set up expected results for the honeycomb mesh created below
        // the following are the edges which do not contain a ghost node
        std::set < std::set < unsigned > > expected_node_pairs;
        unsigned expected_node_pairs_array[] = {5,6,
                                                5,9,
                                                5,10,
                                                6,10,
                                                9,10 };

        for (unsigned i=0; i<10; i=i+2)
        {
            std::set < unsigned > node_pair;
            node_pair.insert(expected_node_pairs_array[i]);
            node_pair.insert(expected_node_pairs_array[i+1]);
            expected_node_pairs.insert(node_pair);
        }

        // Set up simple cell population with honeycomb mesh
        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;
        unsigned thickness_of_ghosts = 1;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghosts);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Check that we can iterate over the set of springs
        std::set< std::set< unsigned > > springs_visited;

        for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator=cell_population.SpringsBegin();
             spring_iterator!=cell_population.SpringsEnd();
             ++spring_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(spring_iterator.GetNodeA()->GetIndex());
            node_pair.insert(spring_iterator.GetNodeB()->GetIndex());

            TS_ASSERT_EQUALS(springs_visited.find(node_pair), springs_visited.end());
            springs_visited.insert(node_pair);

            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(spring_iterator.GetCellA()),
                             spring_iterator.GetNodeA()->GetIndex());

            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(spring_iterator.GetCellB()),
                             spring_iterator.GetNodeB()->GetIndex());
        }

         TS_ASSERT_EQUALS(springs_visited, expected_node_pairs);
    }

    // 3d test with some ghost nodes
    void TestSpringIterator3d()
    {
        // Create a simple mesh
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create vector of cell location indices
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=10; i<mesh.GetNumNodes(); i++)
        {
            cell_location_indices.push_back(i);
        }

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, cell_location_indices.size());

        // Create a cell population, with no ghost nodes at the moment
        MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, cell_location_indices);

        // Check that we can iterate over the set of springs
        std::set< std::set< unsigned > > springs_visited;

        for (MeshBasedCellPopulation<3>::SpringIterator spring_iterator=cell_population.SpringsBegin();
             spring_iterator!=cell_population.SpringsEnd();
             ++spring_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(spring_iterator.GetNodeA()->GetIndex());
            node_pair.insert(spring_iterator.GetNodeB()->GetIndex());

            TS_ASSERT_EQUALS(springs_visited.find(node_pair), springs_visited.end());
            springs_visited.insert(node_pair);

            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(spring_iterator.GetCellA()),
                             spring_iterator.GetNodeA()->GetIndex());

            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(spring_iterator.GetCellB()),
                             spring_iterator.GetNodeB()->GetIndex());
        }

        // Set up expected node pairs
        std::set< std::set<unsigned> > expected_node_pairs;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            Element<3,3>* p_element = mesh.GetElement(i);
            for (unsigned j=0; j<4; j++)
            {
                for (unsigned k=0; k<4; k++)
                {
                    unsigned node_A = p_element->GetNodeGlobalIndex(j);
                    unsigned node_B = p_element->GetNodeGlobalIndex(k);

                    // If nodeA or node_B are <10 they will have been labelled a ghost node above
                    if (node_A != node_B && node_A>=10 && node_B>=10)
                    {
                        std::set<unsigned> node_pair;
                        node_pair.insert(node_A);
                        node_pair.insert(node_B);

                        expected_node_pairs.insert(node_pair);
                    }
                }
            }
        }

        TS_ASSERT_EQUALS(springs_visited, expected_node_pairs);
    }

    void TestCellPopulationWritersIn3dWithGhostNodes()
    {

        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Resetting the Maximum cell Id to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 3D mesh with some ghost nodes
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0,  false,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1,  false,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2,  false,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3,  false,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4,  false,  0.5, 0.5, 0.5));
        nodes.push_back(new Node<3>(5,  true,  -1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(6,  true,   2.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(7,  true,   2.0,  2.0, -1.0));
        nodes.push_back(new Node<3>(8,  true,  -1.0,  2.0, -1.0));
        nodes.push_back(new Node<3>(9,  true,  -1.0, -1.0,  2.0));
        nodes.push_back(new Node<3>(10, true,   2.0, -1.0,  2.0));
        nodes.push_back(new Node<3>(11, true,   2.0,  2.0,  2.0));
        nodes.push_back(new Node<3>(12, true,  -1.0,  2.0,  2.0));
        MutableMesh<3,3> mesh(nodes);

        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<5; index++)
        {
            location_indices.push_back(index);
        }

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        cells[4]->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>()); // coverage
        TS_ASSERT_EQUALS(cells[4]->HasCellProperty<ApoptoticCellProperty>(), true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, location_indices);
        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "MeshBasedCellPopulationWithGhostNodes-3");

        // Test set methods
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddPopulationWriter<CellPopulationAreaWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.AddCellWriter<CellAncestorWriter>();

        // This method is usually called by Update()
        cell_population.CreateVoronoiTessellation();

        std::string output_directory = "TestCellPopulationWritersIn3dWithGhostNodes";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteResultsToFiles(output_directory);
        cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison( results_dir + "results.vizelements", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/results.vizelements").CompareFiles();
        FileComparison( results_dir + "results.viznodes", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/results.viznodes").CompareFiles();
        FileComparison( results_dir + "results.vizcelltypes", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/results.vizcelltypes").CompareFiles();
        FileComparison( results_dir + "results.vizlocationindices", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/results.vizlocationindices").CompareFiles();
        FileComparison( results_dir + "cellpopulationareas.dat", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/cellpopulationareas.dat").CompareFiles();
        FileComparison( results_dir + "cellareas.dat", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/cellareas.dat").CompareFiles();
        FileComparison( results_dir + "voronoi.dat", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/voronoi.dat").CompareFiles();
        FileComparison( results_dir + "results.vizancestors", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/results.vizancestors").CompareFiles();

        // Test the GetCellMutationStateCount function: there should only be healthy cells
        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 5u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 0u);

        // Test the GetCellProliferativeTypeCount() function - we should have 4 stem cells and 1 dead cell (for coverage)
        std::vector<unsigned> cell_types = cell_population.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 4u);
        TS_ASSERT_EQUALS(cell_types[0], 5u);
        TS_ASSERT_EQUALS(cell_types[1], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);
        TS_ASSERT_EQUALS(cell_types[3], 0u);

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare output with saved files of what they should look like
        FileComparison( results_dir + "results.parameters", "cell_based/test/data/TestCellPopulationWritersIn3dWithGhostNodes/results.parameters").CompareFiles();
    }

    void TestVoronoiAreasAndPerimetersWithGhostNodes()
    {
        // Create a small honeycomb mesh surrounded by a single layer of ghost nodes
        HoneycombMeshGenerator generator(2, 2, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create some cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create Voronoi tessellation (normally done in a simulation)
        cell_population.CreateVoronoiTessellation();

        // The Voronoi element corresponding to each real cell should be a regular hexagon
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            if (!cell_population.IsGhostNode(node_index))
            {
                TS_ASSERT_DELTA(cell_population.GetVolumeOfVoronoiElement(node_index), sqrt(3.0)/2, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetSurfaceAreaOfVoronoiElement(node_index), 6/sqrt(3.0), 1e-4);
            }
        }
    }

    void TestVoronoiGhostNodeLabelling2d()
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create 2D mesh with ghost nodes
        MutableMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(6, 6);

        c_vector<double, 2> mesh_centre = zero_vector<double>(2);
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            mesh_centre += mesh.GetNode(node_index)->rGetLocation() / mesh.GetNumNodes();
        }

        // Set up cells by iterating through the nodes
        std::vector<CellPtr> cells;
        std::vector<unsigned> location_indices;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Loop over nodes
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            // If this node is sufficiently close to the centre of the mesh, then create a cell for it
            c_vector<double, 2> node_location;
            node_location = mesh.GetNode(node_index)->rGetLocation();
            if (node_location(0) <= 3)
            {
                FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->SetBirthTime(-1.0);

                cells.push_back(p_cell);
                location_indices.push_back(node_index);
            }
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(mesh, cells, location_indices);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        // Create Voronoi tessellation
        cell_population.CreateVoronoiTessellation();

        // Check the correspondence between ghost nodes is correct
        for (VertexMesh<2,2>::VertexElementIterator elem_iter = cell_population.GetVoronoiTessellation()->GetElementIteratorBegin();
             elem_iter != cell_population.GetVoronoiTessellation()->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            unsigned node_index = cell_population.GetVoronoiTessellation()->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

            c_vector<double, 2> node_location;
            node_location = cell_population.GetNode(node_index)->rGetLocation();
            bool should_be_ghost_node = (node_location(0) > 3);

            TS_ASSERT_EQUALS(cell_population.IsGhostNode(node_index), should_be_ghost_node);
        }
    }

    void TestVoronoiGhostNodeLabelling3d()
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create 3D mesh with ghost nodes
        MutableMesh<3,3> mesh;
        mesh.ConstructCuboid(6, 6, 6);

        c_vector<double, 3> mesh_centre = zero_vector<double>(3);
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            mesh_centre += mesh.GetNode(node_index)->rGetLocation() / mesh.GetNumNodes();
        }

        // Set up cells by iterating through the nodes
        std::vector<CellPtr> cells;
        std::vector<unsigned> location_indices;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Loop over nodes
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            // If this node is sufficiently close to the centre of the mesh, then create a cell for it
            c_vector<double, 3> node_location;
            node_location = mesh.GetNode(node_index)->rGetLocation();
            if (node_location(0) <= 3)
            {
                FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->SetBirthTime(-1.0);

                cells.push_back(p_cell);
                location_indices.push_back(node_index);
            }
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, location_indices);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        // Create Voronoi tessellation
        cell_population.CreateVoronoiTessellation();

        // Check the correspondence between ghost nodes is correct
        for (VertexMesh<3,3>::VertexElementIterator elem_iter = cell_population.GetVoronoiTessellation()->GetElementIteratorBegin();
             elem_iter != cell_population.GetVoronoiTessellation()->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            unsigned node_index = cell_population.GetVoronoiTessellation()->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

            c_vector<double, 3> node_location;
            node_location = cell_population.GetNode(node_index)->rGetLocation();
            bool should_be_ghost_node = (node_location(0) > 3);

            TS_ASSERT_EQUALS(cell_population.IsGhostNode(node_index), should_be_ghost_node);
        }
    }

    void TestAddCellDataToPopulation()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a mesh-based cell population with ghost nodes
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        cell_population.SetDataOnAllCells("variable", 100.0);

        // Check that the data made it there and that copies of the data are independent
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("variable"), 100.0);
            cell_iter->GetCellData()->SetItem("variable", 1.0);
        }

        cell_population.SetDataOnAllCells("added variable", 200.0);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("added variable"), 200.0);
            cell_iter->GetCellData()->SetItem("added variable", 1.0);
        }

        std::vector<std::string> keys = cell_population.Begin()->GetCellData()->GetKeys();
        TS_ASSERT_EQUALS(keys.size(), 2u);
        TS_ASSERT_EQUALS(keys[0], "added variable");
        TS_ASSERT_EQUALS(keys[1], "variable");
    }

    void TestGetTetrahedralMeshForPdeModifier()
    {
        HoneycombMeshGenerator generator(2, 2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        TS_ASSERT_THROWS_THIS(cell_population.GetTetrahedralMeshForPdeModifier(),
            "Currently can't solve PDEs on meshes with ghost nodes");
    }
};

#endif /*TESTMESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_*/
