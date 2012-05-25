/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTNODEBASEDCELLPOPULATION_HPP_
#define TESTNODEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"

class TestNodeBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    template<unsigned DIM>
    void TestSimpleNodeBasedCellPopulation(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        TetrahedralMesh<DIM,DIM> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<DIM> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, DIM> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        unsigned num_cells = cells.size();

        // Create the cell population
        NodeBasedCellPopulation<DIM> node_based_cell_population(mesh, cells);

        // Test we have the correct numbers of nodes and cells
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(node_based_cell_population.rGetCells().size(), num_cells);
        TS_ASSERT_EQUALS(cells.size(), 0u);

        boost::shared_ptr<CellPropertyRegistry> p_population_registry = node_based_cell_population.GetCellPropertyRegistry();

        unsigned counter = 0;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = node_based_cell_population.Begin();
             cell_iter != node_based_cell_population.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            // Test that the cell property registry belonging to the population has made it into the cell's cell property collections.
            CellPropertyRegistry* p_cell_registry = cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry();
            TS_ASSERT_EQUALS(p_cell_registry,p_population_registry.get());
            counter++;
        }

        TS_ASSERT_EQUALS(counter, node_based_cell_population.GetNumRealCells());
    }

public:

    // Test construction, accessors and Iterator
    void TestNodeBasedCellPopulation1d2d3d() throw(Exception)
    {
        TestSimpleNodeBasedCellPopulation<1>("mesh/test/data/1D_0_to_1_10_elements");
        TestSimpleNodeBasedCellPopulation<2>("mesh/test/data/square_4_elements");
        TestSimpleNodeBasedCellPopulation<3>("mesh/test/data/cube_136_elements");
    }

    void TestOtherNodeBasedCellPopulationConstructor()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, generating_mesh.GetNumNodes());

        // Create the cell population
        unsigned num_cells = cells.size();
        std::vector<CellPtr> cells_copy(cells);
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Check box collection gets setup properly.
        c_vector<double, 4> domainSize;
        for (unsigned i=0; i< 2 ;i++)
        {
        	domainSize[2*i] = -10.0;
        	domainSize[2*i+1] = 10.0;
        }
        node_based_cell_population.SplitUpIntoBoxes(1.2, domainSize);

        TS_ASSERT_EQUALS(node_based_cell_population.GetBoxCollection()->GetNumBoxes(), 289u);

        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(node_based_cell_population.rGetCells().size(), num_cells);

        // For coverage, test that cell population constructor with 3rd argument locationIndices throws
        // an exception when the size of locationIndices does not equal the number of cells
        std::vector<unsigned> location_indices;
        location_indices.push_back(0);
        location_indices.push_back(1);
        location_indices.push_back(2);

        TS_ASSERT_THROWS_THIS(NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells_copy, location_indices),
                              "There is not a one-one correspondence between cells and location indices");
    }

    void TestValidateNodeBasedCellPopulation()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes()-1);

        // Fails as no cell corresponding to node 4
        std::vector<CellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(NodeBasedCellPopulation<2> cell_population(mesh, cells_copy),
                              "Node 4 does not appear to have a cell associated with it");

        // Add another cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        p_cell_cycle_model->SetCellProliferativeType(STEM);
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
        double birth_time = -4.0;
        p_cell->SetBirthTime(birth_time);
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 5u);

        // Throws exception as the cut-off length hasn't been set and has its default value of DBL_MAX
        TS_ASSERT_THROWS_THIS(cell_population.Update(),
                "NodeBasedCellPopulation cannot create boxes if the cut-off length has not been set - Call SetMechanicsCutOffLength on the CellPopulation ensuring it is larger than GetCutOffLength() on the force law");
        // Set Cut off length
        cell_population.SetMechanicsCutOffLength(1.2);
        cell_population.Update();

        std::set< std::pair<Node<2>*, Node<2>* > >& r_node_pairs = cell_population.rGetNodePairs();
        r_node_pairs.clear();


        cell_population.SetMechanicsCutOffLength(1e-3);
        cell_population.Update();
        TS_ASSERT(cell_population.rGetNodePairs().empty());

    }

    void TestAddCell()
    {
        // Create two nodes
        ChastePoint<2> point0;
        point0.rGetLocation()[0] = 0.0;
        point0.rGetLocation()[1] = 0.0;
        Node<2>* p_node0 = new Node<2>(0, point0, false);

        ChastePoint<2> point1;
        point1.rGetLocation()[0] = 1.0;
        point1.rGetLocation()[1] = 1.0;
        Node<2>* p_node1 = new Node<2>(1, point1, false);

        std::vector<Node<2>* > nodes;
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        // Create two cells
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model0 = new FixedDurationGenerationBasedCellCycleModel();
        p_model0->SetCellProliferativeType(STEM);
        CellPtr p_cell0(new Cell(p_state, p_model0));
        p_cell0->SetBirthTime(-1);

        FixedDurationGenerationBasedCellCycleModel* p_model1 = new FixedDurationGenerationBasedCellCycleModel();
        p_model1->SetCellProliferativeType(STEM);
        CellPtr p_cell1(new Cell(p_state, p_model1));
        p_cell1->SetBirthTime(-1);

        std::vector<CellPtr> cells;
        cells.push_back(p_cell0);
        cells.push_back(p_cell1);

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // For coverage
        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            TS_ASSERT_EQUALS(node_based_cell_population.IsParticle(i), false);
        }

        // Create a new cell, DON'T set the node index, set birth time=-1
        FixedDurationGenerationBasedCellCycleModel* p_model2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model2->SetCellProliferativeType(STEM);
        CellPtr p_cell2(new Cell(p_state, p_model2));
        p_cell2->SetBirthTime(-1);

        c_vector<double,2> cell2_location;
        cell2_location[0] = 2.0;
        cell2_location[1] = 2.0;

        node_based_cell_population.AddCell(p_cell2, cell2_location);

        // Tidy up
        delete p_node0;
        delete p_node1;
    }

    void TestSetNodeAndAddCell()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Test SetNode() by moving node 0 by a small amount
        AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
        c_vector<double,2> new_location = node_based_cell_population.GetLocationOfCellCentre(*cell_iter);
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        node_based_cell_population.SetNode(node_based_cell_population.GetLocationIndexUsingCell(*cell_iter), new_location_point);

        TS_ASSERT_DELTA(node_based_cell_population.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(node_based_cell_population.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        // Remove a cell so as to populate mDeletedNodeIndices (for coverage)
        node_based_cell_population.GetCellUsingLocationIndex(0)->Kill();
        node_based_cell_population.RemoveDeadCells();

        // Test AddNode
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = 0.51;
        new_point2.rGetLocation()[1] = 0.52;

        unsigned num_nodes = node_based_cell_population.GetNumNodes();
        Node<2>* p_node2 = new Node<2>(num_nodes, new_point2, false);
        unsigned new_node_index = node_based_cell_population.AddNode(p_node2);

        TS_ASSERT_EQUALS(new_node_index, 0u);
        TS_ASSERT_DELTA(node_based_cell_population.GetNode(0)->rGetLocation()[0], 0.51, 1e-12);
        TS_ASSERT_DELTA(node_based_cell_population.GetNode(0)->rGetLocation()[1], 0.52, 1e-12);

        // Test AddCell
        unsigned old_num_nodes = node_based_cell_population.GetNumNodes();
        unsigned old_num_cells = node_based_cell_population.rGetCells().size();

        // Create a new cell, DON'T set the node index, set birth time=-1
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetBirthTime(-1);

        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2.0;
        new_cell_location[1] = 2.0;

        node_based_cell_population.AddCell(p_cell, new_cell_location);

        // CellPopulation should have updated nodes and cells
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(node_based_cell_population.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), old_num_nodes);

        // Check the location of the new node
        TS_ASSERT_DELTA(node_based_cell_population.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(node_based_cell_population.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // Check the index of the new cell
        CellPtr& new_cell = node_based_cell_population.rGetCells().back();
        TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(new_cell), old_num_nodes);
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Make one cell start apoptosis
        cells[27]->StartApoptosis();

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 81u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = node_based_cell_population.RemoveDeadCells();

        node_based_cell_population.SetMechanicsCutOffLength(1.2);
        node_based_cell_population.Update(true);

        // Test that one cell has been removed
        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 80u);

        // Test that one node has been removed
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 80u);

        // Test that each cell's node index has been correctly updated
        unsigned index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
             cell_iter != node_based_cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(*cell_iter), index);
            index++;
        }
    }

    void TestAddAndRemoveAndAddWithOutRemovingDeletedNodes()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh of the domain [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            mesh.SetCellRadius(i, 0.1);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Make one cell start apoptosis
        cells[27]->StartApoptosis();

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 81u);

        // Test GetNeighbouringNodeIndices() method
        std::set<unsigned> node_50_neighbours = node_based_cell_population.GetNeighbouringNodeIndices(50);

        std::set<unsigned> expected_node_50_neighbours;
        expected_node_50_neighbours.insert(10);
        expected_node_50_neighbours.insert(18);
        expected_node_50_neighbours.insert(27);
        expected_node_50_neighbours.insert(34);
        expected_node_50_neighbours.insert(48);
        expected_node_50_neighbours.insert(49);
        expected_node_50_neighbours.insert(64);
        expected_node_50_neighbours.insert(66);

        TS_ASSERT_EQUALS(node_50_neighbours.size(), expected_node_50_neighbours.size());
        TS_ASSERT_EQUALS(node_50_neighbours, expected_node_50_neighbours);

        // Add a cell to the cell population
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetBirthTime(0);
        c_vector<double,2> new_location;
        new_location[0] = 0.3433453454443;
        new_location[1] = 0.3435346344234;
        node_based_cell_population.AddCell(p_new_cell, new_location);

        // Test that the numbers of nodes and cells has been updated
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 82u);

        p_simulation_time->IncrementTimeOneStep();

        // Test that the apoptotic cell has been removed
        unsigned num_removed = node_based_cell_population.RemoveDeadCells();
        node_based_cell_population.SetMechanicsCutOffLength(1.2);
        node_based_cell_population.Update();

        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 81u);

        // Add another cell to the cell population
        FixedDurationGenerationBasedCellCycleModel* p_model2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model2->SetCellProliferativeType(STEM);
        CellPtr p_new_cell2(new Cell(p_state, p_model2));
        p_new_cell2->SetBirthTime(0);

        c_vector<double,2> new_location2;
        new_location2[0] = 0.6433453454443;
        new_location2[1] = 0.6435346344234;
        node_based_cell_population.AddCell(p_new_cell2, new_location2);

        // Test that the numbers of nodes and cells has been updated
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 82u);
    }

    void TestSettingCellAncestors() throw (Exception)
    {
        // Create a small node-based cell population
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Test that the cell population makes each cell fix the corresponding node index as its ancestor
        cell_population.SetCellAncestorsToLocationIndices();

        unsigned counter = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetAncestor(), cell_population.GetLocationIndexUsingCell(*cell_iter));
            counter ++;
        }
        TS_ASSERT_EQUALS(counter, 5u);

        // Test that we can recover the remaining number of ancestors
        std::set<unsigned> remaining_ancestors = cell_population.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 5u);

        // Reallocate ancestors
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Set all cells to have the same ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (1u));
            cell_iter->SetAncestor(p_cell_ancestor);
        }

        // Test that the cell population now shares a common ancestor
        remaining_ancestors = cell_population.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 1u);
    }

    void TestGetLocationOfCellCentreAndGetWidth() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Loop over nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
             cell_iter != node_based_cell_population.End();
             ++cell_iter)
        {
            // Record node location
            c_vector<double, 2> node_location = node_based_cell_population.GetLocationOfCellCentre(*cell_iter);

            // Test GetLocationOfCellCentre()
            TS_ASSERT_DELTA(node_location[0], node_based_cell_population.GetLocationOfCellCentre(*cell_iter)[0], 1e-9);
            TS_ASSERT_DELTA(node_location[1], node_based_cell_population.GetLocationOfCellCentre(*cell_iter)[1], 1e-9);
        }

        // Test GetWidth() method
        double width_x = node_based_cell_population.GetWidth(0);
        TS_ASSERT_DELTA(width_x, 1.0, 1e-6);

        double width_y = node_based_cell_population.GetWidth(1);
        TS_ASSERT_DELTA(width_y, 1.0, 1e-6);

        c_vector<double, 2> size_of_pop = node_based_cell_population.GetSizeOfCellPopulation();
        TS_ASSERT_DELTA(size_of_pop[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(size_of_pop[0], 0.5, 1e-4);

    }

    void TestNodeBasedCellPopulationOutputWriters2d()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // For coverage of GetVolumeOfCell()
        for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
             cell_iter != node_based_cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(node_based_cell_population.GetVolumeOfCell(*cell_iter), M_PI*pow(node_based_cell_population.rGetMesh().GetCellRadius(node_based_cell_population.GetNodeCorrespondingToCell(*cell_iter)->GetIndex()),2), 1e-6);
        }

        // For coverage of WriteResultsToFiles()
        boost::shared_ptr<AbstractCellProperty> p_state(node_based_cell_population.GetCellPropertyRegistry()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc1(node_based_cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(node_based_cell_population.GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(node_based_cell_population.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(node_based_cell_population.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(node_based_cell_population.GetCellPropertyRegistry()->Get<CellLabel>());

        node_based_cell_population.GetCellUsingLocationIndex(0)->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);
        node_based_cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        node_based_cell_population.GetCellUsingLocationIndex(1)->GetCellCycleModel()->SetCellProliferativeType(DIFFERENTIATED);
        node_based_cell_population.GetCellUsingLocationIndex(1)->SetMutationState(p_apc1);
        node_based_cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);
        node_based_cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);
        node_based_cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_apoptotic_state);

        node_based_cell_population.rGetMesh().SetCellRadius(0, 3.0);
        node_based_cell_population.rGetMesh().SetCellRadius(1, 2.0);

        TS_ASSERT_EQUALS(node_based_cell_population.GetOutputCellIdData(), false);
        node_based_cell_population.SetOutputCellIdData(true);
        TS_ASSERT_EQUALS(node_based_cell_population.GetOutputCellIdData(), true);

        // Coverage of writing CellData to VTK
        for (unsigned var=0; var<2; var++)
        {
            for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
                 cell_iter != node_based_cell_population.End();
                 ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem(var, (double) 3.0*var);
            }
        }

        // Test set methods
        std::string output_directory = "TestNodeBasedCellPopulationWriters2d";
        OutputFileHandler output_file_handler(output_directory, false);

        node_based_cell_population.SetOutputCellMutationStates(true);
        node_based_cell_population.SetOutputCellProliferativeTypes(true);
        node_based_cell_population.SetOutputCellCyclePhases(true);
        node_based_cell_population.SetOutputCellAges(true);
        node_based_cell_population.SetOutputCellVolumes(true);

        node_based_cell_population.SetOutputCellAncestors(true);
        node_based_cell_population.SetCellAncestorsToLocationIndices();

        TS_ASSERT_THROWS_NOTHING(node_based_cell_population.CreateOutputFiles(output_directory, false));

        node_based_cell_population.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(node_based_cell_population.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     cell_based/test/data/TestNodeBasedCellPopulationWriters2d/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     cell_based/test/data/TestNodeBasedCellPopulationWriters2d/cellages.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat         cell_based/test/data/TestNodeBasedCellPopulationWriters2d/cellareas.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = node_based_cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 2u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 1u);

        // Test the GetCellProliferativeTypeCount function
        std::vector<unsigned> cell_types = node_based_cell_population.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 3u);
        TS_ASSERT_EQUALS(cell_types[0], 3u);
        TS_ASSERT_EQUALS(cell_types[1], 1u);
        TS_ASSERT_EQUALS(cell_types[2], 1u);

        // Test the Get and set MechanicsCutOfLengthMethods
        node_based_cell_population.SetMechanicsCutOffLength(1.5);
        TS_ASSERT_DELTA(node_based_cell_population.GetMechanicsCutOffLength(),1.5, 1e-9);

        // For coverage
        TS_ASSERT_THROWS_NOTHING(node_based_cell_population.WriteResultsToFiles());

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        node_based_cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare output with saved files of what they should look like
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters         cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.parameters").c_str()), 0);
    }

    void TestNodeBasedCellPopulationOutputWriters3d()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "NodeBasedCellPopulation-3");

        // Test set methods
        std::string output_directory = "TestNodeBasedCellPopulationWriters3d";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.SetOutputCellVolumes(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellCyclePhases(true);

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellAncestors(true);

        TS_ASSERT_THROWS_NOTHING(cell_population.CreateOutputFiles(output_directory, false));

        cell_population.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(cell_population.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     cell_based/test/data/TestNodeBasedCellPopulationWriters3d/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     cell_based/test/data/TestNodeBasedCellPopulationWriters3d/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     cell_based/test/data/TestNodeBasedCellPopulationWriters3d/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     cell_based/test/data/TestNodeBasedCellPopulationWriters3d/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     cell_based/test/data/TestNodeBasedCellPopulationWriters3d/cellages.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat         cell_based/test/data/TestNodeBasedCellPopulationWriters3d/cellareas.dat").c_str()), 0);
    }

    void TestWritingCellCyclePhases()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        cells[0]->SetBirthTime(-23.5);
        cells[1]->SetBirthTime(-0.5);
        cells[2]->SetBirthTime(-1.5);
        cells[3]->SetBirthTime(-15.5);
        cells[4]->SetBirthTime(-23.5);
        cells[0]->GetCellCycleModel()->SetCellProliferativeType(DIFFERENTIATED);

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        TS_ASSERT_EQUALS(node_based_cell_population.GetIdentifier(), "NodeBasedCellPopulation-2");

        // Loop over cells to run to time 0
        for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
             cell_iter != node_based_cell_population.End();
             ++cell_iter)
        {
            cell_iter->ReadyToDivide();
        }

        std::string output_directory = "TestWritingCellCyclePhases";
        OutputFileHandler output_file_handler(output_directory, false);

        node_based_cell_population.SetOutputCellCyclePhases(true);

        node_based_cell_population.CreateOutputFiles(output_directory, false);
        node_based_cell_population.WriteResultsToFiles();
        node_based_cell_population.CloseOutputFiles();

        // Test the rGetCellCyclePhaseCount function
        std::vector<unsigned> cell_cycle_phases = node_based_cell_population.rGetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phases[0], 1u);
        TS_ASSERT_EQUALS(cell_cycle_phases[1], 3u);
        TS_ASSERT_EQUALS(cell_cycle_phases[2], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[3], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[4], 1u);
    }

    void TestArchivingCellPopulation() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "node_based_cell_population.arch";
        ArchiveLocationInfo::SetMeshFilename("node_based_cell_population_mesh");

        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            TetrahedralMesh<2,2> generating_mesh;
            generating_mesh.ConstructFromMeshReader(mesh_reader);

            // Convert this to a NodesOnlyMesh
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(generating_mesh);

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            // Create a cell population
            NodeBasedCellPopulation<2>* const p_cell_population = new NodeBasedCellPopulation<2>(mesh, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                cell_iter != p_cell_population->End();
                ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            p_cell_population->SetMechanicsCutOffLength(1.5);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        {
            // Need to set up time
            unsigned num_steps = 10;

            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            NodeBasedCellPopulation<2>* p_cell_population;

            // Restore the cell population
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_cell_population;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;

            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the cell population has been restored
            TS_ASSERT_EQUALS(p_cell_population->rGetCells().size(), 5u);

            TS_ASSERT_DELTA(p_cell_population->GetMechanicsCutOffLength(), 1.5, 1e-6);

            // Check number of nodes
            TS_ASSERT_EQUALS(p_cell_population->GetNumNodes(), 5u);

            // Check some node positions
            TS_ASSERT_EQUALS(p_cell_population->GetNode(3)->GetIndex(), 3u);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(4)->GetIndex(), 4u);

            TS_ASSERT_DELTA(p_cell_population->GetNode(3)->rGetLocation()[0], 0.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_population->GetNode(3)->rGetLocation()[1], 1.0, 1e-9);
            TS_ASSERT_DELTA(p_cell_population->GetNode(4)->rGetLocation()[0], 0.5, 1e-9);
            TS_ASSERT_DELTA(p_cell_population->GetNode(4)->rGetLocation()[1], 0.5, 1e-9);

            // Check the member variables have been restored
            TS_ASSERT_DELTA(p_cell_population->GetMechanicsCutOffLength(), 1.5, 1e-9);

            // Tidy up
            delete p_cell_population;
        }
    }
};

#endif /*TESTNODEBASEDCELLPOPULATION_HPP_*/
