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

#ifndef TESTNODEBASEDCELLPOPULATION_HPP_
#define TESTNODEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "ArchiveOpener.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellAncestor.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "FixedCentreBasedDivisionRule.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "TetrahedralMesh.hpp"
#include "TransitCellProliferativeType.hpp"
#include "TrianglesMeshReader.hpp"
#include "UblasCustomFunctions.hpp"
#include "WildTypeCellMutationState.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellPopulationAreaWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestNodeBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    template<unsigned DIM>
    void DimensionTestSimpleNodeBasedCellPopulation(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        TetrahedralMesh<DIM,DIM> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<DIM> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, DIM> cells_generator;
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
            TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(*cell_iter), counter*PetscTools::GetNumProcs() + PetscTools::GetMyRank());

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
    void TestNodeBasedCellPopulation1d2d3d()
    {
        DimensionTestSimpleNodeBasedCellPopulation<1>("mesh/test/data/1D_0_to_1_10_elements");
        DimensionTestSimpleNodeBasedCellPopulation<2>("mesh/test/data/square_4_elements");
        DimensionTestSimpleNodeBasedCellPopulation<3>("mesh/test/data/cube_136_elements");
    }

    void TestOtherNodeBasedCellPopulationConstructor()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.2);
        mesh.SetCalculateNodeNeighbours(false);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the cell population
        unsigned num_cells = cells.size();
        std::vector<CellPtr> cells_copy(cells);
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Update the population
        node_based_cell_population.Update();

        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(node_based_cell_population.rGetCells().size(), num_cells);

        // For coverage, test that cell population constructor with 3rd argument locationIndices throws
        // an exception when the size of locationIndices does not equal the number of cells
        std::vector<unsigned> location_indices;
        location_indices.push_back(0);
        location_indices.push_back(1);
        location_indices.push_back(2);

        TS_ASSERT_THROWS_THIS(NodeBasedCellPopulation<2> another_node_based_cell_population(mesh, cells_copy, location_indices),
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
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.2);

        // Create cells
        unsigned num_nodes = (mesh.GetNumNodes() > 0) ? mesh.GetNumNodes() - 1 : 0;

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, num_nodes);

        // Fails as no cell corresponding to a node
        if (mesh.GetNumNodes() > 0)
        {
            std::vector<CellPtr> cells_copy(cells);
            TS_ASSERT_THROWS_CONTAINS(NodeBasedCellPopulation<2> cell_population(mesh, cells_copy),
                                  "does not appear to have a cell associated with it");

            // Add another cell
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            FixedG1GenerationalCellCycleModel* p_cell_cycle_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            double birth_time = -4.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), mesh.GetNumNodes());

        cell_population.Update();

        std::vector< std::pair<Node<2>*, Node<2>* > >& r_node_pairs = cell_population.rGetNodePairs();
        r_node_pairs.clear();

        // Set a new cut-off
        if (PetscTools::IsSequential()) // This causes nodes to jump to different process in parallel
        {
            mesh.Clear();
            mesh.ConstructNodesWithoutMesh(generating_mesh, 1e-3);

            cell_population.Update();
            TS_ASSERT(cell_population.rGetNodePairs().empty());
        }

        // For coverage, test that GetDefaultTimeStep() returns the correct value
        TS_ASSERT_DELTA(cell_population.GetDefaultTimeStep(), 1.0/120.0, 1e-6);
    }

    void TestUpdatingCellLocationMapOnDelete()
    {
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);

        assert(cells.empty());
        cells.resize(mesh.GetNumNodes());
        std::copy(node_based_cell_population.rGetCells().begin(), node_based_cell_population.rGetCells().end(), cells.begin());

        if (PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(cells[0]), 0u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(cells[1]), PetscTools::GetNumProcs());

            cells[0]->Kill();
            node_based_cell_population.RemoveDeadCells();
        }

        node_based_cell_population.Update();

        if (PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 1u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(cells[1]), PetscTools::GetNumProcs());
        }

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestAddCell()
    {
        // Create two nodes (for coverage, give one some node attributes)
        ChastePoint<2> point0;
        point0.rGetLocation()[0] = 0.0;
        point0.rGetLocation()[1] = 0.0;
        Node<2>* p_node0 = new Node<2>(0, point0, false);
        p_node0->AddNodeAttribute(0.0);
        std::vector<double>& attributes = p_node0->rGetNodeAttributes();
        attributes.resize(2);
        attributes[0] = 6.23;
        attributes[1] = 5.91;

        ChastePoint<2> point1;
        point1.rGetLocation()[0] = 1.0;
        point1.rGetLocation()[1] = 1.0;
        Node<2>* p_node1 = new Node<2>(1, point1, false);
        TS_ASSERT_EQUALS(p_node1->HasNodeAttributes(), false);

        std::vector<Node<2>* > nodes;
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create two cells
        std::vector<CellPtr> cells;
        if (PetscTools::AmMaster())
        {
            mesh.GetNode(0)->SetRadius(0.1);
            mesh.GetNode(PetscTools::GetNumProcs())->SetRadius(0.2);

            boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

            FixedG1GenerationalCellCycleModel* p_model0 = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell0(new Cell(p_state, p_model0));
            p_cell0->SetCellProliferativeType(p_stem_type);
            p_cell0->SetBirthTime(-1);

            FixedG1GenerationalCellCycleModel* p_model1 = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell1(new Cell(p_state, p_model1));
            p_cell1->SetCellProliferativeType(p_stem_type);
            p_cell1->SetBirthTime(-1);

            cells.push_back(p_cell0);
            cells.push_back(p_cell1);
        }

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // For coverage
        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
             node_iter != mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            TS_ASSERT_EQUALS(node_iter->IsParticle(), false);
        }

        if (PetscTools::AmMaster())
        {
            boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

            // Create a new cell, DON'T set the node index, set birth time=-1
            FixedG1GenerationalCellCycleModel* p_model2 = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell2(new Cell(p_state, p_model2));
            p_cell2->SetCellProliferativeType(p_stem_type);
            p_cell2->SetBirthTime(-1);

            c_vector<double,2> cell2_location;
            cell2_location[0] = 0.9;
            cell2_location[1] = 1.4;

            node_based_cell_population.AddCell(p_cell2, node_based_cell_population.GetCellUsingLocationIndex(0));

            // Check the radius and node attributes associated with cell 0 are correct
            AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
            TS_ASSERT_DELTA(node_iter->GetRadius(), 0.1, 1e-6);
            TS_ASSERT_EQUALS(node_iter->HasNodeAttributes(), true);
            TS_ASSERT_EQUALS(node_iter->rGetNodeAttributes().size(), 2u);
            TS_ASSERT_DELTA(node_iter->rGetNodeAttributes()[0], 6.23, 1e-4);
            TS_ASSERT_DELTA(node_iter->rGetNodeAttributes()[1], 5.91, 1e-4);

            // Check the radius of cell 1 is correct and it has no associated node attributes
            TS_ASSERT_DELTA((++node_iter)->GetRadius(), 0.2, 1e-6);

            /*
             * Note: since the radius of each node is set to 0.5 in
             * NodesOnlyMesh::ConstructNodesWithoutMesh(), this means
             * that every node in a NodeBasedCellPopulaton has called
             * ConstructNodeAttributes(); however, rGetNodeAttributes()
             * will return an empty vector unless any attributes have
             * been set by the user.
             */
            TS_ASSERT_EQUALS(node_iter->HasNodeAttributes(), true);
            TS_ASSERT_EQUALS(node_iter->rGetNodeAttributes().size(), 0u);

            // Check the radius and node attributes associated with cell 2 are correct (cell 0 divided into 0 and 2)
            TS_ASSERT_DELTA((++node_iter)->GetRadius(), 0.1, 1e-6);
            TS_ASSERT_EQUALS(node_iter->HasNodeAttributes(), true);
            TS_ASSERT_EQUALS(node_iter->rGetNodeAttributes().size(), 2u);
            TS_ASSERT_DELTA(node_iter->rGetNodeAttributes()[0], 6.23, 1e-4);
            TS_ASSERT_DELTA(node_iter->rGetNodeAttributes()[1], 5.91, 1e-4);
        }

        // Avoid memory leak
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
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Test SetNode() by moving node 0 by a small amount
        if (PetscTools::AmMaster())
        {
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

            TS_ASSERT_EQUALS(mesh.SolveNodeMapping(new_node_index), 0u);
            TS_ASSERT_DELTA(node_based_cell_population.GetNode(0)->rGetLocation()[0], 0.51, 1e-12);
            TS_ASSERT_DELTA(node_based_cell_population.GetNode(0)->rGetLocation()[1], 0.52, 1e-12);

            // Test AddCell
            unsigned old_num_nodes = node_based_cell_population.GetNumNodes();
            unsigned old_num_cells = node_based_cell_population.rGetCells().size();

            // Create a new cell, DON'T set the node index, set birth time=-1
            boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->SetBirthTime(-1);

            c_vector<double,2> new_cell_location;
            new_cell_location[0] = 1.4;
            new_cell_location[1] = 1.4;

            typedef FixedCentreBasedDivisionRule<2,2> FixedRule;
            MAKE_PTR_ARGS(FixedRule, p_div_rule, (new_cell_location));
            node_based_cell_population.SetCentreBasedDivisionRule(p_div_rule);

            CellPtr p_parent_cell = node_based_cell_population.GetCellUsingLocationIndex(PetscTools::GetNumProcs());

            CellPtr p_child_cell = node_based_cell_population.AddCell(p_cell, p_parent_cell);

            // CellPopulation should have updated nodes and cells
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), old_num_nodes+1);
            TS_ASSERT_EQUALS(node_based_cell_population.rGetCells().size(), old_num_cells+1);
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), old_num_nodes);

            // Check the location of the new node
            unsigned location_index_new_cell = node_based_cell_population.GetLocationIndexUsingCell(p_child_cell);
            TS_ASSERT_DELTA(node_based_cell_population.GetNode(location_index_new_cell)->rGetLocation()[0], 1.4, 1e-12);
            TS_ASSERT_DELTA(node_based_cell_population.GetNode(location_index_new_cell)->rGetLocation()[1], 1.4, 1e-12);

            // Check the index of the new cell
            CellPtr& new_cell = node_based_cell_population.rGetCells().back();
            TS_ASSERT_EQUALS(node_based_cell_population.GetLocationIndexUsingCell(new_cell), PetscTools::GetNumProcs()*(old_num_nodes + 1));
        }
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
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.2);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        unsigned num_removed = 0;

        if (PetscTools::AmMaster())
        {
            // Make one cell start apoptosis
            cells.resize(mesh.GetNumNodes());
            std::copy(node_based_cell_population.rGetCells().begin(), node_based_cell_population.rGetCells().end(), cells.begin());
            cells[27]->StartApoptosis();

            // Test we have the right numbers of nodes and cells
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 81u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 81u);

            p_simulation_time->IncrementTimeOneStep();

            num_removed = node_based_cell_population.RemoveDeadCells();
        }

        // Update will deadlock if called on a single process.
        node_based_cell_population.Update(true);

        if (PetscTools::AmMaster())
        {
            // Test that one cell has been removed
            TS_ASSERT_EQUALS(num_removed, 1u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 80u);

            // Test that one node has been removed
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 80u);

            // Test that each cell's location index is correct.
            unsigned index = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
                 cell_iter != node_based_cell_population.End();
                 ++cell_iter)
            {
                unsigned global_index = node_based_cell_population.GetLocationIndexUsingCell(*cell_iter);
                unsigned local_index = mesh.SolveNodeMapping(global_index);
                TS_ASSERT_EQUALS(local_index, index);
                index++;
            }
        }
    }

    void TestAddAndRemoveAndAddWithOutRemovingDeletedNodesSmallCutOff()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh of the domain [0,1]x[0,1]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1e-1);
        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->SetRadius(0.1);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Only do this is we have some local nodes in the mesh...
        if (mesh.GetNumNodes() > PetscTools::GetMyRank())
        {
            // Make one cell start apoptosis
            node_based_cell_population.GetCellUsingLocationIndex(PetscTools::GetMyRank())->StartApoptosis();

            // Test we have the right numbers of nodes and cells
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), mesh.GetNumNodes());

            // Test GetNeighbouringNodeIndices() method
            TS_ASSERT_THROWS_THIS(node_based_cell_population.GetNeighbouringNodeIndices(PetscTools::GetMyRank()), "mNodeNeighbours not set up. Call Update() before GetNeighbouringNodeIndices()");
        }

        // Update
        node_based_cell_population.Update();

        // Only do this is we have some local nodes in the mesh...
        if (mesh.GetNumNodes() > PetscTools::GetMyRank())
        {

            TS_ASSERT_THROWS_CONTAINS(node_based_cell_population.GetNeighbouringNodeIndices(PetscTools::GetMyRank()), "mpNodesOnlyMesh::mMaxInteractionDistance is smaller than twice the radius of cell");
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
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.2);

        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->SetRadius(0.1);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Check that the master owns all the nodes (the box size is bigger than the mesh)
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        }
        else
        {
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), 0u);
        }
        // Make one cell start apoptosis
        if (PetscTools::AmMaster())
        {
            node_based_cell_population.GetCellUsingLocationIndex(0)->StartApoptosis();
        }

        node_based_cell_population.Update();

        // Note that the cells are only on the master process
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(node_based_cell_population.IsPdeNodeAssociatedWithNonApoptoticCell(0), false);
            unsigned next_cell_id = PetscTools::GetNumProcs(); // Parallel code uses modular arithmetic to assign cell IDs
            TS_ASSERT_EQUALS(node_based_cell_population.IsPdeNodeAssociatedWithNonApoptoticCell(next_cell_id), true);
        }

        unsigned num_removed;
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

        if (PetscTools::AmMaster())
        {
            std::set<unsigned> node_50_neighbours = node_based_cell_population.GetNeighbouringNodeIndices(50 * PetscTools::GetNumProcs());

            std::set<unsigned> expected_node_50_neighbours;
            expected_node_50_neighbours.insert(10 * PetscTools::GetNumProcs());
            expected_node_50_neighbours.insert(18 * PetscTools::GetNumProcs());
            expected_node_50_neighbours.insert(27 * PetscTools::GetNumProcs());
            expected_node_50_neighbours.insert(34 * PetscTools::GetNumProcs());
            expected_node_50_neighbours.insert(48 * PetscTools::GetNumProcs());
            expected_node_50_neighbours.insert(49 * PetscTools::GetNumProcs());
            expected_node_50_neighbours.insert(64 * PetscTools::GetNumProcs());
            expected_node_50_neighbours.insert(66 * PetscTools::GetNumProcs());

            TS_ASSERT_EQUALS(node_50_neighbours.size(), expected_node_50_neighbours.size());
            TS_ASSERT_EQUALS(node_50_neighbours, expected_node_50_neighbours);

            // Add a cell to the cell population

            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_new_cell(new Cell(p_state, p_model));
            p_new_cell->SetCellProliferativeType(p_stem_type);
            p_new_cell->SetBirthTime(0);
            c_vector<double,2> new_location;
            new_location[0] = 0.3433453454443;
            new_location[1] = 0.3435346344234;

            CellPtr p_parent_cell = node_based_cell_population.GetCellUsingLocationIndex(0u);

            node_based_cell_population.AddCell(p_new_cell, p_parent_cell);

            // Test that the numbers of nodes and cells has been updated
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 82u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 82u);

            p_simulation_time->IncrementTimeOneStep();
        }

        // Test that the apoptotic cell has been removed
        num_removed = node_based_cell_population.RemoveDeadCells();

        node_based_cell_population.Update();

        if (PetscTools::AmMaster())
        {
            TS_ASSERT_EQUALS(num_removed, 1u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 81u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 81u);

            // Add another cell to the cell population
            FixedG1GenerationalCellCycleModel* p_model2 = new FixedG1GenerationalCellCycleModel();
            CellPtr p_new_cell2(new Cell(p_state, p_model2));
            p_new_cell2->SetCellProliferativeType(p_stem_type);
            p_new_cell2->SetBirthTime(0);

            c_vector<double,2> new_location2;
            new_location2[0] = 0.6433453454443;
            new_location2[1] = 0.6435346344234;

            CellPtr p_parent_cell = node_based_cell_population.GetCellUsingLocationIndex(PetscTools::GetNumProcs());

            node_based_cell_population.AddCell(p_new_cell2, p_parent_cell); // Use same parent cell

            // Test that the numbers of nodes and cells has been updated
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), 82u);
            TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), 82u);
        }
    }

    void TestGetNeighbouringNodeIndices()
    {
        EXIT_IF_PARALLEL;    // Doesn't work in parallel yet until halo nodes are updated (#2364)

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a small node-based cell population
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 0.1);

        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->SetRadius(0.55);
        }

        // Make cell 0 smaller
        if (PetscTools::AmMaster())
        {
            mesh.GetNode(0)->SetRadius(0.1);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), mesh.GetNumNodes());

        // Test GetNeighbouringNodeIndices() method
        if (PetscTools::AmMaster())
        {
            // Coverage
            TS_ASSERT_THROWS_THIS(node_based_cell_population.GetNeighbouringNodeIndices(0), "mNodeNeighbours not set up. Call Update() before GetNeighbouringNodeIndices()");
        }

        node_based_cell_population.Update();

        if (PetscTools::AmMaster())
        {
            std::cout << ((node_based_cell_population.GetNode(0))->rGetNeighbours()).size() << std::endl;
            TS_ASSERT_THROWS_THIS(node_based_cell_population.GetNeighbouringNodeIndices(0), "mpNodesOnlyMesh::mMaxInteractionDistance is smaller than twice the radius of cell 0 (0.1) so interactions may be missed. Make the cut-off larger to avoid errors.");
        }

        mesh.Clear();
        mesh.ConstructNodesWithoutMesh(generating_mesh, 0.5);

        // Re set the radii
        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->SetRadius(0.55);
        }

        // Make cell 0 smaller
        if (PetscTools::AmMaster())
        {
            mesh.GetNode(0)->SetRadius(0.1);
        }
        node_based_cell_population.Update();

        TS_ASSERT_THROWS_THIS(node_based_cell_population.GetNeighbouringNodeIndices(0), "mpNodesOnlyMesh::mMaxInteractionDistance is smaller than the sum of radius of cell 0 (0.1) and cell 4 (0.55). Make the cut-off larger to avoid errors.");

        mesh.Clear();
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.1);

        // Re set the radii
        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->SetRadius(0.55);
        }

        // Make cell 0 smaller
        if (PetscTools::AmMaster())
        {
            mesh.GetNode(0)->SetRadius(0.1);
        }

        node_based_cell_population.Update();

        //Test Corner Node should have no neighbours as small cell
        if (PetscTools::AmMaster())
        {
            std::set<unsigned> node_0_neighbours = node_based_cell_population.GetNeighbouringNodeIndices(0);

            std::set<unsigned> expected_node_0_neighbours;

            TS_ASSERT_EQUALS(node_0_neighbours.size(), 0u);
            TS_ASSERT_EQUALS(node_0_neighbours, std::set<unsigned>());

            //Test another corner node
            std::set<unsigned> node_1_neighbours = node_based_cell_population.GetNeighbouringNodeIndices(1);

            std::set<unsigned> expected_node_1_neighbours;
            expected_node_1_neighbours.insert(2);
            expected_node_1_neighbours.insert(4);

            TS_ASSERT_EQUALS(node_1_neighbours.size(), expected_node_1_neighbours.size());
            TS_ASSERT_EQUALS(node_1_neighbours, expected_node_1_neighbours);

            // Test central node
            std::set<unsigned> node_4_neighbours = node_based_cell_population.GetNeighbouringNodeIndices(4);

            std::set<unsigned> expected_node_4_neighbours;
            expected_node_4_neighbours.insert(1);
            expected_node_4_neighbours.insert(2);
            expected_node_4_neighbours.insert(3);

            TS_ASSERT_EQUALS(node_4_neighbours.size(), expected_node_4_neighbours.size());
            TS_ASSERT_EQUALS(node_4_neighbours, expected_node_4_neighbours);
        }
    }

    void TestGetNodesWithinNeighbourhoodRadius()
    {
        EXIT_IF_PARALLEL;    // Doesn't work in parallel yet until halo nodes are updated (#2364)

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a small node-based cell population
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(node_based_cell_population.GetNumRealCells(), mesh.GetNumNodes());

        // Test GetNodesWithinNeighbourhoodRadius() method
        if (PetscTools::AmMaster())
        {
            // Coverage
            TS_ASSERT_THROWS_THIS(node_based_cell_population.GetNodesWithinNeighbourhoodRadius(0,0.1), "mNodeNeighbours not set up. Call Update() before GetNodesWithinNeighbourhoodRadius()");
        }

        node_based_cell_population.Update();

        if (PetscTools::AmMaster())
        {
            //All nodes should have no neighbours with small search radius
            std::set<unsigned> node_0_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(0,0.1);
            std::set<unsigned> node_1_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(1,0.1);
            std::set<unsigned> node_2_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(2,0.1);
            std::set<unsigned> node_3_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(3,0.1);
            std::set<unsigned> node_4_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(4,0.1);

            TS_ASSERT_EQUALS(node_0_neighbours.size(), 0u);
            TS_ASSERT_EQUALS(node_0_neighbours, std::set<unsigned>());
            TS_ASSERT_EQUALS(node_1_neighbours.size(), 0u);
            TS_ASSERT_EQUALS(node_1_neighbours, std::set<unsigned>());
            TS_ASSERT_EQUALS(node_2_neighbours.size(), 0u);
            TS_ASSERT_EQUALS(node_2_neighbours, std::set<unsigned>());
            TS_ASSERT_EQUALS(node_3_neighbours.size(), 0u);
            TS_ASSERT_EQUALS(node_3_neighbours, std::set<unsigned>());
            TS_ASSERT_EQUALS(node_4_neighbours.size(), 0u);
            TS_ASSERT_EQUALS(node_4_neighbours, std::set<unsigned>());


            //For a slightly  larger search radius, the corner nodes will have exactly 1 neighbours and the centre node will have 4
            node_0_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(0,0.72);
            node_4_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(4,0.72);

            std::set<unsigned> expected_node_0_neighbours;
            expected_node_0_neighbours.insert(4);
            std::set<unsigned> expected_node_4_neighbours;
            expected_node_4_neighbours.insert(0);
            expected_node_4_neighbours.insert(1);
            expected_node_4_neighbours.insert(2);
            expected_node_4_neighbours.insert(3);

            TS_ASSERT_EQUALS(node_0_neighbours.size(), 1u);
            TS_ASSERT_EQUALS(node_0_neighbours, expected_node_0_neighbours);
            TS_ASSERT_EQUALS(node_4_neighbours.size(), 4u);
            TS_ASSERT_EQUALS(node_4_neighbours, expected_node_4_neighbours);

            //For a larger search radius, the corner nodes will have exactly three neighbours and the centre node will have 4
            node_0_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(0,1);
            node_4_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(4,1);

            expected_node_0_neighbours.insert(1);
            expected_node_0_neighbours.insert(3);

            TS_ASSERT_EQUALS(node_0_neighbours.size(), 3u);
            TS_ASSERT_EQUALS(node_0_neighbours, expected_node_0_neighbours);
            TS_ASSERT_EQUALS(node_4_neighbours.size(), 4u);
            TS_ASSERT_EQUALS(node_4_neighbours, expected_node_4_neighbours);

            TS_ASSERT_THROWS_THIS(node_based_cell_population.GetNodesWithinNeighbourhoodRadius(0,2.0), "neighbourhoodRadius should be less than or equal to the  the maximum interaction radius defined on the NodesOnlyMesh");

        }

        //Now test with a bigger box
        mesh.Clear();
        mesh.ConstructNodesWithoutMesh(generating_mesh, 2.0);

        node_based_cell_population.Update();

        if (PetscTools::AmMaster())
        {

            //For a very large search radius, all nodes can see all other nodes so each have 4 neighbours
            std::set<unsigned> node_0_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(0,2);
            std::set<unsigned> node_4_neighbours = node_based_cell_population.GetNodesWithinNeighbourhoodRadius(4,2);

            std::set<unsigned> expected_node_0_neighbours;
            expected_node_0_neighbours.insert(1);
            expected_node_0_neighbours.insert(2);
            expected_node_0_neighbours.insert(3);
            expected_node_0_neighbours.insert(4);
            std::set<unsigned> expected_node_4_neighbours;
            expected_node_4_neighbours.insert(0);
            expected_node_4_neighbours.insert(1);
            expected_node_4_neighbours.insert(2);
            expected_node_4_neighbours.insert(3);

            TS_ASSERT_EQUALS(node_0_neighbours.size(), 4u);
            TS_ASSERT_EQUALS(node_0_neighbours, expected_node_0_neighbours);
            TS_ASSERT_EQUALS(node_4_neighbours.size(), 4u);
            TS_ASSERT_EQUALS(node_4_neighbours, expected_node_4_neighbours);

         }
    }

    void TestSettingCellAncestors()
    {
        // Create a small node-based cell population
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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
        TS_ASSERT_EQUALS(counter, mesh.GetNumNodes());

        // Test that we can recover the remaining number of ancestors
        std::set<unsigned> remaining_ancestors = cell_population.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), mesh.GetNumNodes());

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
        TS_ASSERT_EQUALS(remaining_ancestors.size(), (unsigned)(mesh.GetNumNodes() > 0));
    }

    void TestGetLocationOfCellCentreAndGetWidth()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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
        EXIT_IF_PARALLEL;    // Population writers don't work in parallel yet

        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        node_based_cell_population.Update(); // so cell neighbours are calculated
        for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
             cell_iter != node_based_cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(node_based_cell_population.GetVolumeOfCell(*cell_iter), M_PI*0.5*0.5, 1e-6);
        }

        // For coverage of WriteResultsToFiles()
        node_based_cell_population.GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_apc1(node_based_cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(node_based_cell_population.GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(node_based_cell_population.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(node_based_cell_population.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(node_based_cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(node_based_cell_population.GetCellPropertyRegistry()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_diff_type(node_based_cell_population.GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>());

        node_based_cell_population.GetCellUsingLocationIndex(0)->SetCellProliferativeType(p_transit_type);
        node_based_cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        node_based_cell_population.GetCellUsingLocationIndex(1)->SetCellProliferativeType(p_diff_type);
        node_based_cell_population.GetCellUsingLocationIndex(1)->SetMutationState(p_apc1);
        node_based_cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);
        node_based_cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);
        node_based_cell_population.GetCellUsingLocationIndex(3)->AddCellProperty(p_apoptotic_state);

        node_based_cell_population.rGetMesh().GetNode(0)->SetRadius(0.6); // Default is 0.5
        node_based_cell_population.Update(); // To recalculate cell neighbours

        node_based_cell_population.AddCellWriter<CellIdWriter>();

        // Coverage of writing CellData to VTK
        for (AbstractCellPopulation<2>::Iterator cell_iter = node_based_cell_population.Begin();
             cell_iter != node_based_cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("var0", 0.0);
            cell_iter->GetCellData()->SetItem("var1", 3.0);
        }

        // Test set methods
        std::string output_directory = "TestNodeBasedCellPopulationWriters2d";
        OutputFileHandler output_file_handler(output_directory, false);

        node_based_cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        node_based_cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        node_based_cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        node_based_cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        node_based_cell_population.AddCellWriter<CellAgesWriter>();
        node_based_cell_population.AddCellWriter<CellVolumesWriter>();
        node_based_cell_population.AddCellWriter<CellAncestorWriter>();
        node_based_cell_population.AddCellWriter<CellMutationStatesWriter>();
        node_based_cell_population.AddCellWriter<CellAppliedForceWriter>();

        // Set some forces for the applied force writer
        c_vector<double, 2> force_0 = Create_c_vector(1.2, 2.3);
        c_vector<double, 2> force_1 = Create_c_vector(2.3, 3.4);
        c_vector<double, 2> force_2 = Create_c_vector(3.4, 4.5);
        c_vector<double, 2> force_3 = Create_c_vector(4.5, 5.6);
        node_based_cell_population.rGetMesh().GetNode(0)->rGetAppliedForce() = force_0;
        node_based_cell_population.rGetMesh().GetNode(1)->rGetAppliedForce() = force_1;
        node_based_cell_population.rGetMesh().GetNode(2)->rGetAppliedForce() = force_2;
        node_based_cell_population.rGetMesh().GetNode(3)->rGetAppliedForce() = force_3;

        node_based_cell_population.SetCellAncestorsToLocationIndices();

        node_based_cell_population.OpenWritersFiles(output_file_handler);
        node_based_cell_population.WriteResultsToFiles(output_directory);
        node_based_cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.vizcelltypes").CompareFiles();
        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.vizancestors").CompareFiles();
        FileComparison(results_dir + "results.vizmutationstates", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.vizmutationstates").CompareFiles();
        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/cellmutationstates.dat").CompareFiles();
        FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/cellages.dat").CompareFiles();
        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/cellareas.dat").CompareFiles();
        FileComparison(results_dir + "cellappliedforce.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/cellappliedforce.dat").CompareFiles();

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = node_based_cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 1u);

        // Test the GetCellProliferativeTypeCount() function
        std::vector<unsigned> cell_types = node_based_cell_population.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 4u);
        TS_ASSERT_EQUALS(cell_types[0], 2u);
        TS_ASSERT_EQUALS(cell_types[1], 1u);
        TS_ASSERT_EQUALS(cell_types[2], 1u);
        TS_ASSERT_EQUALS(cell_types[3], 0u);

        // Test the Get MechanicsCutOfLengthMethods
        TS_ASSERT_DELTA(node_based_cell_population.GetMechanicsCutOffLength(),1.5, 1e-9);

        // For coverage
        TS_ASSERT_THROWS_NOTHING(node_based_cell_population.WriteResultsToFiles(output_directory));

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        node_based_cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare output with saved files of what they should look like
        FileComparison(results_dir + "results.parameters", "cell_based/test/data/TestNodeBasedCellPopulationWriters2d/results.parameters").CompareFiles();

        // Test VTK output
#ifdef CHASTE_VTK
        node_based_cell_population.WriteVtkResultsToFile(output_directory);

        // Read VTK file and check it doesn't cause any problems
        VtkMeshReader<2,2> vtk_reader(results_dir + "/results_0.vtu");

        // Test that the correct mutation states were recorded
        std::vector<double> mutation_states_data;
        std::vector<double> saved_mutation_states_data;
        saved_mutation_states_data.push_back(5.0);
        saved_mutation_states_data.push_back(3.0);
        saved_mutation_states_data.push_back(4.0);
        saved_mutation_states_data.push_back(4.0);

        vtk_reader.GetPointData("Mutation states", mutation_states_data);
        TS_ASSERT_EQUALS(mutation_states_data.size(), 4u);
        for (unsigned i=0; i<mutation_states_data.size(); i++)
        {
            TS_ASSERT_DELTA(mutation_states_data[i], saved_mutation_states_data[i], 1e-9);
        }

        // Test that the correct cell proliferative types were recorded
        std::vector<double> proliferative_types_data;
        vtk_reader.GetPointData("Cell types", proliferative_types_data);
        TS_ASSERT_EQUALS(proliferative_types_data.size(), 4u);
        TS_ASSERT_DELTA(proliferative_types_data[0], 5.0, 1e-9);
        TS_ASSERT_DELTA(proliferative_types_data[1], 3.0, 1e-9);
        TS_ASSERT_DELTA(proliferative_types_data[2], 4.0, 1e-9);
        TS_ASSERT_DELTA(proliferative_types_data[3], 6.0, 1e-9);

        // Test that the correct cell volumes were recorded
        std::vector<double> cell_volumes_data;
        vtk_reader.GetPointData("Cell volumes", cell_volumes_data);
        TS_ASSERT_EQUALS(cell_volumes_data.size(), 4u);
        TS_ASSERT_DELTA(cell_volumes_data[0], 1.0690, 1e-3);
        TS_ASSERT_DELTA(cell_volumes_data[1], 0.7594, 1e-3);
        TS_ASSERT_DELTA(cell_volumes_data[2], 0.7853, 1e-3);
        TS_ASSERT_DELTA(cell_volumes_data[3], 0.7594, 1e-3);

        // Test that the correct cell cycle phases were recorded
        std::vector<double> cycle_phases_data;
        vtk_reader.GetPointData("Cycle phases", cycle_phases_data);
        TS_ASSERT_EQUALS(cycle_phases_data.size(), 4u);
        for (unsigned i=0; i<cycle_phases_data.size(); i++)
        {
            TS_ASSERT_DELTA(cycle_phases_data[i], 4.0, 1e-9);
        }

        // Test that the correct cell ages were recorded
        std::vector<double> ages_data;
        vtk_reader.GetPointData("Ages", ages_data);
        TS_ASSERT_EQUALS(ages_data.size(), 4u);
        TS_ASSERT_DELTA(ages_data[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(ages_data[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(ages_data[2], 2.0, 1e-9);
        TS_ASSERT_DELTA(ages_data[3], 3.0, 1e-9);

        // Test that the correct ancestors were recorded
        std::vector<double> ancestors_data;
        vtk_reader.GetPointData("Ancestors", ancestors_data);
        TS_ASSERT_EQUALS(ancestors_data.size(), 4u);
        TS_ASSERT_DELTA(ancestors_data[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(ancestors_data[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(ancestors_data[2], 2.0, 1e-9);
        TS_ASSERT_DELTA(ancestors_data[3], 3.0, 1e-9);
#endif
    }

    void TestNodeBasedCellPopulationOutputWriters3d()
    {
        EXIT_IF_PARALLEL;    // Population writers don't work in parallel yet

        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.Update(); // so cell neighbours are calculated when outputting volume

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "NodeBasedCellPopulation-3");

        // Test set methods
        std::string output_directory = "TestNodeBasedCellPopulationWriters3d";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellAppliedForceWriter>();

        // Set some forces for the applied force writer
        c_vector<double, 3> force_0 = Create_c_vector(1.2, 2.3, 3.4);
        c_vector<double, 3> force_1 = Create_c_vector(2.3, 3.4, 4.5);
        c_vector<double, 3> force_2 = Create_c_vector(3.4, 4.5, 5.6);
        c_vector<double, 3> force_3 = Create_c_vector(4.5, 5.6, 6.7);
        cell_population.rGetMesh().GetNode(0)->rGetAppliedForce() = force_0;
        cell_population.rGetMesh().GetNode(1)->rGetAppliedForce() = force_1;
        cell_population.rGetMesh().GetNode(2)->rGetAppliedForce() = force_2;
        cell_population.rGetMesh().GetNode(3)->rGetAppliedForce() = force_3;

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.AddCellWriter<CellAncestorWriter>();

        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteResultsToFiles(output_directory);
        cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/results.vizcelltypes").CompareFiles();
        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/results.vizancestors").CompareFiles();
        FileComparison(results_dir + "results.vizmutationstates", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/results.vizmutationstates").CompareFiles();
        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/cellmutationstates.dat").CompareFiles();
        FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/cellages.dat").CompareFiles();
        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/cellareas.dat").CompareFiles();
        FileComparison(results_dir + "cellappliedforce.dat", "cell_based/test/data/TestNodeBasedCellPopulationWriters3d/cellappliedforce.dat").CompareFiles();

        // Test VTK output
#ifdef CHASTE_VTK
        cell_population.WriteVtkResultsToFile(output_directory);

        // Read VTK file and check it doesn't cause any problems
        VtkMeshReader<3,3> vtk_reader(results_dir + "/results_0.vtu");

        // All cells are wild type
        std::vector<double> mutation_states_data;
        vtk_reader.GetPointData("Mutation states", mutation_states_data);
        TS_ASSERT_EQUALS(mutation_states_data.size(), 51u);
        for (unsigned i=0; i<mutation_states_data.size(); i++)
        {
            TS_ASSERT_DELTA(mutation_states_data[i], 0.0, 1e-9);
        }

        // Test that the correct cell cycle phases were recorded
        std::vector<double> cycle_phases_data;
        vtk_reader.GetPointData("Cycle phases", cycle_phases_data);
        TS_ASSERT_EQUALS(cycle_phases_data.size(), 51u);
        for (unsigned i=0; i<cycle_phases_data.size(); i++)
        {
            TS_ASSERT_DELTA(cycle_phases_data[i], 4.0, 1e-9);
        }

        // Test that the correct cell ages were recorded
        std::vector<double> ages_data;
        vtk_reader.GetPointData("Ages", ages_data);
        TS_ASSERT_EQUALS(ages_data.size(), 51u);
        for (unsigned i=0; i<cycle_phases_data.size(); i++)
        {
            TS_ASSERT_DELTA(ages_data[i], i, 1e-9);
        }

        // Test that the correct ancestors were recorded
        std::vector<double> ancestors_data;
        vtk_reader.GetPointData("Ancestors", ancestors_data);
        TS_ASSERT_EQUALS(ancestors_data.size(), 51u);
        for (unsigned i=0; i<cycle_phases_data.size(); i++)
        {
            TS_ASSERT_DELTA(ancestors_data[i], i, 1e-9);
        }
#endif
    }

    void TestWritingCellCyclePhases()
    {
        EXIT_IF_PARALLEL;    // Population writers dont work in parallel yet.

        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        cells[0]->SetBirthTime(-23.5);
        cells[1]->SetBirthTime(-0.5);
        cells[2]->SetBirthTime(-1.5);
        cells[3]->SetBirthTime(-15.5);
        cells[4]->SetBirthTime(-23.5);
        cells[0]->SetCellProliferativeType(p_diff_type);

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

        node_based_cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        node_based_cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        node_based_cell_population.OpenWritersFiles(output_file_handler);
        node_based_cell_population.WriteResultsToFiles(output_directory);
        node_based_cell_population.CloseWritersFiles();

        // Test the GetCellCyclePhaseCount() function
        std::vector<unsigned> cell_cycle_phases = node_based_cell_population.GetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phases[0], 1u);
        TS_ASSERT_EQUALS(cell_cycle_phases[1], 3u);
        TS_ASSERT_EQUALS(cell_cycle_phases[2], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[3], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[4], 1u);


        // Check throws error when using a non phased based CCM
        std::vector<CellPtr> cells_2;
        CellsGenerator<BernoulliTrialCellCycleModel, 2> cells_generator_2;
        cells_generator_2.GenerateBasic(cells_2, mesh.GetNumNodes());

        // Create a cell population
        NodeBasedCellPopulation<2> node_based_cell_population_2(mesh, cells_2);

        TS_ASSERT_THROWS_THIS(node_based_cell_population_2.GetCellCyclePhaseCount(),"You are trying to record the cell cycle phase of cells with a non phase based cell cycle model.");


    }

    void TestArchivingCellPopulation()
    {
        EXIT_IF_PARALLEL;    // Population archiving doesn't work in parallel yet.

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
            mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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

            p_cell_population->SetUseVariableRadii(true);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Avoid memory leak
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
            TS_ASSERT(p_cell_population->GetUseVariableRadii());

            // Tidy up
            delete p_cell_population;
        }
    }

    void TestAddAndRemoveMovedCell()
    {
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 0.0, 1.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        if (PetscTools::AmMaster())
        {
            unsigned num_initial_cells = cell_population.GetNumNodes();

            boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);

            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));

            boost::shared_ptr<Node<2> > p_node(new Node<2>(PetscTools::GetNumProcs() * (num_initial_cells + 1), false, 0.0, 0.5));
            cell_population.AddMovedCell(p_cell, p_node);

            TS_ASSERT_EQUALS(cell_population.GetNumNodes(), num_initial_cells +1);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_initial_cells + 1);
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell), PetscTools::GetNumProcs() * (num_initial_cells + 1));

            // Now delete the first cell
            cell_population.DeleteMovedCell(0);

            TS_ASSERT_EQUALS(cell_population.GetNumNodes(), num_initial_cells);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_initial_cells);
            //TS_ASSERT_THROWS_THIS(cell_population.GetLocationIndexUsingCell(cells[0]), "");
        }

        // Clean up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestGetTetrahedralMeshForPdeModifier()
    {
        EXIT_IF_PARALLEL;  // The population.GetTetrahedralMeshForPdeModifier() method does not yet work in parallel.

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.2);
        mesh.SetCalculateNodeNeighbours(false);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        TetrahedralMesh<2,2>* p_tet_mesh = cell_population.GetTetrahedralMeshForPdeModifier();

        // Check it has the correct number of nodes and elements
        TS_ASSERT_EQUALS(p_tet_mesh->GetNumNodes(), cell_population.GetNumNodes());
        TS_ASSERT_EQUALS(p_tet_mesh->GetNumElements(), 4u);

        // Check some nodes have the correct locations
        TS_ASSERT_DELTA(p_tet_mesh->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_tet_mesh->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_tet_mesh->GetNode(1)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(p_tet_mesh->GetNode(1)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_tet_mesh->GetNode(2)->rGetLocation()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(p_tet_mesh->GetNode(2)->rGetLocation()[1], 1.0, 1e-6);

        // Tidy up
        delete p_tet_mesh;
    }

    void TestGetCellDataItemAtPdeNode()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        std::string var_name = "foo";
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_THROWS_THIS(cell_population.GetCellDataItemAtPdeNode(0,var_name),
                    "The item foo is not stored");

        cell_population.GetCellUsingLocationIndex(0)->GetCellData()->SetItem(var_name, 3.14);

        TS_ASSERT_DELTA(cell_population.GetCellDataItemAtPdeNode(0,var_name), 3.14, 1e-6);
        }
        else
        {
            TS_ASSERT_THROWS_THIS(cell_population.GetCellDataItemAtPdeNode(0,var_name),
                                "Location index input argument does not correspond to a Cell");
        }
        // Coverage of GetOutputResultsForChasteVisualizer()
        TS_ASSERT_EQUALS(cell_population.GetOutputResultsForChasteVisualizer(), true);
    }
};

#endif /*TESTNODEBASEDCELLPOPULATION_HPP_*/
