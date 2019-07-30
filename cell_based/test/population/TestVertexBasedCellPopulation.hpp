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

#ifndef TESTVERTEXBASEDCELLPOPULATION_HPP_
#define TESTVERTEXBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "FixedVertexBasedDivisionRule.hpp"
#include "ApoptoticCellProperty.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"

// Test warnings
#include "Warnings.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestVertexBasedCellPopulation : public AbstractCellBasedTestSuite
{
public:

    // Test construction, accessors and iterator
    void TestCreateSmallVertexBasedCellPopulationAndGetWidth()
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        unsigned num_cells = cells.size();
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), num_cells);

        unsigned counter = 0;

        // Test VertexBasedCellPopulation::Iterator
        for (VertexBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            TS_ASSERT_EQUALS(counter, p_mesh->GetElement(counter)->GetIndex());

            counter++;
        }

        // Test we have gone through all cells in the for loop
        TS_ASSERT_EQUALS(counter, cell_population.GetNumRealCells());

        // Test GetNumNodes() method
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), p_mesh->GetNumNodes());

        // Test GetWidth() method
        double width_x = cell_population.GetWidth(0);
        TS_ASSERT_DELTA(width_x, 5.5000, 1e-4);

        double width_y = cell_population.GetWidth(1);
        TS_ASSERT_DELTA(width_y, 2.8867, 1e-4);

        // For coverage of GetVolumeOfCell()
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), sqrt(3.0)/2, 1e-6);
        }

        // For coverage of GetRosetteRankOfCell()
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetRosetteRankOfCell(*cell_iter), 3u);
        }

        // Test GetNeighbouringLocationIndices() method
        std::set<unsigned> expected_neighbours_of_cell_0;
        expected_neighbours_of_cell_0.insert(1);
        expected_neighbours_of_cell_0.insert(5);

        std::set<unsigned> neighbours_of_cell_0 = cell_population.GetNeighbouringLocationIndices(*(cell_population.Begin()));
        TS_ASSERT(neighbours_of_cell_0 == expected_neighbours_of_cell_0);

        // For coverage, test that GetDefaultTimeStep() returns the correct value
        TS_ASSERT_DELTA(cell_population.GetDefaultTimeStep(), 0.002, 1e-6);
    }

    void TestValidate()
    {
        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements()-1);

        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<cells.size(); i++)
        {
            cell_location_indices.push_back(i);
        }

        // This should throw an exception as the number of cells does not equal the number of elements
        std::vector<CellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(VertexBasedCellPopulation<2> cell_population(*p_mesh, cells_copy),
                "At time 0, Element 8 does not appear to have a cell associated with it");

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell1(new Cell(p_state, p_model));
        p_cell1->SetCellProliferativeType(p_stem_type);

        double birth_time = 0.0 - p_mesh->GetNumElements()-1;
        p_cell1->SetBirthTime(birth_time);

        cells.push_back(p_cell1);
        cell_location_indices.push_back(p_mesh->GetNumElements()-1);

        // This should pass as the number of cells equals the number of elements
        std::vector<CellPtr> cells_copy2(cells);
        TS_ASSERT_THROWS_NOTHING(VertexBasedCellPopulation<2> cell_population(*p_mesh, cells_copy2));

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Check correspondence between elements and cells
        for (VertexMesh<2,2>::VertexElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            std::set<unsigned> expected_node_indices;
            unsigned expected_index = iter->GetIndex();

            for (unsigned i=0; i<iter->GetNumNodes(); i++)
            {
                expected_node_indices.insert(iter->GetNodeGlobalIndex(i));
            }

            std::set<unsigned> actual_node_indices;
            unsigned elem_index = iter->GetIndex();
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(elem_index);
            VertexElement<2,2>* p_actual_element = cell_population.GetElementCorrespondingToCell(p_cell);
            unsigned actual_index = p_actual_element->GetIndex();

            for (unsigned i=0; i<p_actual_element->GetNumNodes(); i++)
            {
                actual_node_indices.insert(p_actual_element->GetNodeGlobalIndex(i));
            }

            TS_ASSERT_EQUALS(actual_index, expected_index);
            TS_ASSERT_EQUALS(actual_node_indices, expected_node_indices);
        }

        // Create anoter simple vertex-based mesh
        HoneycombVertexMeshGenerator generator2(3, 3);
        MutableVertexMesh<2,2>* p_mesh2 = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells2;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, p_mesh2->GetNumElements()+1);

        std::vector<unsigned> cell_location_indices2;
        for (unsigned i=0; i<cells2.size(); i++)
        {
            cell_location_indices2.push_back(i%p_mesh2->GetNumElements()); // Element 0 will have 2 cells
        }

        // This should throw an exception as the number of cells
        // does not equal the number of elements
        TS_ASSERT_THROWS_THIS(VertexBasedCellPopulation<2> cell_population2(*p_mesh2, cells2, false, true, cell_location_indices2),
                "At time 0, Element 0 appears to have 2 cells associated with it");
    }

    void TestGetDampingConstant()
    {
        // Create mesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApcOneHitCellMutationState, p_apc1);
        MAKE_PTR(ApcTwoHitCellMutationState, p_apc2);
        MAKE_PTR(BetaCateninOneHitCellMutationState, p_bcat1);
        MAKE_PTR(CellLabel, p_label);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        cells[0]->SetMutationState(p_apc1);
        cells[6]->SetMutationState(p_apc2);
        cells[7]->SetMutationState(p_bcat1);
        cells[8]->AddCellProperty(p_label);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.InitialiseCells(); // this method must be called explicitly as there is no simulation

        cell_population.SetDampingConstantMutant(8.0);

        // Test GetDampingConstant()
        double normal_damping_constant = cell_population.GetDampingConstantNormal();
        double mutant_damping_constant = cell_population.GetDampingConstantMutant();

        // Node 0 is contained in cell 0 only, therefore should have mutant damping constant
        double damping_constant_at_node_0 = cell_population.GetDampingConstant(0);
        TS_ASSERT_DELTA(damping_constant_at_node_0, mutant_damping_constant, 1e-6);

        // Node 1 is contained in cell 2 only, therefore should have a normal damping constant
        double damping_constant_at_node_1 = cell_population.GetDampingConstant(1);
        TS_ASSERT_DELTA(damping_constant_at_node_1, normal_damping_constant, 1e-6);

        // Node 4 is contained in cells 0 and 1, therefore should an averaged damping constant
        double damping_constant_at_node_4 = cell_population.GetDampingConstant(4);
        TS_ASSERT_DELTA(damping_constant_at_node_4, (normal_damping_constant+mutant_damping_constant)/2.0, 1e-6);

        // Node 8 is contained in cells 0, 1, 3, therefore should an averaged damping constant
        double damping_constant_at_node_8 = cell_population.GetDampingConstant(8);
        TS_ASSERT_DELTA(damping_constant_at_node_8, (2*normal_damping_constant+mutant_damping_constant)/3.0, 1e-6);

        // Node 27 is contained in cell 6 only, therefore should have a mutant damping constant
        double damping_constant_at_node_27 = cell_population.GetDampingConstant(27);
        TS_ASSERT_DELTA(damping_constant_at_node_27, mutant_damping_constant, 1e-6);

        // Node 25 is contained in cells 6 and 7, therefore should have a mutant damping constant
        double damping_constant_at_node_25 = cell_population.GetDampingConstant(25);
        TS_ASSERT_DELTA(damping_constant_at_node_25, (normal_damping_constant+mutant_damping_constant)/2.0, 1e-6);
    }

    /**
     * This test covers the case where GetDampingConstant() is called on a node
     * that is not contained in any elements.
     *
     * \todo it may be sensible to check that each node is contained in at least
     * one element in a method such as Validate()
     */
    void TestGetDampingConstantException()
    {
        // Create a simple vertex mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(2, true, 0.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(4, true, 0.0, 2.0));

        // Make a square element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem;
        nodes_elem.push_back(nodes[0]);
        nodes_elem.push_back(nodes[1]);
        nodes_elem.push_back(nodes[2]);
        nodes_elem.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Create cells
        MAKE_PTR(ApcOneHitCellMutationState, p_apc1);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements());
        cells[0]->SetMutationState(p_apc1);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);
        cell_population.InitialiseCells(); // this method must be called explicitly as there is no simulation
        cell_population.SetDampingConstantMutant(8.0);

        // Node 0 is contained in element 0 only, therefore should have mutant damping constant
        TS_ASSERT_DELTA(cell_population.GetDampingConstant(0), 8.0, 1e-6);

        // Node 4 is not contained in any elements, so GetDampingConstant() should throw an exception
        TS_ASSERT_THROWS_THIS(cell_population.GetDampingConstant(4),
            "At time 0, Node 4 is not contained in any elements, so GetDampingConstant() returns zero");
    }

    void TestUpdateWithoutBirthOrDeath()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        unsigned num_cells_removed = cell_population.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_cells_removed, 0u);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_THROWS_NOTHING(cell_population.Update());
    }

    void TestAddNode()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 68u);

        Node<2>* p_node = new Node<2>(0, true, -15.0, 23.0);
        cell_population.AddNode(p_node);

        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 69u);
    }

    void TestAddCellWithSimpleMesh()
    {
        // Make some nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, true, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, true, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, true, 0.0, 2.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        // Make a triangular element out of nodes 1,4,2
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[2]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // For coverage, test GetLocationOfCellCentre()

        // Cell 0 is a rectangle with centre of mass (0,0)
        c_vector<double, 2> cell0_location = cell_population.GetLocationOfCellCentre(cell_population.GetCellUsingLocationIndex(0));
        TS_ASSERT_DELTA(cell0_location[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell0_location[1], 0.0, 1e-4);

        // Cell 1 is a triangle with centre of mass (0,4/3)
        c_vector<double, 2> cell1_location = cell_population.GetLocationOfCellCentre(cell_population.GetCellUsingLocationIndex(1));
        TS_ASSERT_DELTA(cell1_location[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell1_location[1], 4.0/3.0, 1e-4);

        unsigned old_num_nodes = vertex_mesh.GetNumNodes();
        unsigned old_num_elements = vertex_mesh.GetNumElements();
        unsigned old_num_cells = cell_population.rGetCells().size();

        // Add new cell by dividing element 0 along short axis (this is the standard division rule)

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_temp_cell(new Cell(p_state, p_model));
        p_temp_cell->SetCellProliferativeType(p_stem_type);
        p_temp_cell->SetBirthTime(-1);

        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule = cell_population.GetVertexBasedDivisionRule();
        c_vector<double, 2> short_axis = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        TS_ASSERT_DELTA(short_axis[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(short_axis[1], 1.0, 1e-9);

        CellPtr p_new_cell = cell_population.AddCell(p_temp_cell, p_cell0);

        // Check that the new cell was successfully added to the cell population
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), old_num_nodes+2);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), old_num_elements+1);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), old_num_cells+1);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), old_num_elements+1);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[1], 1.0, 1e-12);

        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[1], -1.0, 1e-12);

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(cell_population.GetElement(1)->GetNodeGlobalIndex(3), 5u);

        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(cell_population.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);
        expected_elements_containing_node_5.insert(2);

        TS_ASSERT_EQUALS(cell_population.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(2);

        TS_ASSERT_EQUALS(cell_population.GetNode(6)->rGetContainingElementIndices(), expected_elements_containing_node_6);

        // Check the index of the new cell
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_new_cell), old_num_elements);
    }

    void TestAddCellWithGivenDivisionVector()
    {
        // Make a vertex mesh consisting of a single square element
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 2.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));

        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem.push_back(nodes[i]);
        }
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Create a cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(-20.0);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);

        // Add a new cell by dividing element 0 along the axis (1,0)
        c_vector<double,2> cell_division_axis;
        cell_division_axis[0] = 1.0;
        cell_division_axis[1] = 0.0;

        MAKE_PTR_ARGS(FixedVertexBasedDivisionRule<2>, p_div_rule, (cell_division_axis));
        cell_population.SetVertexBasedDivisionRule(p_div_rule);

        FixedG1GenerationalCellCycleModel* p_model2 = new FixedG1GenerationalCellCycleModel();
        CellPtr p_temp_cell(new Cell(p_state, p_model2));
        p_temp_cell->SetCellProliferativeType(p_stem_type);
        p_temp_cell->SetBirthTime(-1.0);

        CellPtr p_new_cell = cell_population.AddCell(p_temp_cell, p_cell);

        // Check that the new cell was successfully added to the cell population
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(cell_population.GetNode(4)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(4)->rGetLocation()[1], 0.5, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetLocation()[1], 0.5, 1e-12);

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_4;
        expected_elements_containing_node_4.insert(0);
        expected_elements_containing_node_4.insert(1);
        TS_ASSERT_EQUALS(cell_population.GetNode(4)->rGetContainingElementIndices(), expected_elements_containing_node_4);

        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);
        TS_ASSERT_EQUALS(cell_population.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        // Check the index of the new cell
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_new_cell), 1u);
    }

    void TestAddCellWithHoneycombMesh()
    {
        // Create a mesh with 9 elements
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));

            double birth_time = 0.0 - elem_index;

            // Cell 4 should divide immediately
            if (elem_index==4)
            {
                p_cell->SetCellProliferativeType(p_transit_type); // as stem cells always divide horizontally
                birth_time = -50.0;
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff_type);
            }

            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Initialise cells (usually called in cell-based simulation constructor)
        cell_population.InitialiseCells();

        unsigned old_num_nodes = p_mesh->GetNumNodes();
        unsigned old_num_elements = p_mesh->GetNumElements();
        unsigned old_num_cells = cell_population.rGetCells().size();

        // Test GetNeighbouringNodeIndices() method
        std::set<unsigned> node_10_neighbours = cell_population.GetNeighbouringNodeIndices(10);

        std::set<unsigned> expected_node_10_neighbours;
        expected_node_10_neighbours.insert(6);
        expected_node_10_neighbours.insert(13);
        expected_node_10_neighbours.insert(14);

        TS_ASSERT_EQUALS(node_10_neighbours.size(), expected_node_10_neighbours.size());
        TS_ASSERT_EQUALS(node_10_neighbours, expected_node_10_neighbours);

        // Add a new cell by dividing cell 4

        cell_population.GetCellUsingLocationIndex(4)->ReadyToDivide();
        CellPtr p_cell4 = cell_population.GetCellUsingLocationIndex(4);

        CellPtr p_new_cell = cell_population.GetCellUsingLocationIndex(4)->Divide();

        // Add new cell to the cell population
        cell_population.AddCell(p_new_cell, p_cell4);

        // Check that the new cell was successfully added to the cell population
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), old_num_nodes+2);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), old_num_elements+1);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), old_num_elements+1);

        // Check the location of the new nodes
        // (sensitive to changes in random number generation)
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[0], 2.2748, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[1], 1.8620, 1e-4);

        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[0], 1.7251, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[1], 1.0247, 1e-4);

        // Now test the nodes in each element
        for (unsigned i=0; i<cell_population.GetNumElements(); i++)
        {
            if (i==4 || i==9)
            {
                // Elements 4 and 9 should each have one less node
                TS_ASSERT_EQUALS(cell_population.GetElement(i)->GetNumNodes(), 5u);
            }
            else if (i==1 || i==8)
            {
                // Elements 3 and 8 should each have one extra node
                TS_ASSERT_EQUALS(cell_population.GetElement(i)->GetNumNodes(), 7u);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_population.GetElement(i)->GetNumNodes(), 6u);
            }
        }

        // Check node ownership for a few elements

        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(2), 8u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(3), 11u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(4), 7u);
        TS_ASSERT_EQUALS(cell_population.GetElement(0)->GetNodeGlobalIndex(5), 3u);

        TS_ASSERT_EQUALS(cell_population.GetElement(4)->GetNodeGlobalIndex(0), 9u);
        TS_ASSERT_EQUALS(cell_population.GetElement(4)->GetNodeGlobalIndex(1), 13u);
        TS_ASSERT_EQUALS(cell_population.GetElement(4)->GetNodeGlobalIndex(2), 17u);
        TS_ASSERT_EQUALS(cell_population.GetElement(4)->GetNodeGlobalIndex(3), 30u);
        TS_ASSERT_EQUALS(cell_population.GetElement(4)->GetNodeGlobalIndex(4), 31u);

        TS_ASSERT_EQUALS(cell_population.GetElement(8)->GetNodeGlobalIndex(0), 17u);
        TS_ASSERT_EQUALS(cell_population.GetElement(8)->GetNodeGlobalIndex(1), 22u);
        TS_ASSERT_EQUALS(cell_population.GetElement(8)->GetNodeGlobalIndex(2), 26u);
        TS_ASSERT_EQUALS(cell_population.GetElement(8)->GetNodeGlobalIndex(3), 29u);
        TS_ASSERT_EQUALS(cell_population.GetElement(8)->GetNodeGlobalIndex(4), 25u);
        TS_ASSERT_EQUALS(cell_population.GetElement(8)->GetNodeGlobalIndex(5), 21u);
        TS_ASSERT_EQUALS(cell_population.GetElement(8)->GetNodeGlobalIndex(6), 30u);

        // Test element ownership for a few nodes

        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(1);
        expected_elements_containing_node_5.insert(2);

        TS_ASSERT_EQUALS(cell_population.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_13;
        expected_elements_containing_node_13.insert(2);
        expected_elements_containing_node_13.insert(4);
        expected_elements_containing_node_13.insert(5);

        TS_ASSERT_EQUALS(cell_population.GetNode(13)->rGetContainingElementIndices(), expected_elements_containing_node_13);

        std::set<unsigned> expected_elements_containing_node_30;
        expected_elements_containing_node_30.insert(8);
        expected_elements_containing_node_30.insert(4);
        expected_elements_containing_node_30.insert(9);

        TS_ASSERT_EQUALS(cell_population.GetNode(30)->rGetContainingElementIndices(), expected_elements_containing_node_30);
    }

    void TestIsCellAssociatedWithADeletedLocation()
    {
        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->GetElement(5)->MarkAsDeleted();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population but do not try to validate
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells, false, false);

        // Test IsCellAssociatedWithADeletedLocation() method
        for (VertexBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            bool is_deleted = cell_population.IsCellAssociatedWithADeletedLocation(*cell_iter);
            bool cell_has_index_5 = (cell_population.GetLocationIndexUsingCell(*cell_iter) == 5);

            TS_ASSERT_EQUALS(is_deleted, cell_has_index_5);
        }
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        cells[5]->StartApoptosis();

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 24u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 24u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 24u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 68u);

        p_simulation_time->IncrementTimeOneStep();

        // This is where the two counters should differ
        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 24u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 23u);

        // Remove dead cells
        unsigned num_cells_removed = cell_population.RemoveDeadCells();

        TS_ASSERT_EQUALS(num_cells_removed, 1u);

        // We should now have one less real cell, since one cell has been
        // marked as dead, so is skipped by the cell population iterator
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 68u);

        cell_population.Update();

        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumAllCells(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 68u);

        // Finally, check the cells' element indices have updated

        // We expect the cell element indices to be {0,11,...,23}
        std::set<unsigned> expected_elem_indices;
        for (unsigned i=0; i<cell_population.GetNumRealCells(); i++)
        {
            expected_elem_indices.insert(i);
        }

        // Get actual cell element indices
        std::set<unsigned> element_indices;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Record element index corresponding to cell
            unsigned element_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            element_indices.insert(element_index);
        }

        TS_ASSERT_EQUALS(element_indices, expected_elem_indices);

        // Test Warnings
        // The case where the warning doesn't occur is tested in TestT2SwapCellKiller
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), "A Cell is removed without performing a T2 swap. This could leave a void in the mesh.");
        Warnings::QuietDestroy();
    }

    void TestVertexBasedCellPopulationWriteResultsToFile()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create a simple vertex-based cell population, comprising various cell types in various cell cycle phases
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        boost::shared_ptr<AbstractCellProperty> p_stem(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_transit(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_diff(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_wildtype(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        std::vector<CellPtr> cells;
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_cell(new Cell(p_wildtype, p_model));
            if (elem_index%3 == 0)
            {
                p_cell->SetCellProliferativeType(p_stem);
            }
            else if (elem_index%3 == 1)
            {
                p_cell->SetCellProliferativeType(p_transit);
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff);
            }

            double birth_time = 0.0 - elem_index;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.InitialiseCells();
        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "VertexBasedCellPopulation-2");

        // Allocate some cells to have a different cell mutation state, cell label or apoptotic cell property
        cell_population.GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_apc1(cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(cell_population.GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(cell_population.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(cell_population.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(5)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(1)->SetMutationState(p_apc1);
        cell_population.GetCellUsingLocationIndex(7)->SetMutationState(p_apc1);
        cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);
        cell_population.GetCellUsingLocationIndex(12)->SetMutationState(p_apc2);
        cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);
        cell_population.GetCellUsingLocationIndex(21)->SetMutationState(p_bcat1);
        cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_apoptotic_state);
        cell_population.GetCellUsingLocationIndex(23)->AddCellProperty(p_apoptotic_state);
        cell_population.SetCellAncestorsToLocationIndices();

        cell_population.AddPopulationWriter<VertexT1SwapLocationsWriter>();
        cell_population.AddPopulationWriter<VertexT3SwapLocationsWriter>();
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellLocationIndexWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Coverage of writing CellData to VTK
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("var0", 0.0);
            cell_iter->GetCellData()->SetItem("var1", 3.0);
        }

        std::string output_directory = "TestVertexBasedCellPopulationWriteResultsToFile";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteResultsToFiles(output_directory);

        SimulationTime::Instance()->IncrementTimeOneStep();
        cell_population.Update();
        cell_population.WriteResultsToFiles(output_directory);

        cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes",       "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.vizelements").CompareFiles();

        FileComparison(results_dir + "cellages.dat",              "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/cellages.dat").CompareFiles();
        FileComparison(results_dir + "results.vizancestors",      "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.vizancestors").CompareFiles();
        FileComparison(results_dir + "loggedcell.dat",            "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/loggedcell.dat").CompareFiles();
        FileComparison(results_dir + "results.vizlabels",         "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.vizlabels").CompareFiles();
        FileComparison(results_dir + "results.vizlocationindices",      "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.vizlocationindices").CompareFiles();
        FileComparison(results_dir + "results.vizmutationstates", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.vizmutationstates").CompareFiles();
        FileComparison(results_dir + "results.vizcellphases", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.vizcellphases").CompareFiles();
        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.vizcelltypes").CompareFiles();
        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/cellareas.dat").CompareFiles();

        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/cellmutationstates.dat").CompareFiles();
        FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/cellcyclephases.dat").CompareFiles();
        FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/celltypes.dat").CompareFiles();

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 18u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 2u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 2u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 2u);

        // Test the GetCellProliferativeTypeCount() function
        std::vector<unsigned> cell_types = cell_population.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 4u);
        TS_ASSERT_EQUALS(cell_types[0], 8u);
        TS_ASSERT_EQUALS(cell_types[1], 8u);
        TS_ASSERT_EQUALS(cell_types[2], 8u);
        TS_ASSERT_EQUALS(cell_types[3], 0u);

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        FileComparison( results_dir + "results.parameters", "cell_based/test/data/TestVertexBasedCellPopulationWriteResultsToFile/results.parameters").CompareFiles();

#ifdef CHASTE_VTK
        // Test that VTK writer has produced some files

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_1.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());

        // PVD file
        FileFinder vtk_file3(results_dir + "results.pvd", RelativeTo::Absolute);
        TS_ASSERT(vtk_file3.Exists());
 #endif //CHASTE_VTK
    }

    void TestArchiving2dVertexBasedCellPopulation()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vertex_cell_population_2d.arch";
        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("vertex_mesh_2d");

        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Archive cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumElements());

            // Create cell population
            VertexBasedCellPopulation<2>* const p_cell_population = new VertexBasedCellPopulation<2>(mesh, cells);

            // Cells have been given birth times of 0 and -1.
            // Loop over them to run to time 0.0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Coverage
            p_cell_population->SetOutputCellRearrangementLocations(false);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&> (*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        // Restore cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore the cell population
            (*p_arch) >> *p_simulation_time;
            VertexBasedCellPopulation<2>* p_cell_population;
            (*p_arch) >> p_cell_population;

            // Check the cell population has been restored correctly
            TS_ASSERT_EQUALS(p_cell_population->rGetCells().size(), 2u);

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

            // Check the mesh has been restored correctly
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumNodes(), 7u);
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumElements(), 2u);
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumFaces(), 0u);

            // Compare the loaded mesh against the original
            MutableVertexMesh<2,2>& loaded_mesh = p_cell_population->rGetMesh();

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), loaded_mesh.GetNumNodes());

            for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
            {
                Node<2>* p_node = mesh.GetNode(node_index);
                Node<2>* p_node2 = loaded_mesh.GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), true);
                TS_ASSERT_EQUALS(p_node2->IsBoundaryNode(), true);

                for (unsigned dimension=0; dimension<2; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(mesh.GetNumElements(), loaded_mesh.GetNumElements());

            for (unsigned elem_index=0; elem_index < mesh.GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(mesh.GetElement(elem_index)->GetNumNodes(),
                                 loaded_mesh.GetElement(elem_index)->GetNumNodes());

                for (unsigned local_index=0; local_index<mesh.GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    unsigned this_index = mesh.GetElement(elem_index)->GetNodeGlobalIndex(local_index);
                    unsigned loaded_index = loaded_mesh.GetElement(elem_index)->GetNodeGlobalIndex(local_index);

                    TS_ASSERT_EQUALS(this_index, loaded_index);
                }
            }

            TS_ASSERT_EQUALS(p_cell_population->GetOutputCellRearrangementLocations(), false);

            // Tidy up
            delete p_cell_population;
        }
    }

    void TestArchiving3dVertexBasedCellPopulation()
    {
        // Create mutable vertex mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 0.5, 0.5, 1.5));

        std::vector<std::vector<Node<3>*> > nodes_faces(10);
        nodes_faces[0].push_back(nodes[0]);
        nodes_faces[0].push_back(nodes[2]);
        nodes_faces[0].push_back(nodes[4]);
        nodes_faces[0].push_back(nodes[1]);
        nodes_faces[1].push_back(nodes[4]);
        nodes_faces[1].push_back(nodes[7]);
        nodes_faces[1].push_back(nodes[5]);
        nodes_faces[1].push_back(nodes[2]);
        nodes_faces[2].push_back(nodes[7]);
        nodes_faces[2].push_back(nodes[6]);
        nodes_faces[2].push_back(nodes[1]);
        nodes_faces[2].push_back(nodes[4]);
        nodes_faces[3].push_back(nodes[0]);
        nodes_faces[3].push_back(nodes[3]);
        nodes_faces[3].push_back(nodes[5]);
        nodes_faces[3].push_back(nodes[2]);
        nodes_faces[4].push_back(nodes[1]);
        nodes_faces[4].push_back(nodes[6]);
        nodes_faces[4].push_back(nodes[3]);
        nodes_faces[4].push_back(nodes[0]);
        nodes_faces[5].push_back(nodes[7]);
        nodes_faces[5].push_back(nodes[6]);
        nodes_faces[5].push_back(nodes[3]);
        nodes_faces[5].push_back(nodes[5]);
        nodes_faces[6].push_back(nodes[6]);
        nodes_faces[6].push_back(nodes[7]);
        nodes_faces[6].push_back(nodes[8]);
        nodes_faces[7].push_back(nodes[6]);
        nodes_faces[7].push_back(nodes[8]);
        nodes_faces[7].push_back(nodes[3]);
        nodes_faces[8].push_back(nodes[3]);
        nodes_faces[8].push_back(nodes[8]);
        nodes_faces[8].push_back(nodes[5]);
        nodes_faces[9].push_back(nodes[5]);
        nodes_faces[9].push_back(nodes[8]);
        nodes_faces[9].push_back(nodes[7]);

        std::vector<VertexElement<2,3>*> faces;
        for (unsigned i=0; i<10; i++)
        {
            faces.push_back(new VertexElement<2,3>(i, nodes_faces[i]));
        }

        std::vector<VertexElement<2,3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;
        for (unsigned i=0; i<6; i++)
        {
            faces_element_0.push_back(faces[i]);
            orientations_0.push_back(true);
        }
        for (unsigned i=6; i<10; i++)
        {
            faces_element_1.push_back(faces[i]);
            orientations_1.push_back(true);
        }
        faces_element_1.push_back(faces[5]);
        orientations_1.push_back(false);

        std::vector<VertexElement<3,3>*> elements;
        elements.push_back(new VertexElement<3,3>(0, faces_element_0, orientations_0));
        elements.push_back(new VertexElement<3,3>(1, faces_element_1, orientations_1));

        MutableVertexMesh<3,3> mesh(nodes, elements);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vertex_cell_population_3d.arch";
        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("vertex_mesh_population_3d");

        // Archive cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumElements());

            // Create cell population
            VertexBasedCellPopulation<3>* const p_cell_population = new VertexBasedCellPopulation<3>(mesh, cells);

            // Cells have been given birth times of 0 and -1.
            // Loop over them to run to time 0.0;
            for (AbstractCellPopulation<3>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&> (*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        // Restore cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore the cell population
            (*p_arch) >> *p_simulation_time;
            VertexBasedCellPopulation<3>* p_cell_population;
            (*p_arch) >> p_cell_population;

            // Check the cell population has been restored correctly
            TS_ASSERT_EQUALS(p_cell_population->rGetCells().size(), 2u);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;
            for (AbstractCellPopulation<3>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the mesh has been restored correctly
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumNodes(), 9u);
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumElements(), 2u);
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumFaces(), 10u);

            // Compare the loaded mesh against the original
            MutableVertexMesh<3,3>& loaded_mesh = p_cell_population->rGetMesh();

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), loaded_mesh.GetNumNodes());

            for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
            {
                Node<3>* p_node = mesh.GetNode(node_index);
                Node<3>* p_node2 = loaded_mesh.GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned dimension=0; dimension<2; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(mesh.GetNumElements(), loaded_mesh.GetNumElements());

            for (unsigned elem_index=0; elem_index < mesh.GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(mesh.GetElement(elem_index)->GetNumNodes(),
                                 loaded_mesh.GetElement(elem_index)->GetNumNodes());

                for (unsigned local_index=0; local_index<mesh.GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(mesh.GetElement(elem_index)->GetNodeGlobalIndex(local_index),
                                     loaded_mesh.GetElement(elem_index)->GetNodeGlobalIndex(local_index));
                }
            }

            // Tidy up
            delete p_cell_population;
        }
    }

    void TestGetTetrahedralMeshForPdeModifier()
    {
        /* Create a simple VertexMesh comprising three VertexElements.
         *
         * Note: central node is not a boundary node.
         *     ____
         *    |\  /|
         *    | \/ |
         *    |  \ |
         *    |___\|
         *
         */

        std::vector<Node<2>*> vertex_nodes;
        vertex_nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        vertex_nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        vertex_nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        vertex_nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
        vertex_nodes.push_back(new Node<2>(4, false, 0.5, 0.5));

        std::vector<std::vector<Node<2>*> > nodes_elements(3);
        nodes_elements[0].push_back(vertex_nodes[0]);
        nodes_elements[0].push_back(vertex_nodes[1]);
        nodes_elements[0].push_back(vertex_nodes[4]);
        nodes_elements[0].push_back(vertex_nodes[3]);
        nodes_elements[1].push_back(vertex_nodes[1]);
        nodes_elements[1].push_back(vertex_nodes[2]);
        nodes_elements[1].push_back(vertex_nodes[4]);
        nodes_elements[2].push_back(vertex_nodes[4]);
        nodes_elements[2].push_back(vertex_nodes[2]);
        nodes_elements[2].push_back(vertex_nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elements[0]));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elements[1]));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elements[2]));

        MutableVertexMesh<2,2>* p_vertex_mesh = new MutableVertexMesh<2,2>(vertex_nodes, vertex_elements);

        // Create cells to correspond with VertexElements
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_vertex_mesh->GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_vertex_mesh, cells);

        // Test TetrahedralMeshForPdeModifier() method
        TetrahedralMesh<2,2>* p_tetrahedral_mesh = cell_population.GetTetrahedralMeshForPdeModifier();

        // The VertexMesh has 5 Nodes and 3 VertexElements, so the TetrahedralMesh has 5+3=8 Nodes
        TS_ASSERT_EQUALS(p_tetrahedral_mesh->GetNumNodes(), 8u);

        // The First VertexElements comprise 4, 3, and 3 Nodes respectively, so the TetrahedralMesh has 4+3+3=10 Elements
        TS_ASSERT_EQUALS(p_tetrahedral_mesh->GetNumElements(), 10u);

        // Check the correct number of boundary elements (the same as the vertex mesh but we don't store boundary edges for Vertexmeshes)
        TS_ASSERT_EQUALS(p_tetrahedral_mesh->GetNumBoundaryElements(), 4u);

        // The first 5 Nodes of the TetrahedralMesh should overlap with those of the VertexMesh and have the same boundaryness.
        for (unsigned i=0; i<5; i++)
        {
            Node<2>* p_tetrahedral_node = p_tetrahedral_mesh->GetNode(i);
            Node<2>* p_vertex_node = p_vertex_mesh->GetNode(i);

            TS_ASSERT_EQUALS(p_tetrahedral_node->GetIndex(), p_vertex_node->GetIndex());
            TS_ASSERT_DELTA(p_tetrahedral_node->rGetLocation()[0], p_vertex_node->rGetLocation()[0], 1e-3);
            TS_ASSERT_DELTA(p_tetrahedral_node->rGetLocation()[1], p_vertex_node->rGetLocation()[1], 1e-3);

            TS_ASSERT_EQUALS(p_tetrahedral_node->IsBoundaryNode(), p_vertex_node->IsBoundaryNode());
        }

        // The last 3 Nodes of the TetrahedralMesh should be located at the centroids of the 3 VertexElements and will all be non boundary nodes.
        for (unsigned i=0; i<3; i++)
        {
            Node<2>* p_tetrahedral_node = p_tetrahedral_mesh->GetNode(i+5);
            c_vector<double,2> tetrahedral_node_location;
            tetrahedral_node_location = p_tetrahedral_node->rGetLocation();
            c_vector<double,2> vertex_element_centroid = p_vertex_mesh->GetCentroidOfElement(i);

            TS_ASSERT_DELTA(tetrahedral_node_location[0], vertex_element_centroid[0], 1e-3);
            TS_ASSERT_DELTA(tetrahedral_node_location[1], vertex_element_centroid[1], 1e-3);

            TS_ASSERT_EQUALS(p_tetrahedral_node->IsBoundaryNode(), false);
        }

        // The first Element of the TetrahedralMesh should contain the following Nodes
        Element<2,2>* p_element_0 = p_tetrahedral_mesh->GetElement(0);

        TS_ASSERT_EQUALS(p_element_0->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(p_element_0->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_element_0->GetNodeGlobalIndex(2), 5u);

        // The Fifth Element of the TetrahedralMesh should contain the following Nodes
        Element<2,2>* p_element_4 = p_tetrahedral_mesh->GetElement(4);

        TS_ASSERT_EQUALS(p_element_4->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(p_element_4->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(p_element_4->GetNodeGlobalIndex(2), 6u);

        // Avoid memory leaks
        delete p_vertex_mesh;
        delete p_tetrahedral_mesh;
    }

    void TestGetCellDataItemAtPdeNode()
    {
        // Create a small cell population
        HoneycombVertexMeshGenerator generator(4, 4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 48u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 16u);

        std::string var_name = "foo";
        TS_ASSERT_THROWS_THIS(cell_population.GetCellDataItemAtPdeNode(0,var_name),
            "The item foo is not stored");

        // Set the cell data item "foo" on each cell
        for (unsigned index=0; index<cell_population.GetNumElements(); index++)
        {
            double value = (double)index;
            cell_population.GetCellUsingLocationIndex(index)->GetCellData()->SetItem(var_name, value);
        }

        // The PDE mesh nodes 48 to 60 coincide with the cell centroids
        for (unsigned index=48; index<60; index++)
        {
            double expected_value = (double)index - 48.0;
            TS_ASSERT_DELTA(cell_population.GetCellDataItemAtPdeNode(index,var_name), expected_value, 1e-6);
        }

        // PDE mesh node 11 interpolates the value of "foo" from cells 1, 2 and 5
        double expected_value_11 = (1.0 + 2.0 + 5.0)/3.0;
        TS_ASSERT_DELTA(cell_population.GetCellDataItemAtPdeNode(11,var_name), expected_value_11, 1e-6);
    }
};

#endif /*TESTVERTEXBASEDCELLPOPULATION_HPP_*/
