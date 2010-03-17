/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTVERTEXBASEDTISSUE_HPP_
#define TESTVERTEXBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexBasedTissue.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "WildTypeCellMutationState.hpp"
#include "LabelledCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

class TestVertexBasedTissue : public AbstractCellBasedTestSuite
{
private:

    /**
     * Set up cells, one for each VertexElement.
     * Give each cell a birth time of -elem_index,
     * so its age is elem_index.
     */
    template<unsigned DIM>
    std::vector<TissueCell> SetUpCells(MutableVertexMesh<DIM,DIM>& rMesh)
    {
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<rMesh.GetNumElements(); i++)
        {
            TissueCell cell(DIFFERENTIATED, p_state, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        return cells;
    }

public:

    // Test construction, accessors and iterator
    void TestCreateSmallVertexBasedTissue() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(*p_mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(tissue.GetNumElements(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), cells.size());

        unsigned counter = 0;

        // Test VertexBasedTissue::Iterator
        for (VertexBasedTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            TS_ASSERT_EQUALS(counter, p_mesh->GetElement(counter)->GetIndex());

            counter++;
        }

        // Test we have gone through all cells in the for loop
        TS_ASSERT_EQUALS(counter, tissue.GetNumRealCells());

        // Test GetNumNodes() method
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), p_mesh->GetNumNodes());
    }


    void TestValidate() throw (Exception)
    {
        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each element.
        // Give each a birth time of -element_index, so the age = element_index.
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumElements()-1; i++)
        {
            TissueCell cell(STEM, p_state, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            cell_location_indices.push_back(i);
        }

        // This should throw an exception as the number of cells
        // does not equal the number of elements
        TS_ASSERT_THROWS_THIS(VertexBasedTissue<2> tissue(*p_mesh, cells),
                "Element 8 does not appear to have a cell associated with it");

        TissueCell cell(STEM, p_state, new FixedDurationGenerationBasedCellCycleModel());
        double birth_time = 0.0 - p_mesh->GetNumElements()-1;
        cell.SetBirthTime(birth_time);
        cells.push_back(cell);
        cell_location_indices.push_back(p_mesh->GetNumElements()-1);

        // This should pass as the number of cells equals the number of elements
        TS_ASSERT_THROWS_NOTHING(VertexBasedTissue<2> tissue(*p_mesh, cells));

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

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
            TissueCell& r_cell = tissue.rGetCellUsingLocationIndex(elem_index);
            VertexElement<2,2>* p_actual_element = tissue.GetElementCorrespondingToCell(r_cell);
            unsigned actual_index = p_actual_element->GetIndex();

            for (unsigned i=0; i<p_actual_element->GetNumNodes(); i++)
            {
                actual_node_indices.insert(p_actual_element->GetNodeGlobalIndex(i));
            }

            TS_ASSERT_EQUALS(actual_index, expected_index);
            TS_ASSERT_EQUALS(actual_node_indices, expected_node_indices);
        }
    }


    void TestGetTargetAreaOfCell() throw (Exception)
    {
        double apoptosis_time = TissueConfig::Instance()->GetApoptosisTime();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.5*apoptosis_time, 3);

        // Create mesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            CellProliferativeType cell_type = STEM;

            if ((i==0) || (i==4))
            {
                cell_type = DIFFERENTIATED;
            }

            TissueCell cell(cell_type, p_state, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - 2*i;

            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            cell_location_indices.push_back(i);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);
        tissue.InitialiseCells(); // this method must be called explicitly as there is no simulation

        // Check GetTargetAreaOfCell()
        for (VertexMesh<2,2>::VertexElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            TissueCell cell = tissue.rGetCellUsingLocationIndex(elem_index);
            double expected_area = TissueConfig::Instance()->GetMatureCellTargetArea();

            if (elem_index!=4 && elem_index<=7u)
            {
                expected_area *= 0.5*(1 + ((double)elem_index)/7.0);
            }

            double actual_area = tissue.GetTargetAreaOfCell(cell);

            TS_ASSERT_DELTA(actual_area, expected_area, 1e-12);
        }

        TissueCell cell_0 = tissue.rGetCellUsingLocationIndex(0);
        TissueCell cell_1 = tissue.rGetCellUsingLocationIndex(1);
        TissueCell cell_4 = tissue.rGetCellUsingLocationIndex(4);

        // Make cell 1 and 4 undergo apoptosis
        cell_1.StartApoptosis();
        cell_4.StartApoptosis();

        double actual_area_0 = tissue.GetTargetAreaOfCell(cell_0);
        double actual_area_1 = tissue.GetTargetAreaOfCell(cell_1);
        double actual_area_4 = tissue.GetTargetAreaOfCell(cell_4);

        double expected_area_0 = 0.5;
        double expected_area_1 = TissueConfig::Instance()->GetMatureCellTargetArea()*0.5*(1.0 + 1.0/7.0);
        double expected_area_4 = TissueConfig::Instance()->GetMatureCellTargetArea();

        TS_ASSERT_DELTA(actual_area_0, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1, expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_4, expected_area_4, 1e-12);

        // Run for one time step
        p_simulation_time->IncrementTimeOneStep();

        double actual_area_0_after_dt = tissue.GetTargetAreaOfCell(cell_0);
        double actual_area_1_after_dt = tissue.GetTargetAreaOfCell(cell_1);
        double actual_area_4_after_dt = tissue.GetTargetAreaOfCell(cell_4);

        // The target areas of cells 1 and 4 should have halved
        expected_area_0 = TissueConfig::Instance()->GetMatureCellTargetArea()*0.5*(1.0 + 0.5*apoptosis_time/2.0);

        TS_ASSERT_DELTA(actual_area_0_after_dt, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_dt, 0.5*expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_dt, 0.5*expected_area_4, 1e-12);

        // Make cell 0 undergo apoptosis
        cell_0.StartApoptosis();

        // Now run on for a further time step
        p_simulation_time->IncrementTimeOneStep();

        double actual_area_0_after_2dt = tissue.GetTargetAreaOfCell(cell_0);
        double actual_area_1_after_2dt = tissue.GetTargetAreaOfCell(cell_1);
        double actual_area_4_after_2dt = tissue.GetTargetAreaOfCell(cell_4);

        // Cells 1 and 4 should now have zero target area and the target area of cell 0 should have halved
        TS_ASSERT_DELTA(actual_area_0_after_2dt, 0.5*expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_2dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_2dt, 0.0, 1e-12);

        // Now run on for even further, for coverage
        p_simulation_time->IncrementTimeOneStep();

        double actual_area_0_after_3dt = tissue.GetTargetAreaOfCell(cell_0);
        double actual_area_1_after_3dt = tissue.GetTargetAreaOfCell(cell_1);
        double actual_area_4_after_3dt = tissue.GetTargetAreaOfCell(cell_4);

        // All apoptotic cells should now have zero target area
        TS_ASSERT_DELTA(actual_area_0_after_3dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_3dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_3dt, 0.0, 1e-12);
    }


    void TestGetDampingConstant()
    {
        TissueConfig::Instance()->SetDampingConstantMutant(8.0);

        // Create mesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_labelled(new LabelledCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            TissueCell cell(STEM, p_state, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            cell_location_indices.push_back(i);
        }
        cells[0].SetMutationState(p_apc1);
        cells[6].SetMutationState(p_apc2);
        cells[7].SetMutationState(p_bcat1);
        cells[8].SetMutationState(p_labelled);

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);
        tissue.InitialiseCells(); // this method must be called explicitly as there is no simulation

        // Test GetDampingConstant()
        double normal_damping_constant = TissueConfig::Instance()->GetDampingConstantNormal();
        double mutant_damping_constant = TissueConfig::Instance()->GetDampingConstantMutant();

        // Node 0 is contained in cell 0 only, therefore should have mutant damping constant
        double damping_constant_at_node_0 = tissue.GetDampingConstant(0);
        TS_ASSERT_DELTA(damping_constant_at_node_0, mutant_damping_constant, 1e-6);

        // Node 1 is contained in cell 2 only, therefore should have a normal damping constant
        double damping_constant_at_node_1 = tissue.GetDampingConstant(1);
        TS_ASSERT_DELTA(damping_constant_at_node_1, normal_damping_constant, 1e-6);

        // Node 4 is contained in cells 0 and 1, therefore should an averaged damping constant
        double damping_constant_at_node_4 = tissue.GetDampingConstant(4);
        TS_ASSERT_DELTA(damping_constant_at_node_4, (normal_damping_constant+mutant_damping_constant)/2.0, 1e-6);

        // Node 8 is contained in cells 0, 1, 3, therefore should an averaged damping constant
        double damping_constant_at_node_8 = tissue.GetDampingConstant(8);
        TS_ASSERT_DELTA(damping_constant_at_node_8, (2*normal_damping_constant+mutant_damping_constant)/3.0, 1e-6);

        // Node 27 is contained in cell 6 only, therefore should have a mutant damping constant
        double damping_constant_at_node_27 = tissue.GetDampingConstant(27);
        TS_ASSERT_DELTA(damping_constant_at_node_27, mutant_damping_constant, 1e-6);

        // Node 25 is contained in cells 6 and 7, therefore should have a mutant damping constant
        double damping_constant_at_node_25 = tissue.GetDampingConstant(25);
        TS_ASSERT_DELTA(damping_constant_at_node_25, mutant_damping_constant, 1e-6);
    }


    void TestUpdateWithoutBirthOrDeath() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(*p_mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        unsigned num_cells_removed = tissue.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_cells_removed, 0u);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_THROWS_NOTHING(tissue.Update());
    }


    void TestAddCellWithSimpleMesh() throw (Exception)
    {
        // Make some nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, false, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 2.0));

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
        std::vector<TissueCell> cells = SetUpCells(vertex_mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(vertex_mesh, cells);

        // For coverage, test GetLocationOfCellCentre()

        // Cell 0 is a rectangle with centre of mass (0,0)
        c_vector<double, 2> cell0_location = tissue.GetLocationOfCellCentre(tissue.rGetCellUsingLocationIndex(0));
        TS_ASSERT_DELTA(cell0_location[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell0_location[1], 0.0, 1e-4);

        // Cell 1 is a triangle with centre of mass (0,4/3)
        c_vector<double, 2> cell1_location = tissue.GetLocationOfCellCentre(tissue.rGetCellUsingLocationIndex(1));
        TS_ASSERT_DELTA(cell1_location[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell1_location[1], 4.0/3.0, 1e-4);

        unsigned old_num_nodes = vertex_mesh.GetNumNodes();
        unsigned old_num_elements = vertex_mesh.GetNumElements();
        unsigned old_num_cells = tissue.rGetCells().size();

        // Add new cell by dividing element 0 along short axis

        c_vector<double,2> cell_division_vector = zero_vector<double>(2);

        TissueCell cell0 = tissue.rGetCellUsingLocationIndex(0);
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TissueCell new_cell(STEM, p_state, new FixedDurationGenerationBasedCellCycleModel());
        new_cell.SetBirthTime(-1);

        TissueCell* p_new_cell = tissue.AddCell(new_cell, cell_division_vector, &cell0);

        // Check that the new cell was successfully added to the tissue
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), old_num_nodes+2);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), old_num_elements+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), old_num_elements+1);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[1], 1.0, 1e-12);

        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[1], -1.0, 1e-12);

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(3), 5u);

        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(3), 6u);

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);
        expected_elements_containing_node_5.insert(2);

        TS_ASSERT_EQUALS(tissue.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(2);

        TS_ASSERT_EQUALS(tissue.GetNode(6)->rGetContainingElementIndices(), expected_elements_containing_node_6);

        // Check the index of the new cell
        TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(*p_new_cell), old_num_elements);
    }


    void TestAddCellWithGivenDivisionVector() throw (Exception)
    {
        // Make a vertex mesh consisting of a single square element
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem.push_back(nodes[i]);
        }
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Create a cell
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TissueCell cell(DIFFERENTIATED, p_state, new FixedDurationGenerationBasedCellCycleModel());
        cell.SetBirthTime(-20.0);
        cells.push_back(cell);

        // Create tissue
        VertexBasedTissue<2> tissue(vertex_mesh, cells);

        TS_ASSERT_EQUALS(tissue.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 1u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 1u);

        // Add a new cell by dividing element 0 along the axis (1,0)
        c_vector<double,2> cell_division_axis;
        cell_division_axis[0] = 1.0;
        cell_division_axis[1] = 0.0;

        TissueCell new_cell(STEM, p_state, new FixedDurationGenerationBasedCellCycleModel());
        new_cell.SetBirthTime(-1.0);

        TissueCell* p_new_cell = tissue.AddCell(new_cell, cell_division_axis, &cell);

        // Check that the new cell was successfully added to the tissue
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 2u);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(tissue.GetNode(4)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(tissue.GetNode(4)->rGetLocation()[1], 0.5, 1e-12);
        TS_ASSERT_DELTA(tissue.GetNode(5)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(tissue.GetNode(5)->rGetLocation()[1], 0.5, 1e-12);

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_4;
        expected_elements_containing_node_4.insert(0);
        expected_elements_containing_node_4.insert(1);
        TS_ASSERT_EQUALS(tissue.GetNode(4)->rGetContainingElementIndices(), expected_elements_containing_node_4);

        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);
        TS_ASSERT_EQUALS(tissue.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        // Check the index of the new cell
        TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(*p_new_cell), 1u);
    }

    void TestAddCellWithHoneycombMesh() throw (Exception)
    {
        // Create a mesh with 9 elements
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            CellProliferativeType cell_type = DIFFERENTIATED;
            double birth_time = 0.0 - elem_index;

            // Cell 4 should divide immediately
            if (elem_index==4)
            {
                cell_type = TRANSIT; // As stem cells always divide horizontally
                birth_time = -50.0;
            }

            TissueCell cell(cell_type, p_state, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Initialise cells (usually called in tissue simulation constructor)
        tissue.InitialiseCells();

        unsigned old_num_nodes = p_mesh->GetNumNodes();
        unsigned old_num_elements = p_mesh->GetNumElements();
        unsigned old_num_cells = tissue.rGetCells().size();

        // Add a new cell by dividing cell 4

        tissue.rGetCellUsingLocationIndex(4).ReadyToDivide();
        TissueCell& cell4 = tissue.rGetCellUsingLocationIndex(4);

        TissueCell new_cell = tissue.rGetCellUsingLocationIndex(4).Divide();

        c_vector<double, 2> new_location = zero_vector<double>(2);

        // Add new cell to the tissue
        tissue.AddCell(new_cell, new_location, &cell4);

        // Check that the new cell was successfully added to the tissue
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), old_num_nodes+2);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), old_num_elements+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), old_num_elements+1);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[0], 2.4721, 1e-4);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[1], 1.7481, 1e-4);

        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[0], 1.5278, 1e-4);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[1], 1.1386, 1e-4);

        // Now test the nodes in each element
        for (unsigned i=0; i<tissue.GetNumElements(); i++)
        {
            if (i==4 || i==9)
            {
                // Elements 4 and 9 should each have one less node
                TS_ASSERT_EQUALS(tissue.GetElement(i)->GetNumNodes(), 5u);
            }
            else if (i==1 || i==8)
            {
                // Elements 3 and 8 should each have one extra node
                TS_ASSERT_EQUALS(tissue.GetElement(i)->GetNumNodes(), 7u);
            }
            else
            {
                TS_ASSERT_EQUALS(tissue.GetElement(i)->GetNumNodes(), 6u);
            }
        }

        // Check node ownership for a few elements

        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(2), 8u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(3), 11u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(4), 7u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(5), 3u);

        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(0), 9u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(1), 13u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(2), 17u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(3), 30u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(4), 31u);

        TS_ASSERT_EQUALS(tissue.GetElement(8)->GetNodeGlobalIndex(0), 17u);
        TS_ASSERT_EQUALS(tissue.GetElement(8)->GetNodeGlobalIndex(1), 22u);
        TS_ASSERT_EQUALS(tissue.GetElement(8)->GetNodeGlobalIndex(2), 26u);
        TS_ASSERT_EQUALS(tissue.GetElement(8)->GetNodeGlobalIndex(3), 29u);
        TS_ASSERT_EQUALS(tissue.GetElement(8)->GetNodeGlobalIndex(4), 25u);
        TS_ASSERT_EQUALS(tissue.GetElement(8)->GetNodeGlobalIndex(5), 21u);
        TS_ASSERT_EQUALS(tissue.GetElement(8)->GetNodeGlobalIndex(6), 30u);

        // Test element ownership for a few nodes

        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(1);
        expected_elements_containing_node_5.insert(2);

        TS_ASSERT_EQUALS(tissue.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_13;
        expected_elements_containing_node_13.insert(2);
        expected_elements_containing_node_13.insert(4);
        expected_elements_containing_node_13.insert(5);

        TS_ASSERT_EQUALS(tissue.GetNode(13)->rGetContainingElementIndices(), expected_elements_containing_node_13);

        std::set<unsigned> expected_elements_containing_node_30;
        expected_elements_containing_node_30.insert(8);
        expected_elements_containing_node_30.insert(4);
        expected_elements_containing_node_30.insert(9);

        TS_ASSERT_EQUALS(tissue.GetNode(30)->rGetContainingElementIndices(), expected_elements_containing_node_30);
    }

    void TestIsCellAssociatedWithADeletedLocation() throw (Exception)
    {
        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->GetElement(5)->MarkAsDeleted();

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(*p_mesh);

        // Create tissue but do not try to validate
        VertexBasedTissue<2> tissue(*p_mesh, cells, false, false);

        // Test IsCellAssociatedWithADeletedLocation() method
        for (VertexBasedTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            bool is_deleted = tissue.IsCellAssociatedWithADeletedLocation(*cell_iter);

            if (tissue.GetLocationIndexUsingCell(*cell_iter) == 5)
            {
                TS_ASSERT_EQUALS(is_deleted, true);
            }
            else
            {
                TS_ASSERT_EQUALS(is_deleted, false);
            }
        }
    }

    void TestRemoveDeadCellsAndUpdate() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(*p_mesh);
        cells[5].StartApoptosis();

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 24u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 24u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 24u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 68u);

        p_simulation_time->IncrementTimeOneStep();

        // Remove dead cells
        unsigned num_cells_removed = tissue.RemoveDeadCells();

        TS_ASSERT_EQUALS(num_cells_removed, 1u);

        // We should now have one less real cell, since one cell has been
        // marked as dead, so is skipped by the tissue iterator
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), 68u);

        tissue.Update();

        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), 23u);
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), 68u);

        // Finally, check the cells' element indices have updated

        // We expect the cell element indices to be {0,11,...,23}
        std::set<unsigned> expected_elem_indices;
        for (unsigned i=0; i<tissue.GetNumRealCells(); i++)
        {
            expected_elem_indices.insert(i);
        }

        // Get actual cell element indices
        std::set<unsigned> element_indices;

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Record element index corresponding to cell
            unsigned element_index = tissue.GetLocationIndexUsingCell(*cell_iter);
            element_indices.insert(element_index);
        }

        TS_ASSERT_EQUALS(element_indices, expected_elem_indices);
    }


    void TestVertexBasedTissueOutputWriters() throw (Exception)
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(*p_mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // For coverage of WriteResultsToFiles()
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_labelled(new LabelledCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc1(new ApcOneHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        tissue.rGetCellUsingLocationIndex(0).SetCellProliferativeType(TRANSIT);
        tissue.rGetCellUsingLocationIndex(0).SetMutationState(p_labelled);
        tissue.rGetCellUsingLocationIndex(1).SetCellProliferativeType(DIFFERENTIATED);
        tissue.rGetCellUsingLocationIndex(1).SetMutationState(p_apc1);
        tissue.rGetCellUsingLocationIndex(2).SetMutationState(p_apc2);
        tissue.rGetCellUsingLocationIndex(3).SetMutationState(p_bcat1);
        tissue.rGetCellUsingLocationIndex(4).SetCellProliferativeType(APOPTOTIC);
        tissue.rGetCellUsingLocationIndex(4).StartApoptosis();
        tissue.rGetCellUsingLocationIndex(5).SetCellProliferativeType(STEM);
        tissue.SetCellAncestorsToLocationIndices();

        std::string output_directory = "TestVertexBasedTissueOutputWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        TissueConfig::Instance()->SetOutputCellMutationStates(true);
        TissueConfig::Instance()->SetOutputCellProliferativeTypes(true);
        TissueConfig::Instance()->SetOutputCellCyclePhases(true);
        TissueConfig::Instance()->SetOutputCellAncestors(true);
        TissueConfig::Instance()->SetOutputCellAges(true);

        TS_ASSERT_THROWS_NOTHING(tissue.CreateOutputFiles(output_directory, false));

        tissue.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(tissue.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/cellages.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellcyclephases.dat     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/cellcyclephases.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "celltypes.dat     notforrelease_cell_based/test/data/TestVertexBasedTissueOutputWriters/celltypes.dat").c_str()), 0);

        // For coverage
        TS_ASSERT_THROWS_NOTHING(tissue.WriteResultsToFiles());
    }


    // At the moment the tissue cannot be properly archived since the mesh cannot be. This test
    // just checks that the cells are correctly archived.
    void TestArchivingVertexBasedTissue() throw(Exception)
    {
        std::string archive_dir = "archive";
        std::string archive_file = "vertex_tissue.arch";
        // The following line is required because the loading of a tissue
        // is usually called by the method TissueSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("vertex_mesh");

        // Archive tissue
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create mesh
            VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh");
            MutableVertexMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells
            std::vector<TissueCell> cells = SetUpCells(mesh);

            // Create tissue
            VertexBasedTissue<2>* const p_tissue = new VertexBasedTissue<2>(mesh, cells);

            // Cells have been given birth times of 0 and -1.
            // Loop over them to run to time 0.0;
            for (AbstractTissue<2>::Iterator cell_iter = p_tissue->Begin();
                 cell_iter != p_tissue->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the tissue to the archive
            (*p_arch) << static_cast<const SimulationTime&> (*p_simulation_time);
            (*p_arch) << p_tissue;

            // Tidy up
            SimulationTime::Destroy();
            delete p_tissue;
        }

        // Restore tissue
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

            // Restore the tissue
            (*p_arch) >> *p_simulation_time;
            VertexBasedTissue<2>* p_tissue;
            (*p_arch) >> p_tissue;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the tissue has been restored
            TS_ASSERT_EQUALS(p_tissue->rGetCells().size(), 2u);

            // Tidy up
            delete p_tissue;
        }
    }


    void TestUpdateNodeLocations() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(*p_mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Make up some forces
        std::vector<c_vector<double, 2> > old_posns(tissue.GetNumNodes());
        std::vector<c_vector<double, 2> > forces_on_nodes(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            old_posns[i][0] = p_mesh->GetNode(i)->rGetLocation()[0];
            old_posns[i][1] = p_mesh->GetNode(i)->rGetLocation()[1];

            forces_on_nodes[i][0] = i*0.01;
            forces_on_nodes[i][1] = 2*i*0.01;
        }

        double time_step = 0.01;

        tissue.UpdateNodeLocations(forces_on_nodes, time_step);

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
        }
    }


    /**
     * Test that post-#878, WntConcentration copes with a VertexBasedTissue.
     * \todo When vertex-based tissue code is added to cell_based folder, move this
     *       test to TestWntConcentration.hpp
     */
    void TestWntConcentrationWithVertexBasedTissue() throw(Exception)
    {
        // Make some nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, false, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 2.0));

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

        // Create cells
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<vertex_mesh.GetNumElements(); i++)
        {
        	WntCellCycleModel* p_cell_cycle_model = new WntCellCycleModel;
        	p_cell_cycle_model->SetDimension(2);
            TissueCell cell(DIFFERENTIATED, p_state, p_cell_cycle_model);
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(vertex_mesh, cells);

        // Set the top of this tissue, for the purposes of computing the WntConcentration
        TissueConfig::Instance()->SetCryptLength(4.0);

        // Set up an instance of the WntConcentration singleton object
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);

        // Check that the singleton can be set up
        p_wnt->SetType(LINEAR);
        p_wnt->SetTissue(tissue);

        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        // Check that the singleton can be destroyed then recreated
        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(NONE);
        WntConcentration<2>::Instance()->SetTissue(tissue);
        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), false); // not fully set up now it is a NONE type

        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(tissue);

        p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true); // set up again

        double wnt_at_cell0 = p_wnt->GetWntLevel(tissue.rGetCellUsingLocationIndex(0));
        double wnt_at_cell1 = p_wnt->GetWntLevel(tissue.rGetCellUsingLocationIndex(1));

        // We have set the top of the tissue to be 4, so the WntConcentration should decrease linearly
        // up the tissue, from one at height 0 to zero at height 4.

        // Cell 0 has centre of mass (0,0)
        TS_ASSERT_DELTA(wnt_at_cell0, 1.0, 1e-4);

        // Cell 1 has centre of mass (0, 4/3)
        TS_ASSERT_DELTA(wnt_at_cell1, 2.0/3.0, 1e-4);
    }

};

#endif /*TESTVERTEXBASEDTISSUE_HPP_*/
