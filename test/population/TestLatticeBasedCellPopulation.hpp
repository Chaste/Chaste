/*

Copyright (C) University of Oxford, 2005-2011

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
#ifndef TESTLATTICEBASEDCELLPOPULATION_HPP_
#define TESTLATTICEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "LatticeBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class TestLatticeBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    /**
     * Helper method. Create a single, wild type, differentiated cell and
     * return as a vector for passing into a cell population constructor.
     *
     * @param cellProliferativeType  the proliferative type of the cell (defaults to DIFFERENTIATED)
     */
    std::vector<CellPtr> CreateSingleCellPtr(CellProliferativeType cellProliferativeType=DIFFERENTIATED)
    {
        std::vector<CellPtr> cells;

        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(cellProliferativeType);
        CellPtr p_cell(new Cell(p_state, p_model));

        cells.push_back(p_cell);
        return cells;
    }

    /**
     * Helper method. Return a given location index corresponding to a single
     * cell as a vector for passing into a cell population constructor.
     *
     * @param locationIndex the locaiton index
     */
    std::vector<unsigned> CreateSingleLocationIndex(unsigned locationIndex)
    {
        std::vector<unsigned> location_indices;
        location_indices.push_back(locationIndex);
        return location_indices;
    }

public:

    void TestConstructorAndBasicMethods() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create two cells
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model_1 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_1->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_cell_1(new Cell(p_state, p_model_1));

        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_cell_2(new Cell(p_state, p_model_2));

        std::vector<CellPtr> cells;
        cells.push_back(p_cell_1);
        cells.push_back(p_cell_2);

        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);
        real_node_indices.push_back(1);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test GetNumNodes()
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 9u);

        // Test GetNode()
        for (unsigned index=0; index<cell_population.GetNumNodes(); index++)
        {
            Node<2>* p_node = cell_population.GetNode(index);
            TS_ASSERT_EQUALS(p_node->GetIndex(), index);

            c_vector<double, 2> node_location = p_node->rGetLocation();
            double expected_x = (double)(index%3);
            double expected_y = (double)(index>2) + (double)(index>5);
            TS_ASSERT_DELTA(node_location[0], expected_x, 1e-3);
            TS_ASSERT_DELTA(node_location[1], expected_y, 1e-3);
        }

        // Test GetNumRealCells()
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);

        // Test GetLocationOfCellCentre() and GetLocationIndexUsingCell()
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 4u);

        c_vector<double, 2> cell_1_location = cell_population.GetLocationOfCellCentre(*cell_iter);
        TS_ASSERT_DELTA(cell_1_location[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(cell_1_location[1], 1.0, 1e-6);

        ++cell_iter;

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 1u);

        c_vector<double, 2> cell_2_location = cell_population.GetLocationOfCellCentre(*cell_iter);
        TS_ASSERT_DELTA(cell_2_location[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(cell_2_location[1], 0.0, 1e-6);

        // Test SetEmptySite(), rGetEmptySites() and GetEmptySiteIndices()
        std::vector<bool>& r_ghost_nodes = cell_population.rGetEmptySites();
        for (unsigned index=0; index<cell_population.GetNumNodes(); index++)
        {
            bool expect_ghost_node = (index != 4) && (index != 1);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(index), expect_ghost_node);
            TS_ASSERT_EQUALS(r_ghost_nodes[index], expect_ghost_node);
        }

        std::set<unsigned> ghost_node_indices = cell_population.GetEmptySiteIndices();
        std::set<unsigned> expected_ghost_node_indices;
        for (unsigned index=0; index<cell_population.GetNumNodes(); index++)
        {
            if ((index != 4) && (index != 1))
            {
                expected_ghost_node_indices.insert(index);
            }
        }

        TS_ASSERT_EQUALS(ghost_node_indices, expected_ghost_node_indices);

        // Test rGetMesh()
        TetrahedralMesh<2,2>& r_mesh = cell_population.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(r_mesh.GetNumElements(), 8u);

        // Test constructor without location indices argument

        std::vector<CellPtr> cells2;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            CellPtr p_cell(new Cell(p_state, p_model));
            cells2.push_back(p_cell);
        }

        LatticeBasedCellPopulation<2> cell_population2(mesh, cells2);

        TS_ASSERT_EQUALS(cell_population2.GetNumRealCells(), 9u);

        r_ghost_nodes = cell_population2.rGetEmptySites();
        for (unsigned index=0; index<cell_population.GetNumNodes(); index++)
        {
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(index), false);
            TS_ASSERT_EQUALS(r_ghost_nodes[index], false);
        }
    }

    void TestCellDivision() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create a single cell and corresponding location index
        std::vector<CellPtr> cells = CreateSingleCellPtr(STEM);
        std::vector<unsigned> real_node_indices = CreateSingleLocationIndex(4);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);

        // Check locations of cells
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*(cell_population.rGetCells().begin())), 4u);

        // Create a new cell
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(STEM);
        CellPtr p_new_cell(new Cell(p_state, p_model_2));

        // Add new cell to the cell population by dividing the cell at node 4
        AbstractCellPopulation<2>::Iterator cell_iter_1 = cell_population.Begin();
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), *cell_iter_1);

        TS_ASSERT_LESS_THAN(cell_population.GetLocationIndexUsingCell(p_new_cell), cell_population.GetNumNodes());
        TS_ASSERT(cell_population.GetLocationIndexUsingCell(p_new_cell) != 4u);

        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);

        // Check locations of parent cell
        AbstractCellPopulation<2>::Iterator cell_iter_2 = cell_population.Begin();
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_2), 4u);
        ++cell_iter_2;
        TS_ASSERT_LESS_THAN(cell_population.GetLocationIndexUsingCell(*cell_iter_2), cell_population.GetNumNodes());
        TS_ASSERT(cell_population.GetLocationIndexUsingCell(*cell_iter_2) != 4u);
    }

    // Checks that this method correctly returns the first degree, i.e. nearest, neighbours
    void TestGetFirstDegreeNeighbouringNodeIndices() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create a single cell and corresponding location index
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> real_node_indices = CreateSingleLocationIndex(4);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        std::set<unsigned> nearest_neighbours = cell_population.GetNthDegreeNeighbouringNodeIndices(4, 1);

        std::set<unsigned> expected_neighbours;
        expected_neighbours.insert(7);
        expected_neighbours.insert(6);
        expected_neighbours.insert(3);
        expected_neighbours.insert(0);
        expected_neighbours.insert(1);
        expected_neighbours.insert(2);
        expected_neighbours.insert(5);
        expected_neighbours.insert(8);

        TS_ASSERT_EQUALS(nearest_neighbours.size(), expected_neighbours.size());
        TS_ASSERT_EQUALS(nearest_neighbours, expected_neighbours);

        std::set<unsigned>::iterator neighbour_iter = nearest_neighbours.begin();
        std::set<unsigned>::iterator expected_iter = expected_neighbours.begin();

        for (unsigned i=0; i<nearest_neighbours.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Now going to check for a corner node
        std::set<unsigned> nearest_neighbours_2 = cell_population.GetNthDegreeNeighbouringNodeIndices(2, 1);

        std::set<unsigned> expected_neighbours_2;
        expected_neighbours_2.insert(5u);
        expected_neighbours_2.insert(4u);
        expected_neighbours_2.insert(1u);

        TS_ASSERT_EQUALS(nearest_neighbours_2.size(), expected_neighbours_2.size());
        TS_ASSERT_EQUALS(nearest_neighbours_2, expected_neighbours_2);

        neighbour_iter = nearest_neighbours_2.begin();
        expected_iter = expected_neighbours_2.begin();

        for (unsigned i=0; i<nearest_neighbours_2.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Now checking a boundary node
        std::set<unsigned> nearest_neighbours_3 = cell_population.GetNthDegreeNeighbouringNodeIndices(3, 1);

        std::set<unsigned> expected_neighbours_3;
        expected_neighbours_3.insert(6u);
        expected_neighbours_3.insert(0u);
        expected_neighbours_3.insert(1u);
        expected_neighbours_3.insert(4u);
        expected_neighbours_3.insert(7u);

        TS_ASSERT_EQUALS(nearest_neighbours_3.size(), expected_neighbours_3.size());
        TS_ASSERT_EQUALS(nearest_neighbours_3, expected_neighbours_3);

        neighbour_iter = nearest_neighbours_3.begin();
        expected_iter = expected_neighbours_3.begin();

        for (unsigned i=0; i<nearest_neighbours_3.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }
    }

    void TestGetNthDegreeNeighbouringNodeIndices() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(6, 6, true); // 7*7 nodes

        // Create a single cell and corresponding location index
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> real_node_indices = CreateSingleLocationIndex(29);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Check the maximum degree possible in each direction
        std::vector<unsigned> max_degrees = cell_population.GetMaximumDegreeInEachDirection(29);
        TS_ASSERT_EQUALS(max_degrees.size(), 8u);
        TS_ASSERT_EQUALS(max_degrees[0], 2u);
        TS_ASSERT_EQUALS(max_degrees[1], 1u);
        TS_ASSERT_EQUALS(max_degrees[2], 1u);
        TS_ASSERT_EQUALS(max_degrees[3], 1u);
        TS_ASSERT_EQUALS(max_degrees[4], 4u);
        TS_ASSERT_EQUALS(max_degrees[5], 4u);
        TS_ASSERT_EQUALS(max_degrees[6], 5u);
        TS_ASSERT_EQUALS(max_degrees[7], 2u);

        // Check the second degree neighbours
        std::set<unsigned> nearest_neighbours_2 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 2);
        std::set<unsigned> expected_neighbours_2;
        expected_neighbours_2.insert(43);
        expected_neighbours_2.insert(15);
        expected_neighbours_2.insert(17);
        expected_neighbours_2.insert(31);
        expected_neighbours_2.insert(45);

        TS_ASSERT_EQUALS(nearest_neighbours_2.size(), expected_neighbours_2.size());
        TS_ASSERT_EQUALS(nearest_neighbours_2, expected_neighbours_2);

        std::set<unsigned>::iterator neighbour_iter = nearest_neighbours_2.begin();
        std::set<unsigned>::iterator expected_iter = expected_neighbours_2.begin();

        for (unsigned i=0; i<nearest_neighbours_2.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Check the third degree neighbours
        std::set<unsigned> nearest_neighbours_3 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 3);

        std::set<unsigned> expected_neighbours_3;
        expected_neighbours_3.insert(8);
        expected_neighbours_3.insert(11);
        expected_neighbours_3.insert(32);

        TS_ASSERT_EQUALS(nearest_neighbours_3.size(), expected_neighbours_3.size());
        TS_ASSERT_EQUALS(nearest_neighbours_3, expected_neighbours_3);

        neighbour_iter = nearest_neighbours_3.begin();
        expected_iter = expected_neighbours_3.begin();

        for (unsigned i=0; i<nearest_neighbours_3.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Check the fourth degree neighbours
        std::set<unsigned> nearest_neighbours_4 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 4);

        std::set<unsigned> expected_neighbours_4;
        expected_neighbours_4.insert(1);
        expected_neighbours_4.insert(5);
        expected_neighbours_4.insert(33);

        TS_ASSERT_EQUALS(nearest_neighbours_4.size(), expected_neighbours_4.size());
        TS_ASSERT_EQUALS(nearest_neighbours_4, expected_neighbours_4);

        neighbour_iter = nearest_neighbours_4.begin();
        expected_iter = expected_neighbours_4.begin();

        for (unsigned i=0; i<nearest_neighbours_4.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Check the fifth degree neighbours
        std::set<unsigned> nearest_neighbours_5 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 5);

        std::set<unsigned> expected_neighbours_5;
        expected_neighbours_5.insert(34);

        TS_ASSERT_EQUALS(nearest_neighbours_5.size(), expected_neighbours_5.size());
        TS_ASSERT_EQUALS(nearest_neighbours_5, expected_neighbours_5);

        neighbour_iter = nearest_neighbours_5.begin();
        expected_iter = expected_neighbours_5.begin();

        for (unsigned i=0; i<nearest_neighbours_5.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Check the sixth degree neighbours
        std::set<unsigned> nearest_neighbours_6 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 6);
        TS_ASSERT_EQUALS(nearest_neighbours_6.empty(), true);
    }

    void TestAddCellWithOneEmptySite() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(6, 6, true); // 7*7 nodes

        // Initially we create a cell for every node except node 45
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (i!=45)
            {
                real_node_indices.push_back(i);
            }
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        TS_ASSERT(cell_population.IsEmptySite(45));
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 49u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 48u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 48u);

        // Create a new cell
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_new_cell(new Cell(p_state, p_model_2));

        // Add new cell to the cell population by dividing the cell at node 24
        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(24);
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_parent_cell);

        // Check the number of cells is correct
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 49u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 49u);

        // Check each cell corresponds to the correct location index
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();

        // The indices of cells 0 to 30 should be unchanged
        for (unsigned i=0; i<31; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // Cell 31 should have been pushed up to index 38
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 38u);
        ++cell_iter;

        // The indices of cells 32 to 37 should be unchanged
        for (unsigned i=32; i<38; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // Cell 38 should have been pushed up to index 45
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 45u);
        ++cell_iter;

        // The indices of cells 39 to 44 should be unchanged
        for (unsigned i=39; i<45; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // The indices of cells 46 to 48 should be unchanged
        for (unsigned i=46; i<49; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // Lastly the new cell should be at index 31
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 31u);
    }

    void TestCellDivisionWithOneEmptySiteOnlySearchingNearestNeighbours() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(3, 3, true); // 4*4 nodes

        // Initially we create a cell for every node except node 13
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (i!=13)
            {
                real_node_indices.push_back(i);
            }
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices, true); // Only searching the nearest neighbours (degree=1)

        TS_ASSERT(cell_population.IsEmptySite(13));
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 15u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 15u);

        // Create a new cell
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_new_cell(new Cell(p_state, p_model_2));

        // Add new cell to the cell population by dividing the cell at node 5
        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(5);

        // Try to divide the parent - should not be possible as there are no free nearest neighbours

        TS_ASSERT_THROWS_THIS(cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_parent_cell),
                              "Cell can not divide as there are no free neighbours at maximum degree in any direction");
    }

    void TestCellDivisionWithOneNearestNeighbourEmptySiteOnlySearchingNearestNeighbours() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(3, 3, true); // 4*4 nodes

        // Initially we create a cell for every node except node 9
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (i!=9)
            {
                real_node_indices.push_back(i);
            }
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices, true); // Only searching the nearest neighbours (degree=1)

        TS_ASSERT(cell_population.IsEmptySite(9));
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 15u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 15u);

        // Create a new cell
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_new_cell(new Cell(p_state, p_model_2));

        // Add new cell to the cell population by dividing the cell at node 5
        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(5);
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_parent_cell);

        // Check the number of cells is correct
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);

        // Check each cell corresponds to the correct location index
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();

        // The indices of cells 0 to 8 should be unchanged
        for (unsigned i=0; i<9; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // The indices of cells 10 to 15 should be unchanged
        for (unsigned i=10; i<16; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // The new cell should be at node index 9
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 9u);

    }

    void TestAddCellWithNoEmptySites() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(6, 6, true); // 7*7 nodes

        // Initially we create a cell for every node except node 9
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
                real_node_indices.push_back(i);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Create a new cell
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_new_cell(new Cell(p_state, p_model));

        // Try adding new cell to the cell population by dividing the cell at node 24
        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(24);

        TS_ASSERT_THROWS_THIS(cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_parent_cell),
                              "Cell can not divide as there are no free neighbours at maximum degree in any direction");
    }

    void TestIterator()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create a single cell and corresponding location index
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> real_node_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<2> lattice_based_cell_population(mesh, cells, real_node_indices);

        // Count the number of cells using the iterator
        unsigned number_of_cells = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = lattice_based_cell_population.Begin();
             cell_iter != lattice_based_cell_population.End();
             ++cell_iter)
        {
            number_of_cells++;
        }

        TS_ASSERT_EQUALS(number_of_cells, lattice_based_cell_population.rGetCells().size());
    }

    void TestGetFreeNeighbouringNodeIndices1d()
    {
        // Create mesh
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(4); // 5 nodes

        /*
         * Numbering the nodes as follows:
         *
         *  0------1------2------3------4
         */

        // Create one differentiated cell, initially corresponding to the far left node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<1> lattice_based_cell_population(mesh, cells, location_indices);

        // Now find the neighbouring available sites
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 1u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-boundary node 2
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 2);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(2);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test eastern end
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 4);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 1u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeNeighbouringNodeIndices2d()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        /*
         * Numbering the nodes as follows:
         *
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         */

        // Create one cell, initially corresponding to the bottom left node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> cell_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<2> lattice_based_cell_population(mesh, cells, cell_indices);

        // Test bottom left node
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner bottom nodes
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 1);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(0);
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test bottom right node
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 2);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(2);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner left nodes
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 3);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(0);
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(6);
        expected_free_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test centre node
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 4);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 8u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<lattice_based_cell_population.GetNumNodes(); i++)
        {
            if (i != 4)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner right nodes
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 5);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(7);
        expected_free_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test top left node
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 6);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(6);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner top nodes
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 7);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(5);
        expected_free_neighbouring_sites.insert(6);
        expected_free_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test top right node
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 8);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(5);
        expected_free_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Now need to check the case where there is more than one real cell in the cell population

        // Create three cells
        std::vector<CellPtr> cells2;
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<3; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            CellPtr p_cell(new Cell(p_state, p_model));
            cells2.push_back(p_cell);
        }

        // Test non-corner bottom nodes, now with some real nodes thrown in
        std::vector<unsigned> location_indices;
        location_indices.push_back(1);
        location_indices.push_back(4);
        location_indices.push_back(5);

        // Create a cell population
        LatticeBasedCellPopulation<2> lattice_based_cell_population2(mesh, cells2, location_indices);

        // Now find the neighbouring available sites
        free_neighbouring_sites = lattice_based_cell_population2.GetFreeNeighbouringNodeIndices(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(0);
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeNeighbouringNodeIndices3d()
    {
        // Create mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2, 2, 2); // 3*3*3 nodes

        /*
         * Numbering the nodes as follows:
         *
         *  6----7----8       15---16---17      24---25---26
         *  |    |    |       |    |    |       |    |    |
         *  3----4----5       12---13---14      21---22---23
         *  |    |    |       |    |    |       |    |    |
         *  0----1----2       9----10---11      18---19---20
         */

        // Create one cell, initially corresponding to the central node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(13);

        // Create a cell population
        LatticeBasedCellPopulation<3> lattice_based_cell_population(mesh, cells, location_indices);

        // Test central node
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(13);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 26u);

        std::set<unsigned> expected_free_neighbouring_sites;
        for (unsigned i=0; i<27; i++)
        {
            if (i != 13)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test bottom left corner node (node 0)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 0);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 7u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<14; i++)
        {
            if ( i==1 || i==3 || i==4 || i==9 || i==10 || i==12 || i==13 )
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeNeighbouringNodeIndices3dOnNonCubeMesh()
    {
        // Create mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2, 3, 4); // 3*4*5 nodes

        /*
         *  Numbering the nodes as follows:
         *
         *      Bottom layer                layer 2                      layer 3
         *     9---10---11               21---22---23                  33---34---35
         *     |    |    |               |    |    |                   |    |    |
         *     6----7----8               18---19---20                  30---31---32
         *     |    |    |               |    |    |                   |    |    |
         *     3----4----5               15---16---17                  27---28---29
         *     |    |    |               |    |    |                   |    |    |
         *     0----1----2               12---13---14                  24---25---26
         *
         *        layer 4                 Top layer
         *     45---46---47             57---58---59
         *     |     |    |             |     |    |
         *     42---43---44             54---55---56
         *     |     |    |             |     |    |
         *     39---40---41             51---52---53
         *     |     |    |             |     |    |
         *     36---37---38             48---49---50
         */

        // Create one cell, initially corresponding to the corner node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<3> lattice_based_cell_population(mesh, cells, location_indices);

        // Test bottom left corner node (node 0)
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 7u);

        std::set<unsigned> expected_free_neighbouring_sites;
        for (unsigned i=0; i<17; i++)
        {
            if (i==1 || i==3 || i==4 || i==12 || i==13 || i==15 || i==16)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test west side node (node 30)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 30);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(30);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 17u);

        expected_free_neighbouring_sites.clear();

        for (unsigned i=0; i<47; i++)
        {
            if (i==15 || i==16 || i==18 || i==19 || i==21 || i==22 || i==27 || i==28 || i==31 || i==33 || i==34 || i==39 || i==40 || i==42 || i==43 || i==45 || i==46)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test east side node (node 44)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 44);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(44);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 17u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<60; i++)
        {
            if ( i==28 || i==29 || i==31 || i==32 || i==34 || i==35 || i==40 || i==41 || i==43 || i==46 || i==47
              || i==52 || i==53 || i==55 || i==56 || i==58 || i==59 )
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test top layer node (node 52)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 52);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(52);

        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 17u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<57; i++)
        {
            if ( i==36 || i==37 || i==38 || i==39 || i==40 || i==41 || i==42 || i==43 || i==44 || i==48 || i==49
              || i==50 || i==51 || i==53 || i==54 || i==55 || i==56 )
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeNeighbouringNodeIndices2dOnNonSquareMesh()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 3, true); // 3*4 nodes

        /*
         * Numbering the nodes as follows:
         *
         *     9---10---11
         *     |    |    |
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         */

        // Create one cell, initially corresponding to the corner node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<2> lattice_based_cell_population(mesh, cells, location_indices);

        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test east side node (node 5)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 5);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(7);
        expected_free_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-boundary node (node 7)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 7);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 8u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<12; i++)
        {
            if (i > 2 && i != 7)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test north boundary node (node 10)

        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 10);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(10);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();

        for (unsigned i=0; i<12; i++)
        {
            if (i > 5 && i != 10)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestArchivingCellPopulation() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "lattice_based_cell_population.arch";
        ArchiveLocationInfo::SetMeshFilename("lattice_based_cell_population");

        // Archive cell population
        {
            // Set up SimulationTime
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRectangularMesh(6, 6, true); // 7*7 nodes

            // Create location indices
            std::vector<unsigned> real_node_indices;
            for (unsigned i=0; i<10; i++)
            {
                real_node_indices.push_back(2*i);
            }

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);
            for (unsigned i=0; i<cells.size(); i++)
            {
                double birth_time = -2.0 * (double)i;
                cells[i]->SetBirthTime(birth_time);
            }

            // Create cell population
            LatticeBasedCellPopulation<2>* const p_cell_population = new LatticeBasedCellPopulation<2>(mesh, cells, real_node_indices);

            // Cells have been given birth times of 0, -2, -4, ...
            // Loop over them to run to time 0
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                cell_iter != p_cell_population->End();
                ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create an output archive
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

            LatticeBasedCellPopulation<2>* p_cell_population;

            // Restore the cell population
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_cell_population;

            // Check that each cell has the correct age
            unsigned counter = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter = counter + 2;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check there is the correct number of cells
            TS_ASSERT_EQUALS(p_cell_population->GetNumRealCells(), 10u);

            // Check there is the correct number of nodes
            TS_ASSERT_EQUALS(p_cell_population->GetEmptySiteIndices().size(), 39u);

            // Check there is the correct number of nodes
            TS_ASSERT_EQUALS(p_cell_population->GetNumNodes(), 49u);

            // Check some node positions
            Node<2>* p_node_3 = p_cell_population->GetNode(3);
            Node<2>* p_node_11 = p_cell_population->GetNode(11);

            TS_ASSERT_EQUALS(p_node_3->GetIndex(), 3u);
            TS_ASSERT_EQUALS(p_node_11->GetIndex(), 11u);

            TS_ASSERT_DELTA(p_node_3->rGetLocation()[0], 3.0, 1e-9);
            TS_ASSERT_DELTA(p_node_3->rGetLocation()[1], 0.0, 1e-9);
            TS_ASSERT_DELTA(p_node_11->rGetLocation()[0], 4.0, 1e-9);
            TS_ASSERT_DELTA(p_node_11->rGetLocation()[1], 1.0, 1e-9);

            TS_ASSERT_EQUALS(p_node_3->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_node_11->IsBoundaryNode(), false);

            // Tidy up
            delete p_cell_population;
        }
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        // Set up SimulationTime singleton
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(3, 3, true); // 4*4 nodes

        // Create location indices
        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(0);
        real_node_indices.push_back(4);
        real_node_indices.push_back(7);
        real_node_indices.push_back(8);
        real_node_indices.push_back(12);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Set the last cell to start apoptosis
        cells[2]->StartApoptosis();

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 5u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = cell_population.RemoveDeadCells();

        // Test that one cell has been removed
        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u);

        // Test that no nodes have been removed
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);

        // Test that each cell's node index has been correctly updated
        unsigned index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 4*index);
            index++;
        }
    }

    void TestLatticeBasedCellPopulationOutputWriters()
    {
        std::string output_directory = "TestLatticeBasedCellPopulationWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create location indices
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<5; i++)
        {
            real_node_indices.push_back(i);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        boost::shared_ptr<AbstractCellProperty> p_wt(new WildTypeCellMutationState);
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(0);
            cells[i]->SetMutationState(p_wt);
        }

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "LatticeBasedCellPopulation-2");

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVariables(true);

        // For coverage of WriteResultsToFiles()
        boost::shared_ptr<AbstractCellProperty> p_state(cell_population.GetCellPropertyRegistry()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc1(cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(cell_population.GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(cell_population.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(cell_population.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_population.GetCellUsingLocationIndex(0)->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);
        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(1)->GetCellCycleModel()->SetCellProliferativeType(DIFFERENTIATED);
        cell_population.GetCellUsingLocationIndex(1)->SetMutationState(p_apc1);
        cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);
        cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);
        cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_apoptotic_state);

        TS_ASSERT_THROWS_NOTHING(cell_population.CreateOutputFiles(output_directory, false));

        cell_population.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(cell_population.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     notforrelease_cell_based/test/data/TestLatticeBasedCellPopulationWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     notforrelease_cell_based/test/data/TestLatticeBasedCellPopulationWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     notforrelease_cell_based/test/data/TestLatticeBasedCellPopulationWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     notforrelease_cell_based/test/data/TestLatticeBasedCellPopulationWriters/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     notforrelease_cell_based/test/data/TestLatticeBasedCellPopulationWriters/cellages.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 0u);
        for (unsigned i=1; i<4; i++)
        {
            TS_ASSERT_EQUALS(cell_mutation_states[i], 1u);
        }

        // Test the GetCellProliferativeTypeCount function
        std::vector<unsigned> cell_types = cell_population.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 3u);
        TS_ASSERT_EQUALS(cell_types[0], 0u);
        TS_ASSERT_EQUALS(cell_types[1], 1u);
        TS_ASSERT_EQUALS(cell_types[2], 4u);

        // For coverage
        TS_ASSERT_THROWS_NOTHING(cell_population.WriteResultsToFiles());

        //Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");
        // Write cell population parameters to file
        TS_ASSERT_THROWS_THIS(cell_population.OutputCellPopulationParameters(parameter_file),"OutputCellPopulationParameters() is not yet implemented for LatticeBasedCellPopulation see #1453");
        parameter_file->close();

    }

    void TestValidateLatticeBasedCellPopulation()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create location indices
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<5; i++)
        {
            real_node_indices.push_back(i);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Set node 2 to be an empty site, to test Validate()
        cell_population.mIsEmptySite[2] = true;

        TS_ASSERT_THROWS_THIS(cell_population.Validate(), "Node 2 is labelled as an empty site and has a cell attached");

        cell_population.mIsEmptySite[2] = false;

        // Set node 7 to not be an empty site, to test Validate()
        cell_population.mIsEmptySite[7] = false;

        TS_ASSERT_THROWS_THIS(cell_population.Validate(), "Node 7 does not appear to be an empty site or have a cell associated with it");
    }

    void TestOtherMethods() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create two cells
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model_1 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_1->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_cell_1(new Cell(p_state, p_model_1));

        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_cell_2(new Cell(p_state, p_model_2));

        std::vector<CellPtr> cells;
        cells.push_back(p_cell_1);
        cells.push_back(p_cell_2);

        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);
        real_node_indices.push_back(1);

        // Create a cell population
        LatticeBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test UpdateNodeLocations
        std::vector< c_vector<double, 2> > unused_vector;
        TS_ASSERT_THROWS_THIS(cell_population.UpdateNodeLocations(unused_vector, 1.0), "UpdateNodeLocations() cannot be called on a LatticeBasedCellPopulation");

        // Test AddNode()
        Node<2> unused_node(10);
        TS_ASSERT_THROWS_THIS(cell_population.AddNode(&unused_node), "AddNode() cannot be called on a LatticeBasedCellPopulation");

        // Test SetNode()
        ChastePoint<2> unused_point;
        TS_ASSERT_THROWS_THIS(cell_population.SetNode(0, unused_point), "SetNode() cannot be called on a LatticeBasedCellPopulation");

        // Test GetDampingConstant()
        TS_ASSERT_THROWS_THIS(cell_population.GetDampingConstant(0), "GetDampingConstant() cannot be called on a LatticeBasedCellPopulation");

        // Test GetWidth() method
        double width_x = cell_population.GetWidth(0);
        TS_ASSERT_DELTA(width_x, 2.0, 1e-4);

        double width_y = cell_population.GetWidth(1);
        TS_ASSERT_DELTA(width_y, 2.0, 1e-4);
    }

    void TestGetFreeVonNeumannNeighbouringNodeIndices2d()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create one cell, initially corresponding to the bottom left node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> cell_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<2> lattice_based_cell_population(mesh, cells, cell_indices);

        // Set the neighbourhoods to be of von Neumann type
        lattice_based_cell_population.SetVonNeumannNeighbourhoods(true);

        /*
         * Numbering the nodes as follows:
         *
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         */

        // Test bottom left node
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner bottom node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 1);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(0);
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test bottom right node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 2);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(2);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner left nodes
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 3);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(0);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(6);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test centre node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 4);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 4u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(5);
        expected_free_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner right nodes
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 5);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test top left node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 6);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(6);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner top nodes
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 7);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(6);
        expected_free_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test top right node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 8);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(5);
        expected_free_neighbouring_sites.insert(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Now check the case where there is more than one cell in the cell population

        TetrahedralMesh<2,2> mesh2;
        mesh2.ConstructRectangularMesh(2, 2, true); // 3*3 nodes
        std::vector<CellPtr> cells2;
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<3; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            CellPtr p_cell(new Cell(p_state, p_model));
            cells2.push_back(p_cell);
        }

        std::vector<unsigned> location_indices;
        location_indices.push_back(1);
        location_indices.push_back(4);
        location_indices.push_back(5);

        // Create a cell population
        LatticeBasedCellPopulation<2> lattice_based_cell_population2(mesh2, cells2, location_indices);

        // Set the neighbourhoods to be of von Neumann type
        lattice_based_cell_population2.SetVonNeumannNeighbourhoods(true);

        free_neighbouring_sites = lattice_based_cell_population2.GetFreeNeighbouringNodeIndices(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(0);
        expected_free_neighbouring_sites.insert(2);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeVonNeumannNeighbouringNodeIndices3d()
    {
        // Create mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2, 2, 2); // 3*3*3 nodes

        /*
         * Numbering the nodes as follows:
         *
         *  6----7----8       15---16---17      24---25---26
         *  |    |    |       |    |    |       |    |    |
         *  3----4----5       12---13---14      21---22---23
         *  |    |    |       |    |    |       |    |    |
         *  0----1----2       9----10---11      18---19---20
         */

        // Create one cell, initially corresponding to the central node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(13);

        // Create a cell population
        LatticeBasedCellPopulation<3> lattice_based_cell_population(mesh, cells, location_indices);

        // Set von Neumann neighbourhoods
        lattice_based_cell_population.SetVonNeumannNeighbourhoods(true);

        // Test central node
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(13);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 6u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(10);
        expected_free_neighbouring_sites.insert(12);
        expected_free_neighbouring_sites.insert(14);
        expected_free_neighbouring_sites.insert(16);
        expected_free_neighbouring_sites.insert(22);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test bottom left corner node (node 0)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 0);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(9);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeVonNeumannNeighbouringNodeIndices2dOnNonSquareMesh()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 3, true); // 3*4 nodes

        /*
         * Numbering the nodes as follows:
         *
         *     9---10---11
         *     |    |    |
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         */

        // Create one cell, initially corresponding to the corner node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<2> lattice_based_cell_population(mesh, cells, location_indices);

        // Set von Neumann neighbourhoods
        lattice_based_cell_population.SetVonNeumannNeighbourhoods(true);

        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test east side node (node 5)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 5);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-boundary node (node 7)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 7);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 4u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(6);
        expected_free_neighbouring_sites.insert(8);
        expected_free_neighbouring_sites.insert(10);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test north boundary node (node 10)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 10);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(10);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(7);
        expected_free_neighbouring_sites.insert(9);
        expected_free_neighbouring_sites.insert(11);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeVonNeumannNeighbouringNodeIndices3dOnNonCubeMesh()
    {
        // Create mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2, 3, 4); // 3*4*5 nodes

        /*
         *  Numbering the nodes as follows:
         *
         *     Bottom layer                layer 2                       layer 3
         *     9---10---11               21---22---23                 33---34---35
         *     |    |    |               |     |    |                  |    |    |
         *     6----7----8               18---19---20                 30---31---32
         *     |    |    |               |     |    |                  |    |    |
         *     3----4----5               15---16---17                 27---28---29
         *     |    |    |               |     |    |                  |    |    |
         *     0----1----2               12---13---14                 24---25---26
         *
         *       layer 4                  Top layer
         *     45---46---47              57---58---59
         *     |     |    |              |     |    |
         *     42---43---44              54---55---56
         *     |     |    |              |     |    |
         *     39---40---41              51---52---53
         *     |     |    |              |     |    |
         *     36---37---38              48---49---50
         */

        // Create one cell, initially corresponding to the corner node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        LatticeBasedCellPopulation<3> lattice_based_cell_population(mesh, cells, location_indices);

        // Set von Neumann neighbourhoods
        lattice_based_cell_population.SetVonNeumannNeighbourhoods(true);

        // Test bottom left corner node (node 0)
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);

        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        expected_free_neighbouring_sites.insert(12);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test west side node (node 30)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 30);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(30);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(27);
        expected_free_neighbouring_sites.insert(31);
        expected_free_neighbouring_sites.insert(33);
        expected_free_neighbouring_sites.insert(18);
        expected_free_neighbouring_sites.insert(42);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test east side node (node 44)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 44);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(44);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(41);
        expected_free_neighbouring_sites.insert(43);
        expected_free_neighbouring_sites.insert(47);
        expected_free_neighbouring_sites.insert(32);
        expected_free_neighbouring_sites.insert(56);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test top layer node (node 52)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 52);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(52);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(40);
        expected_free_neighbouring_sites.insert(49);
        expected_free_neighbouring_sites.insert(51);
        expected_free_neighbouring_sites.insert(53);
        expected_free_neighbouring_sites.insert(55);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }
};

#endif /*TESTLATTICEBASEDCELLPOPULATION_HPP_*/
