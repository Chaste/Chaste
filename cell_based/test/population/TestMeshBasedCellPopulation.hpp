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

#ifndef TESTMESHBASEDCELLPOPULATION_HPP_
#define TESTMESHBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellAncestor.hpp"
#include "CellId.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

class TestMeshBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    void TestSmallMeshBasedCellPopulation(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(meshFilename);
        MutableMesh<ELEMENT_DIM,SPACE_DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, SPACE_DIM> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the cell population
        unsigned num_cells = cells.size();
        MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), num_cells);

        // Test set/get method of member variables
        TS_ASSERT_DELTA(cell_population.GetMeinekeDivisionSeparation(), 0.3, 1e-6);
        cell_population.SetMeinekeDivisionSeparation(0.5);
        TS_ASSERT_DELTA(cell_population.GetMeinekeDivisionSeparation(), 0.5, 1e-6);
        cell_population.SetMeinekeDivisionSeparation(0.3);

        TS_ASSERT_DELTA(cell_population.GetAbsoluteMovementThreshold(), 2.0, 1e-6);
        cell_population.SetAbsoluteMovementThreshold(1.5);
        TS_ASSERT_DELTA(cell_population.GetAbsoluteMovementThreshold(), 1.5, 1e-6);
        cell_population.SetMeinekeDivisionSeparation(0.5);

        unsigned counter = 0;
        for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            counter++;
        }

        TS_ASSERT_EQUALS(counter, cell_population.GetNumRealCells());
    }

public:

    // Test construction, accessors and Iterator
    void TestSmallMeshBasedCellPopulation1d2d3d() throw(Exception)
    {
        TestSmallMeshBasedCellPopulation<1,1>("mesh/test/data/1D_0_to_1_10_elements");
        TestSmallMeshBasedCellPopulation<2,2>("mesh/test/data/square_4_elements");
        TestSmallMeshBasedCellPopulation<3,3>("mesh/test/data/cube_136_elements");
    }

    // Test construction, accessors and Iterator
    void TestSmallMeshBasedCellPopulation2dIn3d() throw(Exception)
    {
        TestSmallMeshBasedCellPopulation<2,3>("cell_based/test/data/Simple2dMeshIn3d/Simple2dMeshIn3d");
        TestSmallMeshBasedCellPopulation<2,3>("mesh/test/data/disk_in_3d");
    }

    // Test get centroid
    void TestGetCentroidOfCellPopulation() throw(Exception)
    {
        // Create a simple mesh
        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node.
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Check position of centroid
        c_vector<double, 2> expected_centroid_position;
        expected_centroid_position(0) = 0.75;
        expected_centroid_position(1) = 0.25*pow(3,0.5);
        TS_ASSERT_DELTA(cell_population.GetCentroidOfCellPopulation()(0), expected_centroid_position(0), 1e-4);
        TS_ASSERT_DELTA(cell_population.GetCentroidOfCellPopulation()(1), expected_centroid_position(1), 1e-4)
    }

    void TestValidateMeshBasedCellPopulation()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node apart from one.
        // Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Fails as no cell corresponding to node 4
        std::vector<CellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(MeshBasedCellPopulation<2> cell_population2(mesh, cells_copy),
                              "Node 4 does not appear to have a cell associated with it");

        // Add another cell
        AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        double birth_time = -4.0;
        p_cell->SetBirthTime(birth_time);
        cells.push_back(p_cell);

        std::vector<CellPtr> cells_copy2(cells);
        TS_ASSERT_THROWS_NOTHING(MeshBasedCellPopulation<2> cell_population2(mesh, cells_copy2));

        // A bit of Northern compatibility testing hidden here (not relevant to this test!)
        TS_ASSERT_THROWS_NOWT(MeshBasedCellPopulation<2> cell_population2(mesh, cells));
        TS_ASSERT_CHAMPION(true);
    }

    void TestCreateCellPair()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Create two cell pairs
        std::pair<CellPtr,CellPtr> cell_pair1 = cell_population.CreateCellPair(cells[0], cells[1]);
        std::pair<CellPtr,CellPtr> cell_pair2 = cell_population.CreateCellPair(cells[1], cells[0]);
        TS_ASSERT_EQUALS(cell_pair1, cell_pair2);
        TS_ASSERT_EQUALS((cell_pair1.first == cells[0] || cell_pair1.first == cells[1]), true);
        TS_ASSERT_EQUALS((cell_pair1.second == cells[0] || cell_pair1.second == cells[1]), true);
        TS_ASSERT_DIFFERS(cell_pair1.first,  cell_pair1.second);
    }

    void TestGetDampingConstant()
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Bestow mutations on some cells
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApcOneHitCellMutationState, p_apc1);
        MAKE_PTR(ApcTwoHitCellMutationState, p_apc2);
        MAKE_PTR(BetaCateninOneHitCellMutationState, p_bcat1);
        MAKE_PTR(CellLabel, p_label);

        cells[0]->SetMutationState(p_state);
        cells[1]->SetMutationState(p_apc1);
        cells[2]->SetMutationState(p_apc2);
        cells[3]->SetMutationState(p_bcat1);
        cells[4]->AddCellProperty(p_label);

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Change the mutant damping constant to be different from the normal
        cell_population.SetDampingConstantMutant(23.57);

        TS_ASSERT_EQUALS(cell_population.UseAreaBasedDampingConstant(), false);

        double damping_const_0 = cell_population.GetDampingConstant(0);
        double damping_const_1 = cell_population.GetDampingConstant(1);
        double damping_const_2 = cell_population.GetDampingConstant(2);
        double damping_const_3 = cell_population.GetDampingConstant(3);
        double damping_const_4 = cell_population.GetDampingConstant(4);

        // Check that each mutation state gives the correct damping constant
        TS_ASSERT_DELTA(damping_const_0, cell_population.GetDampingConstantNormal(), 1e-6);
        TS_ASSERT_DELTA(damping_const_1, cell_population.GetDampingConstantMutant(), 1e-6);
        TS_ASSERT_DELTA(damping_const_2, cell_population.GetDampingConstantMutant(), 1e-6);
        TS_ASSERT_DELTA(damping_const_3, cell_population.GetDampingConstantMutant(), 1e-6);
        TS_ASSERT_DELTA(damping_const_4, cell_population.GetDampingConstantMutant(), 1e-6);

        // Coverage
        TS_ASSERT_DELTA(cell_population.GetAreaBasedDampingConstantParameter(), 0.1, 1e-6);
        cell_population.SetAreaBasedDampingConstantParameter(0.5);
        TS_ASSERT_DELTA(cell_population.GetAreaBasedDampingConstantParameter(), 0.5, 1e-6);

        //test Get and Set methods for DampingConstantNormal and DampingConstantNormalMutant
        TS_ASSERT_DELTA(cell_population.GetDampingConstantNormal(),1.0, 1e-6);
        TS_ASSERT_DELTA(cell_population.GetDampingConstantMutant(),23.57, 1e-6);

        cell_population.SetDampingConstantNormal(2.0);
        cell_population.SetDampingConstantMutant(3.0);

        TS_ASSERT_DELTA(cell_population.GetDampingConstantNormal(),2.0, 1e-6);
        TS_ASSERT_DELTA(cell_population.GetDampingConstantMutant(),3.0, 1e-6);
    }

    void TestAreaBasedDampingConstant()
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        MAKE_PTR(ApcTwoHitCellMutationState, p_apc2);
        cells[9]->SetMutationState(p_apc2);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.UseAreaBasedDampingConstant(), false);

        double damping_const = cell_population.GetDampingConstant(8);

        TS_ASSERT_DELTA(damping_const, cell_population.GetDampingConstantNormal(), 1e-6);

        double mutant_damping_const = cell_population.GetDampingConstant(9);

        TS_ASSERT_DELTA(mutant_damping_const, cell_population.GetDampingConstantMutant(), 1e-6);

        cell_population.SetAreaBasedDampingConstant(true);

        TS_ASSERT_EQUALS(cell_population.UseAreaBasedDampingConstant(), true);

        // Note that this method is usually called by OffLatticeSimulation::Solve()
        cell_population.CreateVoronoiTessellation();

        double area_based_damping_const = cell_population.GetDampingConstant(8);

        // Since the cell population is in mechanical equilibrium, we should get the same damping constant as before
        TS_ASSERT_DELTA(area_based_damping_const, cell_population.GetDampingConstantNormal(), 1e-6);
    }

    void TestSetNodeAndAddCell()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population, with no ghost nodes at the moment
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Test SetNode() by moving node 0 by a small amount
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        c_vector<double,2> new_location = cell_population.GetLocationOfCellCentre(*cell_iter);
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        cell_population.SetNode(cell_population.GetLocationIndexUsingCell(*cell_iter), new_location_point);

        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        // Test AddCell()
        unsigned old_num_nodes = mesh.GetNumNodes();
        unsigned old_num_cells = cell_population.rGetCells().size();

        // Create a new cell, DON'T set the node index, set birth time=-1
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->SetBirthTime(-1);
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2;
        new_cell_location[1] = 2;

        cell_population.AddCell(p_cell, new_cell_location, cells[3] /*random choice of parent*/);

        // CellPopulation should have updated mesh and cells
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), old_num_nodes+1);

        // Same test via cell population class
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), old_num_cells+1);

        // Check the location of the new node
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // Check the index of the new cell
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(cell_population.rGetCells().back()), old_num_nodes);
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        cells[27]->StartApoptosis();

        // Create a cell population without ghost nodes
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 81u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 81u);

        // Test GetNeighbouringNodeIndices() method
        std::set<unsigned> node_50_neighbours = cell_population.GetNeighbouringNodeIndices(50);

        std::set<unsigned> expected_node_50_neighbours;
        expected_node_50_neighbours.insert(10);
        expected_node_50_neighbours.insert(18);
        expected_node_50_neighbours.insert(27);
        expected_node_50_neighbours.insert(34);

        TS_ASSERT_EQUALS(node_50_neighbours.size(), expected_node_50_neighbours.size());
        TS_ASSERT_EQUALS(node_50_neighbours, expected_node_50_neighbours);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = cell_population.RemoveDeadCells();

        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 80u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 80u);
        TS_ASSERT_DIFFERS(cell_population.rGetCells().size(), cells.size()); // CellPopulation now copies cells

        cell_population.Update();

        // For coverage
        NodeMap map(mesh.GetNumAllNodes());
        map.ResetToIdentity();
        cell_population.UpdateGhostNodesAfterReMesh(map);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 80u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumAllNodes());

        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            TS_ASSERT_EQUALS(cell_population.IsGhostNode(i), false);
        }

        // Finally, check the cells node indices have updated

        // We expect the cell node indices to be {0,11,...,79}
        std::set<unsigned> expected_node_indices;
        for (unsigned i=0; i<cell_population.GetNumRealCells(); i++)
        {
            expected_node_indices.insert(i);
        }

        // Get actual cell node indices
        std::set<unsigned> node_indices;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Record node index corresponding to cell
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            node_indices.insert(node_index);
        }

        TS_ASSERT_EQUALS(node_indices, expected_node_indices);
    }

    void TestUpdateNodeLocations()
    {
        // Test MeshBasedCellPopulation::UpdateNodeLocations()

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population, with no ghost nodes at the moment
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Make up some forces
        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
        std::vector<c_vector<double, 2> > forces_on_nodes(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            old_posns[i][0] = cell_population.GetNode(i)->rGetLocation()[0];
            old_posns[i][1] = cell_population.GetNode(i)->rGetLocation()[1];

            forces_on_nodes[i][0] = i*0.01;
            forces_on_nodes[i][1] = 2*i*0.01;
        }

        // Call method
        double time_step = 0.01;
        cell_population.UpdateNodeLocations(forces_on_nodes, time_step);

        // Check that node locations were correctly updated
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
        }
    }

    void TestVoronoiMethods()
    {
        // First test the 2D case...

        // Create 2D mesh
        std::vector<Node<2> *> nodes2d;
        nodes2d.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes2d.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes2d.push_back(new Node<2>(3, true, 0.0, 1.0));
        nodes2d.push_back(new Node<2>(4, false, 0.5, 0.5));
        MutableMesh<2,2> mesh2d(nodes2d);

        // Create cells
        std::vector<CellPtr> cells2d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2d;
        cells_generator2d.GenerateBasic(cells2d, mesh2d.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population2d(mesh2d, cells2d);

        // For coverage of GetVolumeOfCell()
        for (AbstractCellPopulation<2>::Iterator iter = cell_population2d.Begin();
             iter != cell_population2d.End();
             ++iter)
        {
            if (cell_population2d.GetLocationIndexUsingCell(*iter) == 4)
            {
                TS_ASSERT_DELTA(cell_population2d.GetVolumeOfCell(*iter), 0.5, 1e-6);
            }
            else
            {
                TS_ASSERT_DELTA(cell_population2d.GetVolumeOfCell(*iter), 0.0, 1e-6);
            }
        }

        // The Voronoi tessellation has been created by this point, since GetVolumeOfCell() has been called

        // Test element areas
        TS_ASSERT_DELTA(cell_population2d.GetVolumeOfVoronoiElement(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetVolumeOfVoronoiElement(1), 0.0, 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetVolumeOfVoronoiElement(2), 0.0, 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetVolumeOfVoronoiElement(3), 0.0, 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetVolumeOfVoronoiElement(4), 0.5, 1e-6);

        // Test element perimeters
        TS_ASSERT_DELTA(cell_population2d.GetSurfaceAreaOfVoronoiElement(0), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetSurfaceAreaOfVoronoiElement(1), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetSurfaceAreaOfVoronoiElement(2), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetSurfaceAreaOfVoronoiElement(3), sqrt(2), 1e-6);
        TS_ASSERT_DELTA(cell_population2d.GetSurfaceAreaOfVoronoiElement(4), 2*sqrt(2), 1e-6);

        // ...now test the 3D case

        // Create 3D mesh
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes3d.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes3d.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh3d(nodes3d);

        // Create cells
        std::vector<CellPtr> cells3d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator3d;
        cells_generator3d.GenerateBasic(cells3d, mesh3d.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<3> cell_population3d(mesh3d, cells3d);

        // Create Voronoi tessellation
        cell_population3d.CreateVoronoiTessellation();

        // Test element volumes and surface areas
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_THROWS_THIS(cell_population3d.GetVolumeOfVoronoiElement(i),
                                  "This index does not correspond to a VertexElement");

            TS_ASSERT_THROWS_THIS(cell_population3d.GetSurfaceAreaOfVoronoiElement(i),
                                  "This index does not correspond to a VertexElement");
        }

        // The Voronoi tessellation should comprise a single tetrahedral VertexElement
        TS_ASSERT_EQUALS(cell_population3d.GetVoronoiTessellation()->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(cell_population3d.GetVoronoiTessellation()->GetNumFaces(), 4u);
        TS_ASSERT_EQUALS(cell_population3d.GetVoronoiTessellation()->GetNumElements(), 1u);

        // The faces are not all equal
        for (unsigned face_index=0; face_index<4; face_index++)
        {
            VertexElement<2,3>* p_face = cell_population3d.GetVoronoiTessellation()->GetFace(face_index);

            if (face_index == 1)
            {
                TS_ASSERT_DELTA(cell_population3d.GetVoronoiTessellation()->GetAreaOfFace(p_face), 1.9485, 1e-4);
            }
            else
            {
                TS_ASSERT_DELTA(cell_population3d.GetVoronoiTessellation()->GetAreaOfFace(p_face), 1.125, 1e-4);
            }
        }

        TS_ASSERT_DELTA(cell_population3d.GetVolumeOfVoronoiElement(4), 0.6495, 1e-4);
        TS_ASSERT_DELTA(cell_population3d.GetSurfaceAreaOfVoronoiElement(4), 3*1.125 + 1.9485, 1e-4);

        // Check that the Voronoi tessellation can be returned successfully as a reference
        VertexMesh<3,3>* p_tessellation1 = cell_population3d.GetVoronoiTessellation();
        TS_ASSERT_EQUALS(p_tessellation1->GetNumNodes(), 4u);

        // Move node 0 by a small amount
        AbstractCellPopulation<3>::Iterator cell_iter = cell_population3d.Begin();
        c_vector<double,3> new_location = cell_population3d.GetLocationOfCellCentre(*cell_iter);
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        new_location[2] += 1e-2;
        ChastePoint<3> new_location_point(new_location);
        cell_population3d.SetNode(cell_population3d.GetLocationIndexUsingCell(*cell_iter), new_location_point);

        // Re-create Voronoi tessellation
        cell_population3d.CreateVoronoiTessellation();

        VertexMesh<3,3>* p_tessellation2 = cell_population3d.GetVoronoiTessellation();
        TS_ASSERT_EQUALS(p_tessellation2->GetNumNodes(), 4u);

//        TS_ASSERT_EQUALS(p_tessellation1->GetNumNodes(), 4u);
    }

    void TestCellPopulationWritersIn2d()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Resetting the Maximum cell Id to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Cover mutation state reporting
        boost::shared_ptr<AbstractCellProperty> p_apc1(CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(CellPropertyRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        cells[0]->AddCellProperty(p_apoptotic_state);
        cells[1]->SetMutationState(p_apc1);
        cells[2]->SetMutationState(p_apc2);
        cells[3]->SetMutationState(p_bcat1);
        cells[4]->AddCellProperty(p_label);

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "MeshBasedCellPopulation-2-2");

        // Test set/get methods
        TS_ASSERT_EQUALS(cell_population.GetOutputVoronoiData(), false);
        TS_ASSERT_EQUALS(cell_population.GetOutputCellIdData(), false);
        TS_ASSERT_EQUALS(cell_population.GetWriteVtkAsPoints(), false);

        cell_population.SetOutputVoronoiData(true);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetWriteVtkAsPoints(true);

        TS_ASSERT_EQUALS(cell_population.GetOutputVoronoiData(), true);
        TS_ASSERT_EQUALS(cell_population.GetOutputCellIdData(), true);
        TS_ASSERT_EQUALS(cell_population.GetWriteVtkAsPoints(), true);

        // Coverage of writing CellData to VTK
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("var1", (double) 0.0);
            cell_iter->GetCellData()->SetItem("var2", (double) 3.0);
        }

        // Test set methods
        cell_population.SetOutputResultsForChasteVisualizer(true);
        cell_population.SetOutputCellPopulationVolumes(true);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellAncestors(true);

        // This method is usually called by Update()
        cell_population.CreateVoronoiTessellation();

        std::string output_directory = "TestCellPopulationWritersIn2d";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.CreateOutputFiles(output_directory, false);
        cell_population.WriteResultsToFiles();
        cell_population.CloseOutputFiles();

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 2u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 1u);

        // Test the GetCellProliferativeTypeCount() function - we should have five stem cells
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

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("results.vizelements");
        files_to_compare.push_back("results.viznodes");
        files_to_compare.push_back("results.vizcelltypes");
        files_to_compare.push_back("cellpopulationareas.dat");
        files_to_compare.push_back("cellareas.dat");
        files_to_compare.push_back("cellmutationstates.dat");
        files_to_compare.push_back("voronoi.dat");
        files_to_compare.push_back("results.vizancestors");
        files_to_compare.push_back("results.parameters");

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            FileComparison comparer(results_dir + files_to_compare[i],"cell_based/test/data/TestCellPopulationWritersIn2d/" + files_to_compare[i]);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestCellPopulationWritersIn3d()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Resetting the Maximum cell Id to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 3D mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        cells[4]->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>()); // coverage

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "MeshBasedCellPopulation-3-3");

        // Coverage of writing CellData to VTK
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("a variable", 0.0);
            cell_iter->GetCellData()->SetItem("another", 100.0);
        }

        // Test set methods
        cell_population.SetOutputVoronoiData(true);
        cell_population.SetOutputCellPopulationVolumes(true);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellCyclePhases(true);

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellAncestors(true);

        // This method is usually called by Update()
        cell_population.CreateVoronoiTessellation();

        std::string output_directory = "TestCellPopulationWritersIn3d";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.CreateOutputFiles(output_directory, false);
        cell_population.WriteResultsToFiles();
        cell_population.CloseOutputFiles();

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 5u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 0u);

        // Test the GetCellProliferativeTypeCount function
        std::vector<unsigned> cell_types = cell_population.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 4u);
        TS_ASSERT_EQUALS(cell_types[0], 5u);
        TS_ASSERT_EQUALS(cell_types[1], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("results.vizelements");
        files_to_compare.push_back("results.viznodes");
        files_to_compare.push_back("results.vizcelltypes");
        files_to_compare.push_back("cellpopulationareas.dat");
        files_to_compare.push_back("cellareas.dat");
        files_to_compare.push_back("cellmutationstates.dat");
        files_to_compare.push_back("voronoi.dat");
        files_to_compare.push_back("results.vizancestors");

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            FileComparison comparer(results_dir + files_to_compare[i],"cell_based/test/data/TestCellPopulationWritersIn3d/" + files_to_compare[i]);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestGetLocationOfCellCentreAndGetNodeCorrespondingToCellAndGetWidth() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Loop over nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Record node location
            c_vector<double,2> node_location = cell_population.GetLocationOfCellCentre(*cell_iter);

            // Test GetLocationOfCellCentre()
            TS_ASSERT_DELTA(node_location[0], cell_population.GetLocationOfCellCentre(*cell_iter)[0], 1e-9);
            TS_ASSERT_DELTA(node_location[1], cell_population.GetLocationOfCellCentre(*cell_iter)[1], 1e-9);
        }

        // Test GetWidth() method
        double width_x = cell_population.GetWidth(0);
        TS_ASSERT_DELTA(width_x, 1.0, 1e-6);

        double width_y = cell_population.GetWidth(1);
        TS_ASSERT_DELTA(width_y, 1.0, 1e-6);
    }

    void TestAddCellDataToPopulation()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index, so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        cell_population.SetDataOnAllCells("variable", 100.0);

        //Check that the data made it there and that copies of the data are independent
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

    // This test checks that the cells and nodes are correctly archived.
    void TestArchivingMeshBasedCellPopulation() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mesh_based_cell_population.arch";
        ArchiveLocationInfo::SetMeshFilename("mesh_based_cell_population_mesh");

        std::vector<c_vector<double,2> > cell_locations;

        // Archive a cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells, one for each node. Give each a birth time of -node_index,
            // so the age = node_index
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            // Create the cell population
            MeshBasedCellPopulation<2>* const p_cell_population = new MeshBasedCellPopulation<2>(mesh, cells);

            //Add cell data
            p_cell_population->SetDataOnAllCells("nothing", DOUBLE_UNSET);

            TS_ASSERT_THROWS_THIS(p_cell_population->Begin()->GetCellData()->GetItem("nothing"),
                                "The item nothing has not yet been set");

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            unsigned index_for_data = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter=p_cell_population->Begin();
                 cell_iter!=p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
                cell_locations.push_back(p_cell_population->GetLocationOfCellCentre(*cell_iter));
                cell_iter->GetCellData()->SetItem("data", (double) index_for_data);
                index_for_data++;
            }

            std::pair<CellPtr,CellPtr> cell_pair_0_1 = p_cell_population->CreateCellPair(p_cell_population->GetCellUsingLocationIndex(0), p_cell_population->GetCellUsingLocationIndex(1));
            p_cell_population->MarkSpring(cell_pair_0_1);

            // Set area-based viscosity
            p_cell_population->SetAreaBasedDampingConstant(true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&> (*p_simulation_time);
            (*p_arch) << p_cell_population;
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        // Restore cell population
        {
            // Need to set up time
            unsigned num_steps=10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            MeshBasedCellPopulation<2>* p_cell_population;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> *p_simulation_time;

            (*p_arch) >> p_cell_population;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter=p_cell_population->Begin();
                 cell_iter!=p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(),(double)(counter),1e-7);
                TS_ASSERT_DELTA(p_cell_population->GetLocationOfCellCentre(*cell_iter)[0], cell_locations[counter][0], 1e-9);
                TS_ASSERT_DELTA(p_cell_population->GetLocationOfCellCentre(*cell_iter)[1], cell_locations[counter][1], 1e-9);
                TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("data"), (double) p_cell_population->GetLocationIndexUsingCell(*cell_iter), 1e-12);
                counter++;
            }

            TS_ASSERT_EQUALS(p_cell_population->GetNode(0)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(1)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(2)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(3)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(4)->IsBoundaryNode(), false);

            // Check the marked spring
            std::pair<CellPtr,CellPtr> cell_pair_0_1 = p_cell_population->CreateCellPair(p_cell_population->GetCellUsingLocationIndex(0), p_cell_population->GetCellUsingLocationIndex(1));
            TS_ASSERT_EQUALS(p_cell_population->IsMarkedSpring(cell_pair_0_1), true);

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the cell population has been restored
            TS_ASSERT_EQUALS(p_cell_population->rGetCells().size(), 5u);

            // Check area-based viscosity is still true
            TS_ASSERT_EQUALS(p_cell_population->UseAreaBasedDampingConstant(), true);

            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumNodes(), 5u);

            delete p_cell_population;
        }
    }

    void TestSpringMarking()
    {
        // Create a small cell population
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));

        MutableMesh<2,2> mesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        std::pair<CellPtr,CellPtr> cell_pair_1_2 = cell_population.CreateCellPair(cell_population.GetCellUsingLocationIndex(1), cell_population.GetCellUsingLocationIndex(2));
        std::pair<CellPtr,CellPtr> cell_pair_3_4 = cell_population.CreateCellPair(cell_population.GetCellUsingLocationIndex(3), cell_population.GetCellUsingLocationIndex(4));
        std::pair<CellPtr,CellPtr> cell_pair_1_4 = cell_population.CreateCellPair(cell_population.GetCellUsingLocationIndex(1), cell_population.GetCellUsingLocationIndex(4));
        std::pair<CellPtr,CellPtr> cell_pair_0_2 = cell_population.CreateCellPair(cell_population.GetCellUsingLocationIndex(0), cell_population.GetCellUsingLocationIndex(2));

        // Mark some springs
        cell_population.MarkSpring(cell_pair_1_2);

        // Unmark and re-mark spring (for coverage)
        cell_population.UnmarkSpring(cell_pair_1_2);
        cell_population.MarkSpring(cell_pair_1_2);

        cell_population.MarkSpring(cell_pair_3_4);

        // Check if springs are marked
        TS_ASSERT_EQUALS(cell_population.IsMarkedSpring(cell_pair_1_2), true);
        TS_ASSERT_EQUALS(cell_population.IsMarkedSpring(cell_pair_3_4), true);

        TS_ASSERT_EQUALS(cell_population.IsMarkedSpring(cell_pair_1_4), false);
        TS_ASSERT_EQUALS(cell_population.IsMarkedSpring(cell_pair_0_2), false);

        // Delete cell 4
        cell_population.GetCellUsingLocationIndex(4)->Kill();
        cell_population.RemoveDeadCells();

        // Check springs with non-deleted cells are still marked
        TS_ASSERT_EQUALS(cell_population.IsMarkedSpring(cell_pair_1_2), true);
        cell_population.CheckCellPointers();

        // Move cell 2
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        ++cell_iter;
        ++cell_iter;
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 2u);
        ChastePoint<2> new_location(1, 10);
        cell_population.SetNode(cell_population.GetLocationIndexUsingCell(*cell_iter), new_location);

        // Update cell population
        cell_population.Update();

        cell_population.CheckCellPointers();

        // Check there is no marked spring between nodes 1 & 2
        TS_ASSERT_EQUALS(cell_population.IsMarkedSpring(cell_pair_1_2), false);
    }

    void TestSettingCellAncestors() throw (Exception)
    {
        // Create a small mesh-based cell population
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));
        MutableMesh<2,2> mesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Test that the cell population makes each cell fix the corresponding node index as its ancestor
        cell_population.SetCellAncestorsToLocationIndices();

        unsigned counter = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetAncestor(), cell_population.GetLocationIndexUsingCell(*cell_iter));
            counter++;
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

    void TestIsCellAssociatedWithADeletedLocation() throw (Exception)
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(4, 4, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->GetNode(0)->MarkAsDeleted();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create cell population but do not try to validate
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, std::vector<unsigned>(), false, false);

        // Test IsCellAssociatedWithADeletedLocation() method
        for (MeshBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            bool is_deleted = cell_population.IsCellAssociatedWithADeletedLocation(*cell_iter);

            if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 0)
            {
                TS_ASSERT_EQUALS(is_deleted, true);
            }
            else
            {
                TS_ASSERT_EQUALS(is_deleted, false);
            }
        }
    }
};

#endif /*TESTMESHBASEDCELLPOPULATION_HPP_*/
