/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTNODEBASEDCELLPOPULATIONWITHPARTICLES_HPP_
#define TESTNODEBASEDCELLPOPULATIONWITHPARTICLES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellwiseData.hpp"
#include "SmartPointers.hpp"

class TestNodeBasedCellPopulationWithParticles : public AbstractCellBasedTestSuite
{
public:
	 /*
	 * Here we set up a test with 5 nodes, make a cell for each. We then set cell
	 * 0 to be associated with node 1 instead of node 0, and Validate() throws an
	 * exception. We then set node 0 to be a particle node, and Validate() passes.
	 */
	void TestValidateNodeBasedCellPopulationWithParticles() throw(Exception)
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

		std::vector<unsigned> cell_location_indices;
		for (unsigned i=0; i<cells.size(); i++)
		{
			cell_location_indices.push_back(i);
		}

		// Fails as the cell population constructor is not given the location indices
		// corresponding to real cells, so cannot work out which nodes are
		// particles
		std::vector<CellPtr> cells_copy(cells);
		TS_ASSERT_THROWS_THIS(NodeBasedCellPopulationWithParticles<2> dodgy_cell_population(mesh, cells_copy),
				"Node 4 does not appear to be a particle or has a cell associated with it");

		// Passes as the cell population constructor automatically works out which
		// cells are particles using the mesh and cell_location_indices
		NodeBasedCellPopulationWithParticles<2> cell_population(mesh, cells, cell_location_indices);

		// Here we set the particles to what they already are
		std::set<unsigned> particle_indices;
		particle_indices.insert(mesh.GetNumNodes()-1u);
		cell_population.SetParticles(particle_indices);

		// So validate passes at the moment
		cell_population.Validate();

		// Test GetCellUsingLocationIndex()

		TS_ASSERT_THROWS_NOTHING(cell_population.GetCellUsingLocationIndex(0)); // real cell
		TS_ASSERT_THROWS_THIS(cell_population.GetCellUsingLocationIndex(mesh.GetNumNodes()-1u),"Location index input argument does not correspond to a Cell"); // particles

		// Now we label a real cell's node as particle
		particle_indices.insert(1u);

		// Validate detects this inconsistency
		TS_ASSERT_THROWS_THIS(cell_population.SetParticles(particle_indices),"Node 1 is labelled as a particle and has a cell attached");
	}

	// Test with particles, checking that the Iterator doesn't loop over particles
	void TestNodeBasedCellPopulationWithParticlesSetup() throw(Exception)
	{
		unsigned num_cells_depth = 11;
		unsigned num_cells_width = 6;
		HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		// Set up cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel,2> cells_generator;
		cells_generator.GenerateGivenLocationIndices(cells, location_indices);

		// Create a cell population
		NodeBasedCellPopulationWithParticles<2> cell_population(mesh, cells, location_indices);

		// Create a set of node indices corresponding to particles
		std::set<unsigned> node_indices;
		std::set<unsigned> location_indices_set;
		std::set<unsigned> particle_indices;

		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			node_indices.insert(mesh.GetNode(i)->GetIndex());
		}
		for (unsigned i=0; i<location_indices.size(); i++)
		{
			location_indices_set.insert(location_indices[i]);
		}

		std::set_difference(node_indices.begin(), node_indices.end(),
							location_indices_set.begin(), location_indices_set.end(),
							std::inserter(particle_indices, particle_indices.begin()));

		std::vector<bool> is_particle(mesh.GetNumNodes(), false);
		for (std::set<unsigned>::iterator it=particle_indices.begin();
			 it!=particle_indices.end();
			 it++)
		{
			is_particle[*it] = true;
		}

		TS_ASSERT_EQUALS(cell_population.rGetParticles(), is_particle);

		// Test the GetParticleIndices method
		std::set<unsigned> particle_indices2 = cell_population.GetParticleIndices();
		TS_ASSERT_EQUALS(particle_indices, particle_indices2);

		// Check the iterator doesn't loop over particles
		unsigned counter = 0;
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			 cell_iter != cell_population.End();
			 ++cell_iter)
		{
			unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
			TS_ASSERT_EQUALS(is_particle[node_index], false);
			counter++;
		}

		TS_ASSERT_EQUALS(counter, cell_population.GetNumRealCells());

		// Check counter = num_nodes - num_particles_nodes
		TS_ASSERT_EQUALS(counter + particle_indices.size(), mesh.GetNumNodes());

		TS_ASSERT_EQUALS(cell_population.rGetParticles().size(), mesh.GetNumNodes());
	}

	void TestCellPopulationIteratorWithNoCells() throw(Exception)
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

	        // Create vector of cell location indices
	        std::vector<unsigned> cell_location_indices;
	        cell_location_indices.push_back(80);

	        // Create a single cell
	        std::vector<CellPtr> cells;
	        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
	        cells_generator.GenerateBasic(cells, cell_location_indices.size());
	        cells[0]->StartApoptosis();

	        // Create a cell population
	        NodeBasedCellPopulationWithParticles<2> cell_population(mesh, cells, cell_location_indices);
	        cell_population.SetMechanicsCutOffLength(1.2);

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

	void TestRemoveDeadCellsAndReMeshWithParticles()
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
	        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
	        cells_generator.GenerateBasic(cells, cell_location_indices.size());
	        cells[27]->StartApoptosis();

	        // Create a cell population, with some random particles
	        NodeBasedCellPopulationWithParticles<2> cell_population_with_particles(mesh, cells, cell_location_indices);
	        cell_population_with_particles.SetMechanicsCutOffLength(1.2);


	        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);

	        // Num real cells should be num_nodes (81) - num_particles (11) = 70
	        TS_ASSERT_EQUALS(cell_population_with_particles.GetNumRealCells(), 70u);

	        p_simulation_time->IncrementTimeOneStep();

	        unsigned num_removed_with_particles = cell_population_with_particles.RemoveDeadCells();

	        TS_ASSERT_EQUALS(num_removed_with_particles, 1u);
	        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
	        TS_ASSERT_DIFFERS(cell_population_with_particles.rGetCells().size(), cells.size()); // CellPopulation now copies cells

	        // Num real cells should be num_nodes (81) - num_particle (11) - 1 deleted node = 69
	        TS_ASSERT_EQUALS(cell_population_with_particles.GetNumRealCells(), 69u);
	        TS_ASSERT_EQUALS(cell_population_with_particles.rGetParticles().size(), mesh.GetNumAllNodes());

	        cell_population_with_particles.Update();

	        // For coverage
	        NodeMap map(mesh.GetNumAllNodes());
	        map.ResetToIdentity();
	        cell_population_with_particles.UpdateParticlesAfterReMesh(map);

	        // Num real cells should be new_num_nodes (80) - num_particles (11)
	        TS_ASSERT_EQUALS(cell_population_with_particles.GetNumRealCells(), 69u);
	        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumAllNodes());
	        TS_ASSERT_EQUALS(cell_population_with_particles.rGetParticles().size(), mesh.GetNumNodes());

	        // Nodes 0-9 should not been renumbered so are still particles.
	        // the particle at node 80 is now at 79 as node 27 was deleted..
	        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
	        {
	            // True (ie should be a particle) if i<10 or i==79, else false
	            TS_ASSERT_EQUALS(cell_population_with_particles.IsParticle(i), ((i<10)||(i==79)));
	        }

	        // Finally, check the cells node indices have updated

	        // We expect the cell node indices to be {10,11,...,79}
	        std::set<unsigned> expected_node_indices;
	        for (unsigned i=0; i<cell_population_with_particles.GetNumRealCells(); i++)
	        {
	            expected_node_indices.insert(i+10);
	        }

	        // Get actual cell node indices
	        std::set<unsigned> node_indices_with_particles;

	        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population_with_particles.Begin();
	             cell_iter != cell_population_with_particles.End();
	             ++cell_iter)
	        {
	            // Record node index corresponding to cell
	            unsigned node_index_with_particles = cell_population_with_particles.GetLocationIndexUsingCell(*cell_iter);
	            node_indices_with_particles.insert(node_index_with_particles);
	        }

	        TS_ASSERT_EQUALS(node_indices_with_particles, expected_node_indices);
	    }

	void TestAddAndRemoveAndAddWithOutUpdate()
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
	        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
	        cells_generator.GenerateBasic(cells, cell_location_indices.size());
	        cells[27]->StartApoptosis();

	        // Create a cell population, with some random particles
	        NodeBasedCellPopulationWithParticles<2> cell_population(mesh, cells, cell_location_indices);
	        cell_population.SetMechanicsCutOffLength(1.2);

	        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
	        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 70u);

	        MAKE_PTR(WildTypeCellMutationState, p_state);

	        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
	        p_model->SetCellProliferativeType(STEM);
	        CellPtr p_new_cell(new Cell(p_state, p_model));
	        p_new_cell->SetBirthTime(0);

	        c_vector<double,2> new_location;
	        new_location[0] = 0.3433453454443;
	        new_location[0] = 0.3435346344234;
	        cell_population.AddCell(p_new_cell, new_location, cells[0] /*random choice of parent*/);

	        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
	        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 71u);

	        p_simulation_time->IncrementTimeOneStep();

	        unsigned num_removed = cell_population.RemoveDeadCells();
	        TS_ASSERT_EQUALS(num_removed, 1u);

	        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
	        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 70u);

	        FixedDurationGenerationBasedCellCycleModel* p_model2 = new FixedDurationGenerationBasedCellCycleModel();
	        p_model2->SetCellProliferativeType(STEM);
	        CellPtr p_new_cell2(new Cell(p_state, p_model2));
	        p_new_cell2->SetBirthTime(0);

	        c_vector<double,2> new_location2;
	        new_location2[0] = 0.6433453454443;
	        new_location2[0] = 0.6435346344234;
	        cell_population.AddCell(p_new_cell2, new_location2, cells[1] /*random choice of parent*/);

	        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
	        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 71u);
	    }

	void TestUpdateNodeLocations() throw(Exception)
	    {
			HoneycombMeshGenerator generator(3, 3, 1);
			TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

			// Convert this to a NodesOnlyMesh
			MAKE_PTR(NodesOnlyMesh<2>, p_mesh);
			p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

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
	            CellProliferativeType cell_type;
	            unsigned generation;
	            double y = 0.0;

	            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
	            {
	                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
	            }

	            FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel;
	            p_cell_cycle_model->SetDimension(2);

	            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
	            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

	            double birth_time = 0.0;
	            birth_time = -p_random_num_gen->ranf();

	            if (y <= 0.3)
	            {
	                cell_type = STEM;
	                generation = 0;
	                birth_time *= typical_stem_cycle_time; // hours
	            }
	            else if (y < 2.0)
	            {
	                cell_type = TRANSIT;
	                generation = 1;
	                birth_time *= typical_transit_cycle_time; // hours
	            }
	            else if (y < 3.0)
	            {
	                cell_type = TRANSIT;
	                generation = 2;
	                birth_time *= typical_transit_cycle_time; // hours
	            }
	            else if (y < 4.0)
	            {
	                cell_type = TRANSIT;
	                generation = 3;
	                birth_time *= typical_transit_cycle_time; // hours
	            }
	            else
	            {
	                cell_type = p_cell_cycle_model->CanCellTerminallyDifferentiate() ? DIFFERENTIATED : TRANSIT;
	                generation = 4;
	                birth_time *= typical_transit_cycle_time; // hours
	            }

	            p_cell_cycle_model->SetGeneration(generation);
	            p_cell_cycle_model->SetCellProliferativeType(cell_type);

	            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

	            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
	            p_cell->SetBirthTime(birth_time);

	            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
	            {
	                cells.push_back(p_cell);
	            }
	        }

	        NodeBasedCellPopulationWithParticles<2> cell_population(*p_mesh, cells, location_indices);
	        cell_population.SetMechanicsCutOffLength(2.0);

	        // Make up some forces
	        std::vector<c_vector<double, 2> > old_posns(cell_population.GetNumNodes());
	        std::vector<c_vector<double, 2> > forces(cell_population.GetNumNodes());
	        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
	        {
	            old_posns[i][0] = cell_population.GetNode(i)->rGetLocation()[0];
	            old_posns[i][1] = cell_population.GetNode(i)->rGetLocation()[1];

	            forces[i][0] = i*0.01;
	            forces[i][1] = 2*i*0.01;
	        }

	        // Call method
	        double time_step = 0.01;
	        cell_population.UpdateNodeLocations(forces, time_step);

	        // Check that cells locations were correctly updated
	        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
	             cell_iter != cell_population.End();
	             ++cell_iter)
	        {
	            unsigned i = cell_population.GetLocationIndexUsingCell(*cell_iter);
	            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
	            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
	        }

	        // Check that particles were correctly updated
	        // First, create a set of node indices corresponding to particles
			std::set<unsigned> node_indices;
			std::set<unsigned> particle_indices;

			for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
			{
				node_indices.insert(p_mesh->GetNode(i)->GetIndex());
			}

			std::set_difference(node_indices.begin(), node_indices.end(),
								location_indices.begin(), location_indices.end(),
								std::inserter(particle_indices, particle_indices.begin()));

			// Second, loop over all particles
			for (std::set<unsigned>::iterator it=particle_indices.begin();
				 it!=particle_indices.end();
				 it++)
			{
				TS_ASSERT_DELTA(cell_population.GetNode(*it)->rGetLocation()[0], old_posns[*it][0] +   (*it)*0.01*0.01, 1e-9);
				TS_ASSERT_DELTA(cell_population.GetNode(*it)->rGetLocation()[1], old_posns[*it][1] + 2*(*it)*0.01*0.01, 1e-9);
			}
	    }

	void TestCellPopulationWritersIn3dWithParticles() throw(Exception)
	    {
	        // Set up SimulationTime (needed if VTK is used)
	        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

	        // Resetting the Maximum cell Id to zero (to account for previous tests)
	        Cell::ResetMaxCellId();

	        // Create a simple 3D mesh with some particles
	        std::vector<Node<3>*> nodes;
	        nodes.push_back(new Node<3>(0,  true,  0.0, 0.0, 0.0));
	        nodes.push_back(new Node<3>(1,  true,  1.0, 1.0, 0.0));
	        nodes.push_back(new Node<3>(2,  true,  1.0, 0.0, 1.0));
	        nodes.push_back(new Node<3>(3,  true,  0.0, 1.0, 1.0));
	        nodes.push_back(new Node<3>(4,  false, 0.5, 0.5, 0.5));
	        nodes.push_back(new Node<3>(5,  false, -1.0, -1.0, -1.0));
	        nodes.push_back(new Node<3>(6,  false,  2.0, -1.0, -1.0));
	        nodes.push_back(new Node<3>(7,  false,  2.0,  2.0, -1.0));
	        nodes.push_back(new Node<3>(8,  false, -1.0,  2.0, -1.0));
	        nodes.push_back(new Node<3>(9,  false, -1.0, -1.0,  2.0));
	        nodes.push_back(new Node<3>(10, false,  2.0, -1.0,  2.0));
	        nodes.push_back(new Node<3>(11, false,  2.0,  2.0,  2.0));
	        nodes.push_back(new Node<3>(12, false, -1.0,  2.0,  2.0));

			// Convert this to a NodesOnlyMesh
			MAKE_PTR(NodesOnlyMesh<3>, p_mesh);
			p_mesh->ConstructNodesWithoutMesh(nodes);

	        std::vector<unsigned> location_indices;
	        for (unsigned index=0; index<5; index++)
	        {
	            location_indices.push_back(index);
	        }

	        // Set up cells
	        std::vector<CellPtr> cells;
	        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
	        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

	        cells[4]->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>()); // coverage

	        TS_ASSERT_EQUALS(cells[4]->HasCellProperty<ApoptoticCellProperty>(), true);

	        // Create cell population
	        NodeBasedCellPopulationWithParticles<3> cell_population(*p_mesh, cells, location_indices);

	        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "NodeBasedCellPopulationWithParticles-3");

	        // Test set methods
	        cell_population.SetOutputCellVolumes(true);
	        cell_population.SetOutputCellMutationStates(true);
	        cell_population.SetOutputCellProliferativeTypes(true);
	        cell_population.SetOutputCellAges(true);
	        cell_population.SetOutputCellCyclePhases(true);

	        cell_population.SetCellAncestorsToLocationIndices();
	        cell_population.SetOutputCellAncestors(true);


	        std::string output_directory = "TestCellPopulationWritersIn3dWithParticles";
	        OutputFileHandler output_file_handler(output_directory, false);

	        cell_population.CreateOutputFiles(output_directory, false);
	        cell_population.WriteResultsToFiles();
	        cell_population.CloseOutputFiles();

	        // Compare output with saved files of what they should look like
	        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

	        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes      cell_based/test/data/TestCellPopulationWritersIn3dWithParticles/results.viznodes").c_str()), 0);
	        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes  cell_based/test/data/TestCellPopulationWritersIn3dWithParticles/results.vizcelltypes").c_str()), 0);
	        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors   cell_based/test/data/TestCellPopulationWritersIn3dWithParticles/results.vizancestors").c_str()), 0);

	        // Test the GetCellMutationStateCount function: there should only be healthy cells
	        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
	        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
	        TS_ASSERT_EQUALS(cell_mutation_states[0], 5u);
	        TS_ASSERT_EQUALS(cell_mutation_states[1], 0u);
	        TS_ASSERT_EQUALS(cell_mutation_states[2], 0u);
	        TS_ASSERT_EQUALS(cell_mutation_states[3], 0u);

	        // Test the GetCellProliferativeTypeCount function - we should have 4 stem cells and 1 dead cell (for coverage)
	        std::vector<unsigned> cell_types = cell_population.rGetCellProliferativeTypeCount();
	        TS_ASSERT_EQUALS(cell_types.size(), 3u);
	        TS_ASSERT_EQUALS(cell_types[0], 5u);
	        TS_ASSERT_EQUALS(cell_types[1], 0u);
	        TS_ASSERT_EQUALS(cell_types[2], 0u);

	        // Test that the cell population parameters are output correctly
	        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

	        // Write cell population parameters to file
	        cell_population.OutputCellPopulationParameters(parameter_file);
	        parameter_file->close();

	        // Compare output with saved files of what they should look like
	        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters         cell_based/test/data/TestCellPopulationWritersIn3dWithParticles/results.parameters").c_str()), 0);

	        // Tidy up
	        CellwiseData<3>::Destroy();
	    }
};

#endif /*TESTNODEBASEDCELLPOPULATIONWITHPARTICLES_HPP_*/
