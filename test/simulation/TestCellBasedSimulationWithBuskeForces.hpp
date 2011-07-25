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
#ifndef TESTCELLBASEDSIMULATIONWITHBUSKEFORCES_HPP_
#define TESTCELLBASEDSIMULATIONWITHBUSKEFORCES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "CellBasedSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"

class TestCellBasedSimulationWithBuskeForces : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    /**
     * Create a simulation of a NodeBasedCellPopulation with a BuskeInteractionForce system.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayerWithBuskeAdhesiveForce() throw (Exception)
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithBuskeAdhesiveForce");
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        BuskeAdhesiveForce<2> buske_adhesive_force;
        buske_adhesive_force.SetAdhesionEnergyParameter(0.002);
        simulator.AddForce(&buske_adhesive_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with a BuskeElasticForce system.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayerWithBuskeElasticForce() throw (Exception)
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithBuskeElasticForce");
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        BuskeElasticForce<2> buske_elastic_force;
        simulator.AddForce(&buske_elastic_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with a BuskeCompressionForce system.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayerWithBuskeCompressionForce() throw (Exception)
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), TRANSIT);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithBuskeCompressionForce");
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        BuskeCompressionForce<2> buske_compression_force;
        simulator.AddForce(&buske_compression_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with all Buske forces.
     * Test that no exceptions are thrown.
     */
    void TestAllBuskeForces() throw (Exception)
    {
		// Create a simple mesh
		unsigned num_cells_depth = 5;
		unsigned num_cells_width = 5;
		HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
		TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

		// Convert this to a NodesOnlyMesh
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), TRANSIT);

		// Create a node-based cell population
		NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
		node_based_cell_population.SetMechanicsCutOffLength(1.5);

		// Set up cell-based simulation
		CellBasedSimulation<2> simulator(node_based_cell_population);
		simulator.SetOutputDirectory("TestAllBuskeForces");
		simulator.SetEndTime(5.0);

		// Create a force law and pass it to the simulation
		BuskeCompressionForce<2> buske_compression_force;
		BuskeElasticForce<2> buske_elastic_force;
		BuskeAdhesiveForce<2> buske_adhesive_force;
		simulator.AddForce(&buske_compression_force);
		simulator.AddForce(&buske_elastic_force);
		simulator.AddForce(&buske_adhesive_force);

		simulator.Solve();
    }

    /**
	 * Test that 2 nodes relax to the equilibrium distance.
     */

	void TestBuskeRelaxationForces() throw (Exception)
	{
		// Create a simple mesh with two nodes
		unsigned num_cells_depth = 1;
		unsigned num_cells_width = 2;
		HoneycombMeshGenerator generator_buske(num_cells_width, num_cells_depth, 0, 1.0);
		TetrahedralMesh<2,2>* p_generating_mesh_buske = generator_buske.GetMesh();

		// Convert this to a NodesOnlyMesh
		NodesOnlyMesh<2> mesh_buske;
		mesh_buske.ConstructNodesWithoutMesh(*p_generating_mesh_buske);

		// Create cells
		std::vector<CellPtr> cells_buske;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator_buske;
		cells_generator_buske.GenerateBasicRandom(cells_buske, mesh_buske.GetNumNodes(), DIFFERENTIATED);

		// Create a node-based cell population
		NodeBasedCellPopulation<2> node_based_cell_population_buske(mesh_buske, cells_buske);
		node_based_cell_population_buske.SetMechanicsCutOffLength(1.5);

		// Set up cell-based simulation
		CellBasedSimulation<2> simulator(node_based_cell_population_buske);
		simulator.SetOutputDirectory("TestBuskeRelaxation");
		simulator.SetEndTime(1.0);

		// Create all three Buske force laws and pass to the simulation
		BuskeCompressionForce<2> buske_compression_force;
		BuskeElasticForce<2> buske_elastic_force;
		BuskeAdhesiveForce<2> buske_adhesive_force;
		simulator.AddForce(&buske_compression_force);
		simulator.AddForce(&buske_elastic_force);
		simulator.AddForce(&buske_adhesive_force);

		// Solve
		simulator.Solve();

		// The nodes should be about 1.7 apart as this is the minimum of the sum of the energies.
		TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()[0], -0.3596,  1e-4);
		TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(1)->rGetLocation()[0], 1.3596,  1e-4);
	}
};

#endif /*TESTCELLBASEDSIMULATIONWITHBUSKEFORCES_HPP_*/
