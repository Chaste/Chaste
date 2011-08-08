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
#ifndef TESTPERIODICFORCES_HPP_
#define TESTPERIODICFORCES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "GeneralisedPeriodicLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellBasedSimulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MutableMesh.hpp"
#include "Debug.hpp"

class TestPeriodicForces : public AbstractCellBasedTestSuite
{
public:

    void TestBoundaryConditions() throw (Exception)
    {
    	// Create a simple mesh
		unsigned num_cells_depth = 2;
		unsigned num_cells_width = 2;
		HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

		// Create cell population
		MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

		// Get initial width of the cell population
		double initial_width = cell_population.GetWidth(0);

        // Create force
        GeneralisedPeriodicLinearSpringForce<2> linear_force;
        linear_force.SetInitialWidth(initial_width);

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        linear_force.AddForceContribution(node_forces, cell_population);

        TS_ASSERT_DELTA(node_forces[0][0], 7.5, 1e-3);
		TS_ASSERT_DELTA(node_forces[0][1], -2.00962, 1e-3);
		TS_ASSERT_DELTA(node_forces[1][0], -7.5, 1e-3);
		TS_ASSERT_DELTA(node_forces[1][1], -2.88444e-15, 1e-3);
		TS_ASSERT_DELTA(node_forces[2][0], 7.5, 1e-3);
		TS_ASSERT_DELTA(node_forces[2][1], 2.88444e-15, 1e-3);
		TS_ASSERT_DELTA(node_forces[3][0], -7.5, 1e-3);
		TS_ASSERT_DELTA(node_forces[3][1], 2.00962, 1e-3);

        // Set up simulation
		CellBasedSimulation<2> simulator(cell_population);
		simulator.AddForce(&linear_force);
		simulator.SetOutputDirectory("TestPeriodicBoundaryConditions");
		simulator.SetEndTime(0.01);

		simulator.Solve();

		SimulationTime::Destroy();
    }

    void TestPeriodicForceOnNonUniformMesh() throw(Exception)
	{
		// Create 2D mesh
		std::vector<Node<2> *> nodes;
		nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
		nodes.push_back(new Node<2>(1, true, 4.0, 0.0));
		nodes.push_back(new Node<2>(2, true, 1.5, 0.5));
		nodes.push_back(new Node<2>(3, true, 3.0, 0.5));
		nodes.push_back(new Node<2>(4, false, 4.5, 0.5));
		nodes.push_back(new Node<2>(5, false, 1.0, 1.0));
		nodes.push_back(new Node<2>(6, false, 3.5, 1.0));
		MutableMesh<2,2> mesh(nodes);

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

		// Create cell population
		MeshBasedCellPopulation<2> cell_population(mesh, cells);

		// Create Voronoi tessellation
		cell_population.CreateVoronoiTessellation();

		// Set up simulation
		CellBasedSimulation<2> simulator(cell_population);

		simulator.SetOutputDirectory("TestNewPeriodicForceClass");
		simulator.SetEndTime(0.01);

        // Create force
        GeneralisedPeriodicLinearSpringForce<2> linear_force;
        linear_force.SetInitialWidth(4.5);

		simulator.AddForce(&linear_force);

		// Initialise a vector of node forces
		std::vector<c_vector<double, 2> > node_forces;
		node_forces.reserve(cell_population.GetNumNodes());

		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			 node_forces.push_back(zero_vector<double>(2));
		}

		linear_force.AddForceContribution(node_forces, cell_population);

		// Close enough
		TS_ASSERT_DELTA(node_forces[0][0], 60.7698, 1e-3);
		TS_ASSERT_DELTA(node_forces[0][1], -4.74342, 1e-3);
		TS_ASSERT_DELTA(node_forces[1][0], -80.7733, 1e-3);
		TS_ASSERT_DELTA(node_forces[1][1], 3.82704, 1e-3);
		TS_ASSERT_DELTA(node_forces[2][0], 17.6281, 1e-3);
		TS_ASSERT_DELTA(node_forces[2][1], -10.4214, 1e-3);
		TS_ASSERT_DELTA(node_forces[3][0], -24.4709, 1e-3);
		TS_ASSERT_DELTA(node_forces[3][1], -0.0364322, 1e-3);
		TS_ASSERT_DELTA(node_forces[4][0], 10.6066, 1e-3);
		TS_ASSERT_DELTA(node_forces[4][1], 12.1902, 1e-3);
		TS_ASSERT_DELTA(node_forces[5][0], 18.2577, 1e-3);
		TS_ASSERT_DELTA(node_forces[5][1], -1.54716, 1e-3);
		TS_ASSERT_DELTA(node_forces[6][0], -2.01801, 1e-3);
		TS_ASSERT_DELTA(node_forces[6][1], 0.731214, 1e-3);

		// Run simulation
		simulator.Solve();

		SimulationTime::Destroy();

	}

};

#endif /* TESTPERIODICFORCES_HPP_ */
