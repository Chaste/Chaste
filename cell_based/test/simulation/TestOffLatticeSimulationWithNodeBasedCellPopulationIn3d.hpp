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

#ifndef TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATIONIN3D_HPP_
#define TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATIONIN3D_HPP_

#include <cxxtest/TestSuite.h>


// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"

class TestOffLatticeSimulationWithNodeBasedCellPopulationIn3d : public AbstractCellBasedTestSuite
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

    void Test3dNodeBasedRestrictedToSphere() throw (Exception)
    {
        // Create mesh
        std::vector<Node<3>*> nodes;

        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), TRANSIT);

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);
        node_based_cell_population.SetOutputCellProliferativeTypes(true);
        node_based_cell_population.SetOutputCellAncestors(true);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("NodeBased3dOnSphere");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,3> centre = zero_vector<double>(3);
        centre(2) = 1.0;
        double radius = 1.0;
        MAKE_PTR_ARGS(SphereGeometryBoundaryCondition<3>, p_boundary_condition, (&node_based_cell_population, centre, radius)); // Circle radius 1 centre (0,0,1)
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        // Run simulation
        simulator.Solve();

        // Check some results
        for (AbstractCellPopulation<3>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
                      cell_iter != simulator.rGetCellPopulation().End();
                      ++cell_iter)
        {
            c_vector<double,3> node_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);

            TS_ASSERT_DELTA(norm_2(node_location-centre), radius, 1e-3);
        }

        // clean up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATIONIN3D_HPP_*/
