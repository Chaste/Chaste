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
#include "PlaneBoundaryCondition.hpp"
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
        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(*p_mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);
        node_based_cell_population.SetOutputCellProliferativeTypes(true);

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

        // Avoid memory leak
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    void Test3dNodeBasedPlaneBoundary() throw (Exception)
    {
        // Create mesh
        std::vector<Node<3>*> nodes;

        nodes.push_back(new Node<3>(0u,  false,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -1.0, 0.0, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3>* p_mesh = new NodesOnlyMesh<3>;
        p_mesh->ConstructNodesWithoutMesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(*p_mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);
        node_based_cell_population.SetOutputCellProliferativeTypes(true);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("NodeBased3dPlaneBoundary");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create a plane boundary and pass them to the simulation
        c_vector<double,3> point_on_plane = zero_vector<double>(3);
        c_vector<double,3> normal_to_plane = zero_vector<double>(3);
        point_on_plane(0) = 0.5;
        normal_to_plane(0) = 1.0;

        // Restrict to x<1/2
        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition, (&node_based_cell_population, point_on_plane, normal_to_plane));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        point_on_plane(0) = -0.5;
        normal_to_plane(0) = -1.0;

        // Restrict to x>-1/2
        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition2, (&node_based_cell_population, point_on_plane, normal_to_plane));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition2);

        // Run simulation
        simulator.Solve();

        // Check some results
        for (AbstractCellPopulation<3>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
             cell_iter != simulator.rGetCellPopulation().End();
             ++cell_iter)
        {
            c_vector<double,3> node_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);

            TS_ASSERT_LESS_THAN(-0.5, node_location[0]);
            TS_ASSERT_LESS_THAN(node_location[0], 0.5);
        }

        // Avoid memory leak
        delete p_mesh;
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

};

#endif /*TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATIONIN3D_HPP_*/
