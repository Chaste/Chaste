/*

Copyright (c) 2005-2023, University of Oxford.
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
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "SphereGeometryBoundaryCondition.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "PeriodicNodesOnlyMesh.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

// Cell population writers
#include "NodeVelocityWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestOffLatticeSimulationWithNodeBasedCellPopulationIn3d : public AbstractCellBasedWithTimingsTestSuite
{

private:

    /**
     * Helper method to generate a simple mesh that works in parallel
     *
     * @param nx number of nodes in the x direction
     * @param ny number of nodes in the y direction
     * @param nz number of nodes in the z direction
     * @return the nodes for the nodes only mesh
     */
    std::vector<Node<3>*> GenerateMesh( unsigned nx, unsigned ny, unsigned nz )
    {
        std::vector<Node<3>*> nodes(nx*ny*nz);
        for ( unsigned k = 0; k < nz; k++ )
        {
            for ( unsigned j = 0; j < ny; j++ )
            {
                for ( unsigned i = 0; i < nx; i++ )
                {
                    double x = (double)i + 0.5*(double)(j%2) + 0.5*(double)(k%2) - (double)( ((j+2)*k)%2);
                    double y = sqrt(3.0)/2.0 * (double)j;
                    double z = sqrt(3.0)/2.0 * (double)k;
                    nodes[ k*nx*ny + j*nx + i ] = new Node<3>(i+j*nx+k*nx*ny, false, x, y, z );
                }
            }
        }
        return nodes;
    }

public:

    void Test3dNodeBasedRestrictedToSphere()
    {
        // Create a simple 3D NodeBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);

        // Set output options
        node_based_cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        node_based_cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("NodeBased3dOnSphere");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(10.0); // 50.0

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,3> centre = zero_vector<double>(3);
        centre(2) = 1.0;
        double radius = 2.0;
        MAKE_PTR_ARGS(SphereGeometryBoundaryCondition<3>, p_boundary_condition, (&node_based_cell_population, centre, radius)); // Circle radius 1 centre (0,0,1)
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        // Kill all cells moving past z=1;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer,(&node_based_cell_population, unit_vector<double>(3,2), unit_vector<double>(3,2)));
        simulator.AddCellKiller(p_cell_killer);

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

        std::vector<unsigned> mutation_state_count_after_solve = simulator.rGetCellPopulation().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(mutation_state_count_after_solve.size(), 4u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[0], 4u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[1], 0u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[2], 0u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[3], 0u);

        std::vector<unsigned> prolif_type_count_after_solve = simulator.rGetCellPopulation().GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(prolif_type_count_after_solve.size(), 4u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[0], 0u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[1], 2u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[2], 0u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[3], 0u);

        std::vector<unsigned> cell_cycle_phase_count = simulator.rGetCellPopulation().GetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phase_count.size(), 5u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count[0], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count[1], 1u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count[2], 3u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count[3], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phase_count[4], 0u);

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void Test3dNodeBasedPlaneBoundary()
    {
        EXIT_IF_PARALLEL;    // Output doesn't work in parallel so we cannot solve a simulation #2365

        // Create a simple 3D NodeBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, -1.0, 0.0, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);
        node_based_cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

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

            TS_ASSERT_LESS_THAN_EQUALS(-0.5, node_location[0]);
            TS_ASSERT_LESS_THAN_EQUALS(node_location[0], 0.5);
        }

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    // Now Test Periodicity
    void TestSimple3dTissueXPeriodic()
    {
        // Set up the node positions
        std::vector<Node<3>*> nodes = GenerateMesh(3,3,3);

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,3> periodic_width = zero_vector<double>(3);
        periodic_width[0] = 6.0; //periodic in x with width 6.0
        PeriodicNodesOnlyMesh<3> mesh(periodic_width);
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("Test3dOffLatticeSimulationWithXPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic bounday influencing the cells
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't outside the domain
        // Use node iterators to make sure we only get the nodes from this process for parallel implementation
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_begin = mesh.GetNodeIteratorBegin();
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_end = mesh.GetNodeIteratorEnd();
        for ( PeriodicNodesOnlyMesh<3>::NodeIterator node_iter = node_iter_begin; node_iter != node_iter_end; ++node_iter )
        {
            double x_node = (*node_iter).rGetLocation()[0];
            TS_ASSERT_LESS_THAN_EQUALS(0,x_node);
            TS_ASSERT_LESS_THAN_EQUALS(x_node, periodic_width[0]);
        }

        // Now run the simulation again with the periodic boundary in a different place and check its the same

        // First reset the singletons
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);

        // Convert this to a PeriodicNodesOnlyMesh
        PeriodicNodesOnlyMesh<3> mesh_2(periodic_width);
        mesh_2.ConstructNodesWithoutMesh(nodes, 1.5);

        // Add an offset
        double x_offset = periodic_width[0]/2.0;
        mesh_2.Translate(-x_offset,0.0);

        // Create cells
        std::vector<CellPtr> cells_2;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_2;
        cells_generator_2.GenerateBasicRandom(cells_2, mesh_2.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population_2(mesh_2, cells_2);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator_2(node_based_cell_population_2);
        simulator_2.SetOutputDirectory("Test3dOffLatticeSimulationWith2ndXPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic boundary influencing the cells
        simulator_2.SetSamplingTimestepMultiple(120);
        simulator_2.SetEndTime(5.0);

        // Pass the same force law to the simulation
        simulator_2.AddForce(p_linear_force);

        simulator_2.Solve();

        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check that cells are in the same place in each simulation
            // Note the way we do this only works on 1 processor
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                double x_1 = simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0];
                double x_2 = simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[0];

                if (x_1 < x_offset)
                {
                    TS_ASSERT_DELTA(x_1+x_offset, x_2, 1e-6)
                }
                else
                {
                    TS_ASSERT_DELTA(x_1-x_offset, x_2, 1e-6)
                }

                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[2],1e-6);
            }
        }


        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check with a different interaction distance
            // This can only be done on a single processor as the random number generator
            // is set for each process and gives different results with different
            // numbers of nodes on a single process

            // First reset the singletons
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            RandomNumberGenerator::Instance()->Reseed(0);

            // Convert this to a Cylindrical2dNodesOnlyMesh
            PeriodicNodesOnlyMesh<3> mesh_3(periodic_width);
            mesh_3.ConstructNodesWithoutMesh(nodes, 2.0);

            // Create cells
            std::vector<CellPtr> cells_3;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator_3;
            cells_generator_3.GenerateBasicRandom(cells_3, mesh_3.GetNumNodes());

            // Create a node-based cell population
            NodeBasedCellPopulation<3> node_based_cell_population_3(mesh_3, cells_3);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator_3(node_based_cell_population_3);
            simulator_3.SetOutputDirectory("Test3dOffLatticeSimulationWith3rdXPeriodicNodeBasedCellPopulation");

            // Run for long enough to see the periodic boundary influencing the cells
            simulator_3.SetSamplingTimestepMultiple(120);
            simulator_3.SetEndTime(5.0);

            // Pass the same force law to the simulation
            simulator_3.AddForce(p_linear_force);

            simulator_3.Solve();

            // Check that cells are in the same place in each simulation
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[2],1e-6);
            }
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestSimple3dTissueYPeriodic()
    {
        // Set up the node positions
        std::vector<Node<3>*> nodes = GenerateMesh(3,3,3);

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,3> periodic_width = zero_vector<double>(3);
        periodic_width[1] = 6.0; //periodic in y with width 6.0
        PeriodicNodesOnlyMesh<3> mesh(periodic_width);
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("Test3dOffLatticeSimulationWithYPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic bounday influencing the cells
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't outside the domain
        // Use node iterators to make sure we only get the nodes from this process for parallel implementation
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_begin = mesh.GetNodeIteratorBegin();
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_end = mesh.GetNodeIteratorEnd();
        for ( PeriodicNodesOnlyMesh<3>::NodeIterator node_iter = node_iter_begin; node_iter != node_iter_end; ++node_iter )
        {
            double y_node = (*node_iter).rGetLocation()[1];
            TS_ASSERT_LESS_THAN_EQUALS(0,y_node);
            TS_ASSERT_LESS_THAN_EQUALS(y_node, periodic_width[1]);
        }

        // Now run the simulation again with the periodic boundary in a different place and check its the same

        // First reset the singletons
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);

        // Convert this to a Cylindrical2dNodesOnlyMesh
        PeriodicNodesOnlyMesh<3> mesh_2(periodic_width);
        mesh_2.ConstructNodesWithoutMesh(nodes, 1.5);

        // Add an offset
        double y_offset = periodic_width[0]/2.0;
        mesh_2.Translate(0.0,-y_offset);

        // Create cells
        std::vector<CellPtr> cells_2;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_2;
        cells_generator_2.GenerateBasicRandom(cells_2, mesh_2.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population_2(mesh_2, cells_2);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator_2(node_based_cell_population_2);
        simulator_2.SetOutputDirectory("Test3dOffLatticeSimulationWith2ndYPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic boundary influencing the cells
        simulator_2.SetSamplingTimestepMultiple(120);
        simulator_2.SetEndTime(5.0);

        // Pass the same force law to the simulation
        simulator_2.AddForce(p_linear_force);

        simulator_2.Solve();

        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check that cells are in the same place in each simulation
            // Note the way we do this only works on 1 processor
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                double y_1 = simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1];
                double y_2 = simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[1];

                if (y_1 < y_offset)
                {
                    TS_ASSERT_DELTA(y_1+y_offset, y_2, 1e-6)
                }
                else
                {
                    TS_ASSERT_DELTA(y_1-y_offset, y_2, 1e-6)
                }

                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[2],1e-6);
            }
        }

        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check with a different interaction distance
            // This can only be done on a single processor as the random number generator
            // is set for each process and gives different results with different
            // numbers of nodes on a single process

            // First reset the singletons
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            RandomNumberGenerator::Instance()->Reseed(0);

            // Convert this to a Cylindrical2dNodesOnlyMesh
            PeriodicNodesOnlyMesh<3> mesh_3(periodic_width);
            mesh_3.ConstructNodesWithoutMesh(nodes, 2.0);

            // Create cells
            std::vector<CellPtr> cells_3;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator_3;
            cells_generator_3.GenerateBasicRandom(cells_3, mesh_3.GetNumNodes());

            // Create a node-based cell population
            NodeBasedCellPopulation<3> node_based_cell_population_3(mesh_3, cells_3);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator_3(node_based_cell_population_3);
            simulator_3.SetOutputDirectory("Test3dOffLatticeSimulationWith3rdPeriodicNodeBasedCellPopulation");

            // Run for long enough to see the periodic boundary influencing the cells
            simulator_3.SetSamplingTimestepMultiple(120);
            simulator_3.SetEndTime(5.0);

            // Pass the same force law to the simulation
            simulator_3.AddForce(p_linear_force);

            simulator_3.Solve();

            // Check that cells are in the same place in each simulation
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[2],1e-6);
            }
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestSimple3dTissueZPeriodic()
    {
        // Set up the node positions
        std::vector<Node<3>*> nodes = GenerateMesh(3,3,3);

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,3> periodic_width = zero_vector<double>(3);
        periodic_width[2] = 6.0; //periodic in y with width 6.0
        PeriodicNodesOnlyMesh<3> mesh(periodic_width);
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("Test3dOffLatticeSimulationWithZPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic bounday influencing the cells
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

         simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't outside the domain
        // Use node iterators to make sure we only get the nodes from this process for parallel implementation
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_begin = mesh.GetNodeIteratorBegin();
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_end = mesh.GetNodeIteratorEnd();
        for ( PeriodicNodesOnlyMesh<3>::NodeIterator node_iter = node_iter_begin; node_iter != node_iter_end; ++node_iter )
        {
            double z_node = (*node_iter).rGetLocation()[2];
            TS_ASSERT_LESS_THAN_EQUALS(0,z_node);
            TS_ASSERT_LESS_THAN_EQUALS(z_node, periodic_width[2]);
        }

        // Now run the simulation again with the periodic boundary in a different place and check its the same

        // First reset the singletons
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);

        // Convert this to a Cylindrical2dNodesOnlyMesh
        PeriodicNodesOnlyMesh<3> mesh_2(periodic_width);
        mesh_2.ConstructNodesWithoutMesh(nodes, 1.5);

        // Add an offset
        double z_offset = periodic_width[0]/2.0;
        mesh_2.Translate(0.0,0.0,-z_offset);

        // Create cells
        std::vector<CellPtr> cells_2;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_2;
        cells_generator_2.GenerateBasicRandom(cells_2, mesh_2.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population_2(mesh_2, cells_2);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator_2(node_based_cell_population_2);
        simulator_2.SetOutputDirectory("Test3dOffLatticeSimulationWith2ndZPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic boundary influencing the cells
        simulator_2.SetSamplingTimestepMultiple(120);
        simulator_2.SetEndTime(5.0);

        // Pass the same force law to the simulation
        simulator_2.AddForce(p_linear_force);

        simulator_2.Solve();

        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check that cells are in the same place in each simulation
            // Note the way we do this only works on 1 processor
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                double z_1 = simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2];
                double z_2 = simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[2];

                if (z_1 < z_offset)
                {
                    TS_ASSERT_DELTA(z_1+z_offset, z_2, 1e-6)
                }
                else
                {
                    TS_ASSERT_DELTA(z_1-z_offset, z_2, 1e-6)
                }

                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
            }
        }

        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check with a different interaction distance
            // This can only be done on a single processor as the random number generator
            // is set for each process and gives different results with different
            // numbers of nodes on a single process

            // First reset the singletons
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            RandomNumberGenerator::Instance()->Reseed(0);

            // Convert this to a Cylindrical2dNodesOnlyMesh
            PeriodicNodesOnlyMesh<3> mesh_3(periodic_width);
            mesh_3.ConstructNodesWithoutMesh(nodes, 2.0);

            // Create cells
            std::vector<CellPtr> cells_3;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator_3;
            cells_generator_3.GenerateBasicRandom(cells_3, mesh_3.GetNumNodes());

            // Create a node-based cell population
            NodeBasedCellPopulation<3> node_based_cell_population_3(mesh_3, cells_3);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator_3(node_based_cell_population_3);
            simulator_3.SetOutputDirectory("Test3dOffLatticeSimulationWith3rdZPeriodicNodeBasedCellPopulation");

            // Run for long enough to see the periodic boundary influencing the cells
            simulator_3.SetSamplingTimestepMultiple(120);
            simulator_3.SetEndTime(5.0);

            // Pass the same force law to the simulation
            simulator_3.AddForce(p_linear_force);

            simulator_3.Solve();

            // Check that cells are in the same place in each simulation
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[2],1e-6);
            }
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

    }

    void TestSimple3dTissueXZPeriodic()
    {
        // Set up the node positions
        std::vector<Node<3>*> nodes = GenerateMesh(3,3,3);

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,3> periodic_width = zero_vector<double>(3);
        periodic_width[0] = 6.0;//periodic in x with width 6.0
        periodic_width[2] = 6.0;//periodic in z with width 6.0
        PeriodicNodesOnlyMesh<3> mesh(periodic_width);
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("Test3dOffLatticeSimulationWithXZPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic bounday influencing the cells
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't outside the domain
        // Use node iterators to make sure we only get the nodes from this process for parallel implementation
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_begin = mesh.GetNodeIteratorBegin();
        PeriodicNodesOnlyMesh<3>::NodeIterator node_iter_end = mesh.GetNodeIteratorEnd();
        for ( PeriodicNodesOnlyMesh<3>::NodeIterator node_iter = node_iter_begin; node_iter != node_iter_end; ++node_iter )
        {
            const c_vector<double,3>& node_location = (*node_iter).rGetLocation();
            TS_ASSERT_LESS_THAN_EQUALS(0,node_location[0]);
            TS_ASSERT_LESS_THAN_EQUALS(0,node_location[2]);
            TS_ASSERT_LESS_THAN_EQUALS(node_location[0], periodic_width[0]);
            TS_ASSERT_LESS_THAN_EQUALS(node_location[2], periodic_width[2]);
        }

        // Now run the simulation again with the periodic boundary in a different place and check its the same

        // First reset the singletons
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);

        // Convert this to a Cylindrical2dNodesOnlyMesh
        PeriodicNodesOnlyMesh<3> mesh_2(periodic_width);
        mesh_2.ConstructNodesWithoutMesh(nodes, 1.5);

        // Add an offset
        double offset = periodic_width[0]/2.0;
        mesh_2.Translate(-offset,0.0,-offset);

        // Create cells
        std::vector<CellPtr> cells_2;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_2;
        cells_generator_2.GenerateBasicRandom(cells_2, mesh_2.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population_2(mesh_2, cells_2);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator_2(node_based_cell_population_2);
        simulator_2.SetOutputDirectory("Test3dOffLatticeSimulationWith2ndXZPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic boundary influencing the cells
        simulator_2.SetSamplingTimestepMultiple(120);
        simulator_2.SetEndTime(5.0);

        // Pass the same force law to the simulation
        simulator_2.AddForce(p_linear_force);

        simulator_2.Solve();

        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check that cells are in the same place in each simulation
            // Note the way we do this only works on 1 processor
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                double x_1 = simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0];
                double x_2 = simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[0];
                double z_1 = simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2];
                double z_2 = simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[2];

                if (x_1 < offset)
                {
                    TS_ASSERT_DELTA(x_1+offset, x_2, 1e-6)
                }
                else
                {
                    TS_ASSERT_DELTA(x_1-offset, x_2, 1e-6)
                }

                if (z_1 < offset)
                {
                    TS_ASSERT_DELTA(z_1+offset, z_2, 1e-6)
                }
                else
                {
                    TS_ASSERT_DELTA(z_1-offset, z_2, 1e-6)
                }

                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
            }
        }

        if ( PetscTools::GetNumProcs() == 1 )
        {
            // Check with a different interaction distance
            // This can only be done on a single processor as the random number generator
            // is set for each process and gives different results with different
            // numbers of nodes on a single process

            // First reset the singletons
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            RandomNumberGenerator::Instance()->Reseed(0);

            // Convert this to a Cylindrical2dNodesOnlyMesh
            PeriodicNodesOnlyMesh<3> mesh_3(periodic_width);
            mesh_3.ConstructNodesWithoutMesh(nodes, 2.0);

            // Create cells
            std::vector<CellPtr> cells_3;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator_3;
            cells_generator_3.GenerateBasicRandom(cells_3, mesh_3.GetNumNodes());

            // Create a node-based cell population
            NodeBasedCellPopulation<3> node_based_cell_population_3(mesh_3, cells_3);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator_3(node_based_cell_population_3);
            simulator_3.SetOutputDirectory("Test3dOffLatticeSimulationWith3rdXZPeriodicNodeBasedCellPopulation");

            // Run for long enough to see the periodic boundary influencing the cells
            simulator_3.SetSamplingTimestepMultiple(120);
            simulator_3.SetEndTime(5.0);

            // Pass the same force law to the simulation
            simulator_3.AddForce(p_linear_force);

            simulator_3.Solve();

            // Check that cells are in the same place in each simulation
            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
                TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[2],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[2],1e-6);
            }

        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

};

#endif /*TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATIONIN3D_HPP_*/
