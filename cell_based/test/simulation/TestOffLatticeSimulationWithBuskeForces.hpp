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

#ifndef TESTOFFLATTICESIMULATIONWITHBUSKEFORCES_HPP_
#define TESTOFFLATTICESIMULATIONWITHBUSKEFORCES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestOffLatticeSimulationWithBuskeForces : public AbstractCellBasedWithTimingsTestSuite
{
public:

    /**
     * Create a simulation of a NodeBasedCellPopulation with a BuskeInteractionForce system.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayerWithBuskeAdhesiveForce()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithBuskeAdhesiveForce");
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(BuskeAdhesiveForce<2>, p_buske_adhesive_force);
        p_buske_adhesive_force->SetAdhesionEnergyParameter(0.002);
        simulator.AddForce(p_buske_adhesive_force);

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

        TS_ASSERT_LESS_THAN_EQUALS(1e-3, min_distance_between_cells);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with a BuskeElasticForce system.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayerWithBuskeElasticForce()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithBuskeElasticForce");
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(BuskeElasticForce<2>, p_buske_elastic_force);
        simulator.AddForce(p_buske_elastic_force);

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

        TS_ASSERT_LESS_THAN_EQUALS(1e-3, min_distance_between_cells);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with a BuskeCompressionForce system.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayerWithBuskeCompressionForce()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithBuskeCompressionForce");
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(BuskeCompressionForce<2>, buske_compression_force);
        buske_compression_force->SetCompressionEnergyParameter(0.01);
        simulator.AddForce(buske_compression_force);

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
    void TestAllBuskeForces()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestAllBuskeForces");
        simulator.SetEndTime(5.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(BuskeCompressionForce<2>, p_buske_compression_force);
        MAKE_PTR(BuskeElasticForce<2>, p_buske_elastic_force);
        MAKE_PTR(BuskeAdhesiveForce<2>, p_buske_adhesive_force);
        p_buske_compression_force->SetCompressionEnergyParameter(0.01);
        simulator.AddForce(p_buske_compression_force);
        simulator.AddForce(p_buske_elastic_force);
        simulator.AddForce(p_buske_adhesive_force);

        simulator.Solve();
    }

    /**
     * Test that two nodes relax to the equilibrium distance.
     */
    void TestBuskeRelaxationForces()
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh with two nodes
        HoneycombMeshGenerator generator(2, 1, 0, 0.9);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_buske;
        cells_generator_buske.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population_buske(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population_buske);
        simulator.SetOutputDirectory("TestBuskeRelaxation");
        simulator.SetEndTime(1.0);

        // Create all three Buske force laws and pass to the simulation
        MAKE_PTR(BuskeCompressionForce<2>, p_buske_compression_force);
        MAKE_PTR(BuskeElasticForce<2>, p_buske_elastic_force);
        MAKE_PTR(BuskeAdhesiveForce<2>, p_buske_adhesive_force);
        simulator.AddForce(p_buske_compression_force);
        simulator.AddForce(p_buske_elastic_force);
        simulator.AddForce(p_buske_adhesive_force);

        // Solve
        simulator.Solve();

        // The nodes should be about 0.85 apart as this is the minimum of the sum of the energies.
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()[0], 0.0155,  1e-4);
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(1)->rGetLocation()[0], 0.8844,  1e-4);
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHBUSKEFORCES_HPP_*/
