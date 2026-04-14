/*

Copyright (c) 2005-2024, University of Oxford.
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

#ifndef TESTOFFLATTICESIMULATION_HPP_
#define TESTOFFLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <cmath>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ChemotacticForce.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulationWithMyStoppingEvent.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellIdWriter.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "GPUModifier.cuh"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "flamegpu/flamegpu.h"

#include "PetscSetupAndFinalize.hpp"

FLAMEGPU_AGENT_FUNCTION(test_do_nothing, flamegpu::MessageNone, flamegpu::MessageNone) {
    return flamegpu::ALIVE;
}

FLAMEGPU_INIT_FUNCTION(test_simple_force_create_agents) {
  // Retrieve the host agent tools for agent sheep in the default state
  flamegpu::HostAgentAPI cell = FLAMEGPU->agent("cell");

  // Create 10 new cell agents
  for (int i = 0; i < 3; ++i) {
      flamegpu::HostNewAgentAPI new_cell = cell.newAgent();
      new_cell.setVariable<float>("x", i * 0.5f);
      new_cell.setVariable<float>("y", i * 0.5f);
      new_cell.setVariable<float>("x_force", 0.0f);
      new_cell.setVariable<float>("y_force", 0.0f);
      new_cell.setVariable<float>("radius", 0.5f);
  }
}

FLAMEGPU_AGENT_FUNCTION(test_output_location, flamegpu::MessageNone, flamegpu::MessageBruteForce) {
    FLAMEGPU->message_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x"));
    FLAMEGPU->message_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y"));
    FLAMEGPU->message_out.setVariable<float>("radius", FLAMEGPU->getVariable<float>("radius"));
    return flamegpu::ALIVE;
}

// Models repulsion force without division/apoptosis
FLAMEGPU_AGENT_FUNCTION(test_compute_force_meineke_spring, flamegpu::MessageBruteForce, flamegpu::MessageNone) {
    const float x = FLAMEGPU->getVariable<float>("x");
    const float y = FLAMEGPU->getVariable<float>("y");
    float x_force = 0.0;
    float y_force = 0.0;
    float radius = FLAMEGPU->getVariable<float>("radius");

    for (const auto& message : FLAMEGPU->message_in) {
        float other_x = message.getVariable<float>("x");
        float other_y = message.getVariable<float>("y");
        float other_radius = message.getVariable<float>("radius");
        
        // Compute unit distance
        float x_dist = other_x - x;
        float y_dist = other_y - y;
        float distance_between_nodes = sqrt(x_dist * x_dist + y_dist * y_dist);

        float unit_x = x_dist / distance_between_nodes;
        float unit_y = y_dist / distance_between_nodes;
        
        // Only compute force if within cutoff distance and for positive distance
        const float cutoff_length = 1.5f;
        if (distance_between_nodes < cutoff_length && distance_between_nodes > 0.0f) {

            // Compute rest length
            const float rest_length = radius + other_radius; 
            const float rest_length_final = rest_length;
            
            // TODO: Should check here if newly divided or apoptosis happening


            // Compute the force
            float overlap = distance_between_nodes - rest_length;
            bool is_closer_than_rest_length = (overlap <= 0);
            const float spring_stiffness = 15.0f;
            const float multiplication_factor = 1.0f;

            
            // A reasonably stable simple force law
            if (is_closer_than_rest_length) //overlap is negative
            {
                //assert(overlap > -rest_length_final);
                x_force += multiplication_factor * spring_stiffness * unit_x * rest_length_final* log(1.0 + overlap/rest_length_final);
                y_force  = multiplication_factor * spring_stiffness * unit_y * rest_length_final* log(1.0 + overlap/rest_length_final);
            }
            else
            {
                double alpha = 5.0;
                x_force += multiplication_factor * spring_stiffness * unit_x * overlap * exp(-alpha * overlap/rest_length_final);
                y_force += multiplication_factor * spring_stiffness * unit_y * overlap * exp(-alpha * overlap/rest_length_final);
            }
        }

        
    }

    FLAMEGPU->setVariable<float>("x_force", x_force);        
    FLAMEGPU->setVariable<float>("y_force", y_force);        

    return flamegpu::ALIVE;
}

class TestGPUOffLatticeSimulation : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestFlameGPU2Simulation()
    {

        //flamegpu::ModelDescription model("Chaste Test");
        //
        //// Define an agent
        //flamegpu::AgentDescription agent = model.newAgent("cell"); 
        //agent.newVariable<float>("x");

        //// Agent functions
        //flamegpu::AgentFunctionDescription func = agent.newFunction("do_nothing", test_do_nothing);
        //
        //// Set execution root
        //model.addExecutionRoot(func);
        //
        ////model.addInitFunction(create_agents);
        //
        //model.generateLayers();

        //flamegpu::CUDASimulation cuda_model(model);
        //cuda_model.simulate();
        //cuda_model.reset(true);
    }
    

    // This test aims to perform a simple force calculation using both chaste and fgpu2
    // and show that the results are comparable
    void TestSimpleForceCalculation()
    {
        /* 
         * Chaste computation 
         */
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        // Create a NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.5, 0.5));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 100.0);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update(); //Needs to be called separately as not in a simulation

        GeneralisedLinearSpringForce<2> generalsiedLinearSpringForce;

        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }
        generalsiedLinearSpringForce.AddForceContribution(cell_population);

        /* 
         * Flame computation 
         */
        flamegpu::ModelDescription model("TestSimpleForceCalculation");
        
        // Define an agent
        flamegpu::AgentDescription cell_agent = model.newAgent("cell"); 
        cell_agent.newVariable<float>("x");
        cell_agent.newVariable<float>("y");
        cell_agent.newVariable<float>("radius");
        cell_agent.newVariable<float>("x_force");
        cell_agent.newVariable<float>("y_force");
        
        // Define the location message
        flamegpu::MessageBruteForce::Description location_message = model.newMessage<flamegpu::MessageBruteForce>("location_message");
        location_message.newVariable<float>("x");
        location_message.newVariable<float>("y");
        location_message.newVariable<float>("radius");

        // Agent functions
        flamegpu::AgentFunctionDescription output_location_desc = cell_agent.newFunction("test_output_location", test_output_location);
        output_location_desc.setMessageOutput("location_message");
        
        flamegpu::AgentFunctionDescription compute_force_desc = cell_agent.newFunction("test_compute_force_meineke_spring", test_compute_force_meineke_spring);
        compute_force_desc.setMessageInput("location_message");

        compute_force_desc.dependsOn(output_location_desc);
        
        // Set execution root
        model.addExecutionRoot(output_location_desc);
        
        model.addInitFunction(test_simple_force_create_agents);
        
        model.generateLayers();

        flamegpu::CUDASimulation cuda_model(model);
        cuda_model.SimulationConfig().steps = 1;
        cuda_model.simulate();
        
        // Get results
        flamegpu::AgentVector out_pop(cell_agent);
        cuda_model.getPopulationData(out_pop);
        
        /*
         * Compare forces
         */
        
        for (int i = 0; i < 3; i++) {
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], out_pop[i].getVariable<float>("x_force"), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], out_pop[i].getVariable<float>("y_force"), 1e-4);
        }
    }
    
    void TestGPUModifier() {
        
        double size_of_box = 250.0;
        unsigned cells_across = 380;
        double scaling = size_of_box/(double(cells_across-1));

        // Create a simple 3D NodeBasedCellPopulation consisting of cells evenly spaced in a regular grid
        std::vector<Node<2>*> nodes;
        unsigned index = 0;
        for (unsigned i=0; i<cells_across; i++)
        {
            for (unsigned j=0; j<cells_across; j++)
            {
                nodes.push_back(new Node<2>(index, false,  (double) i * scaling , (double) j * scaling));
                index++;
            }
        }

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        //node_based_cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("GPUNodeBased");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(1.0);

        MAKE_PTR(GPUModifier<2>, gpuModifier);
        simulator.AddSimulationModifier(gpuModifier);

        // Run simulation
        simulator.Solve();

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCPUPathway() {
        
        double size_of_box = 250.0;
        unsigned cells_across = 380;
        double scaling = size_of_box/(double(cells_across-1));

        // Create a simple 3D NodeBasedCellPopulation consisting of cells evenly spaced in a regular grid
        std::vector<Node<2>*> nodes;
        unsigned index = 0;
        for (unsigned i=0; i<cells_across; i++)
        {
            for (unsigned j=0; j<cells_across; j++)
            {
                nodes.push_back(new Node<2>(index, false,  (double) i * scaling , (double) j * scaling));
                index++;
            }
        }

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        //node_based_cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("GPUNodeBased");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, springForce);
        simulator.AddForce(springForce);

        // Run simulation
        simulator.Solve();

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
  
};

#endif /*TESTGPUOFFLATTICESIMULATION_HPP_*/
