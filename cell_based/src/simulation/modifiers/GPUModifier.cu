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

#include "GPUModifier.cuh"
#include "MeshBasedCellPopulation.hpp"

FLAMEGPU_AGENT_FUNCTION(output_location, flamegpu::MessageNone, flamegpu::MessageSpatial2D) {
    FLAMEGPU->message_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x"));
    FLAMEGPU->message_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y"));
    FLAMEGPU->message_out.setVariable<float>("radius", FLAMEGPU->getVariable<float>("radius"));
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(euler_integrate, flamegpu::MessageNone, flamegpu::MessageNone) {
    float timestep = 1.0 / 120.0;
    float x = FLAMEGPU->getVariable<float>("x");
    float y = FLAMEGPU->getVariable<float>("y");
    float x_force = FLAMEGPU->getVariable<float>("x_force");
    float y_force = FLAMEGPU->getVariable<float>("y_force");
    FLAMEGPU->setVariable<float>("x", x + x_force*timestep);
    FLAMEGPU->setVariable<float>("y", y + y_force*timestep);
}

// Models repulsion force without division/apoptosis
FLAMEGPU_AGENT_FUNCTION(compute_force_meineke_spring, flamegpu::MessageSpatial2D, flamegpu::MessageNone) {
    const double x = FLAMEGPU->getVariable<float>("x");
    const double y = FLAMEGPU->getVariable<float>("y");
    float x_force = 0.0;
    float y_force = 0.0;
    float radius = FLAMEGPU->getVariable<float>("radius");

    for (const auto& message : FLAMEGPU->message_in(x, y)) {
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
                y_force += multiplication_factor * spring_stiffness * unit_y * rest_length_final* log(1.0 + overlap/rest_length_final);
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

FLAMEGPU_AGENT_FUNCTION(output_location_3D, flamegpu::MessageNone, flamegpu::MessageSpatial3D) {
    FLAMEGPU->message_out.setVariable<float>("x", FLAMEGPU->getVariable<float>("x"));
    FLAMEGPU->message_out.setVariable<float>("y", FLAMEGPU->getVariable<float>("y"));
    FLAMEGPU->message_out.setVariable<float>("z", FLAMEGPU->getVariable<float>("z"));
    FLAMEGPU->message_out.setVariable<float>("radius", FLAMEGPU->getVariable<float>("radius"));
    return flamegpu::ALIVE;
}

FLAMEGPU_AGENT_FUNCTION(euler_integrate_3D, flamegpu::MessageNone, flamegpu::MessageNone) {
    float timestep = 1.0 / 120.0;
    float x = FLAMEGPU->getVariable<float>("x");
    float y = FLAMEGPU->getVariable<float>("y");
    float z = FLAMEGPU->getVariable<float>("z");
    float x_force = FLAMEGPU->getVariable<float>("x_force");
    float y_force = FLAMEGPU->getVariable<float>("y_force");
    float z_force = FLAMEGPU->getVariable<float>("z_force");
    FLAMEGPU->setVariable<float>("x", x + x_force*timestep);
    FLAMEGPU->setVariable<float>("y", y + y_force*timestep);
    FLAMEGPU->setVariable<float>("z", z + z_force*timestep);
}

// Models repulsion force without division/apoptosis
FLAMEGPU_AGENT_FUNCTION(compute_force_meineke_spring_3D, flamegpu::MessageSpatial3D, flamegpu::MessageNone) {
    const double x = FLAMEGPU->getVariable<float>("x");
    const double y = FLAMEGPU->getVariable<float>("y");
    const double z = FLAMEGPU->getVariable<float>("z");
    float x_force = 0.0;
    float y_force = 0.0;
    float z_force = 0.0;
    float radius = FLAMEGPU->getVariable<float>("radius");

    for (const auto& message : FLAMEGPU->message_in(x, y, z)) {
        float other_x = message.getVariable<float>("x");
        float other_y = message.getVariable<float>("y");
        float other_z = message.getVariable<float>("z");
        float other_radius = message.getVariable<float>("radius");
        
        // Compute unit distance
        float x_dist = other_x - x;
        float y_dist = other_y - y;
        float z_dist = other_z - z;
        float distance_between_nodes = sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);

        float unit_x = x_dist / distance_between_nodes;
        float unit_y = y_dist / distance_between_nodes;
        float unit_z = z_dist / distance_between_nodes;
        
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
                y_force += multiplication_factor * spring_stiffness * unit_y * rest_length_final* log(1.0 + overlap/rest_length_final);
                z_force += multiplication_factor * spring_stiffness * unit_z * rest_length_final* log(1.0 + overlap/rest_length_final);
            }
            else
            {
                double alpha = 5.0;
                x_force += multiplication_factor * spring_stiffness * unit_x * overlap * exp(-alpha * overlap/rest_length_final);
                y_force += multiplication_factor * spring_stiffness * unit_y * overlap * exp(-alpha * overlap/rest_length_final);
                z_force += multiplication_factor * spring_stiffness * unit_z * overlap * exp(-alpha * overlap/rest_length_final);
            }
        }
    }

    FLAMEGPU->setVariable<float>("x_force", x_force);
    FLAMEGPU->setVariable<float>("y_force", y_force);
    FLAMEGPU->setVariable<float>("z_force", z_force);
    return flamegpu::ALIVE;
}

template<unsigned DIM>
GPUModifier<DIM>::GPUModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mpFlameGPUModel(nullptr),
    mpCellAgentDescription(nullptr),
    mpFlameGPUSimulation(nullptr)
{
}

template<unsigned DIM>
GPUModifier<DIM>::~GPUModifier()
{
}

template<unsigned DIM>
void GPUModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    auto cpu_duration = std::chrono::duration_cast<std::chrono::milliseconds>(start_time - mTimePoint);

    // Reset the simulation
    mpFlameGPUSimulation->resetStepCounter();

    // Extract cell locations from chaste
    // Get number of cells & resize agent vector to match
    unsigned int numCells = rCellPopulation.rGetMesh().GetNumNodes();
    mpCellAgentVector->resize(numCells);

    // Set the positions and clear the forces
    auto& rMesh = rCellPopulation.rGetMesh();
    auto& cellVector = *mpCellAgentVector; // Grab ref to vector for easier indexing
    unsigned int i = 0;
    for (auto iter = rMesh.GetNodeIteratorBegin(); iter != rMesh.GetNodeIteratorEnd(); ++iter) {
      auto cell = cellVector[i];
      auto& loc = iter->rGetLocation();
      cell.setVariable<float>("x", iter->rGetLocation()[0]);
      cell.setVariable<float>("y", iter->rGetLocation()[1]);
      cell.setVariable<float>("radius", 0.5f);
      cell.setVariable<float>("x_force", 0.0f);
      cell.setVariable<float>("y_force", 0.0f);

      if constexpr (DIM == 3) {
        cell.setVariable<float>("z", iter->rGetLocation()[2]);
        cell.setVariable<float>("z_force", 0.0f);
      }
      
      i++;
    }

    auto host_device_wrangle_complete_time = std::chrono::high_resolution_clock::now();
    auto host_device_wrangle_duration = std::chrono::duration_cast<std::chrono::milliseconds>(host_device_wrangle_complete_time - start_time);

    // Create cell population for FlameGPU simulation
    mpFlameGPUSimulation->setPopulationData(*mpCellAgentVector);

    auto host_device_transfer_complete_time = std::chrono::high_resolution_clock::now();
    auto host_device_duration = std::chrono::duration_cast<std::chrono::milliseconds>(host_device_transfer_complete_time - host_device_wrangle_complete_time);

    // Run the simulation
    mpFlameGPUSimulation->simulate();

    auto simulation_complete_time = std::chrono::high_resolution_clock::now();
    auto simulation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(simulation_complete_time - host_device_transfer_complete_time);

    // Extract results
    flamegpu::AgentVector out_pop(*mpCellAgentDescription);
    mpFlameGPUSimulation->getPopulationData(*mpCellAgentVector);

    auto device_host_transfer_complete_time = std::chrono::high_resolution_clock::now();
    auto device_host_duration = std::chrono::duration_cast<std::chrono::milliseconds>(device_host_transfer_complete_time - simulation_complete_time);

    // Apply results to chaste - TODO: Assumes no change in pop size. Should always be true for force resolution?
    i = 0;
    for (auto iter = rMesh.GetNodeIteratorBegin(); iter != rMesh.GetNodeIteratorEnd(); ++iter) {
        iter->rGetModifiableLocation()[0] = cellVector[i].getVariable<float>("x");
        iter->rGetModifiableLocation()[1] = cellVector[i].getVariable<float>("y");
        if constexpr (DIM == 3) {
            iter->rGetModifiableLocation()[2] = cellVector[i].getVariable<float>("z");
        }
        i++;
    }
    auto device_host_wrangle_complete_time = std::chrono::high_resolution_clock::now();
    auto device_host_wrangle_duration = std::chrono::duration_cast<std::chrono::milliseconds>(device_host_wrangle_complete_time - device_host_transfer_complete_time);

    mTimePoint = device_host_wrangle_complete_time;
    auto total_duration = host_device_wrangle_duration + host_device_duration + simulation_duration + device_host_duration + device_host_wrangle_duration;
    mTimingInfo.push_back(std::array<long, 7>{host_device_wrangle_duration.count(), host_device_duration.count(), simulation_duration.count(), device_host_duration.count(), device_host_wrangle_duration.count(), total_duration.count(), cpu_duration.count()});
    
}

template<unsigned DIM>
void GPUModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    mpFlameGPUModel = std::make_unique<flamegpu::ModelDescription>("ForceResolutionModel");
    
    // Define an agent
    mpCellAgentDescription = std::make_unique<flamegpu::AgentDescription>(mpFlameGPUModel->newAgent("cell"));
    mpCellAgentDescription->newVariable<float>("x");
    mpCellAgentDescription->newVariable<float>("y");
    mpCellAgentDescription->newVariable<float>("radius");
    mpCellAgentDescription->newVariable<float>("x_force");
    mpCellAgentDescription->newVariable<float>("y_force");
    if constexpr (DIM == 3) {
        mpCellAgentDescription->newVariable<float>("z");
        mpCellAgentDescription->newVariable<float>("z_force");
    }
    
    // Define the location message
    flamegpu::MessageSpatial2D::Description location_message = mpFlameGPUModel->newMessage<flamegpu::MessageSpatial2D>("location_message");
    flamegpu::MessageSpatial3D::Description location_message_3D = mpFlameGPUModel->newMessage<flamegpu::MessageSpatial3D>("location_message_3D");
    location_message.newVariable<float>("radius");
    location_message_3D.newVariable<float>("radius");

    location_message.setMin(-500.0, -500.0);
    location_message.setMax(500.0, 500.0);
    location_message.setRadius(1.5);

    location_message_3D.setMin(-500.0, -500.0, -500.0);
    location_message_3D.setMax(500.0, 500.0, 500.0);
    location_message_3D.setRadius(1.5);

    // Agent functions
    flamegpu::AgentFunctionDescription output_location_desc = mpCellAgentDescription->newFunction("output_location", output_location);
    flamegpu::AgentFunctionDescription compute_force_desc = mpCellAgentDescription->newFunction("compute_force_meineke_spring", compute_force_meineke_spring);
    compute_force_desc.dependsOn(output_location_desc);
    flamegpu::AgentFunctionDescription integrate_desc = mpCellAgentDescription->newFunction("euler_integrate", euler_integrate);
    integrate_desc.dependsOn(compute_force_desc);
    output_location_desc.setMessageOutput("location_message");
    compute_force_desc.setMessageInput("location_message");
    
    // Agent functions
    flamegpu::AgentFunctionDescription output_location_desc_3D = mpCellAgentDescription->newFunction("output_location_3D", output_location_3D);
    flamegpu::AgentFunctionDescription compute_force_desc_3D = mpCellAgentDescription->newFunction("compute_force_meineke_spring_3D", compute_force_meineke_spring_3D);
    compute_force_desc_3D.dependsOn(output_location_desc_3D);
    flamegpu::AgentFunctionDescription integrate_desc_3D = mpCellAgentDescription->newFunction("euler_integrate_3D", euler_integrate_3D);
    integrate_desc_3D.dependsOn(compute_force_desc_3D);
    output_location_desc_3D.setMessageOutput("location_message_3D");
    compute_force_desc_3D.setMessageInput("location_message_3D");

    // Set execution root
    if constexpr (DIM == 2) {
        mpFlameGPUModel->addExecutionRoot(output_location_desc);
    } else {
        mpFlameGPUModel->addExecutionRoot(output_location_desc_3D);
    }
    
    // Generate execution plan
    mpFlameGPUModel->generateLayers();
      
    // Construct a simulation object from the model and configure it to run for a single step
    mpFlameGPUSimulation = std::make_unique<flamegpu::CUDASimulation>(*mpFlameGPUModel);
    mpFlameGPUSimulation->SimulationConfig().steps = 1;
    
    // Allocate a vector for transferring agent data between host & device
    mpCellAgentVector = std::make_unique<flamegpu::AgentVector>(*mpCellAgentDescription);
}


template<unsigned DIM>
void GPUModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class GPUModifier<1>;
template class GPUModifier<2>;
template class GPUModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GPUModifier)

