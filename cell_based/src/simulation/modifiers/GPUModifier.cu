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
      cellVector[i].setVariable<float>("x", iter->rGetLocation()[0]);
      cellVector[i].setVariable<float>("y", iter->rGetLocation()[1]);
      cellVector[i].setVariable<float>("radius", 1.5f);
      cellVector[i].setVariable<float>("x_force", 0.0f);
      cellVector[i].setVariable<float>("y_force", 0.0f);
      i++;
    }

    // Create cell population for FlameGPU simulation
    mpFlameGPUSimulation->setPopulationData(*mpCellAgentVector);

    // Run the simulation
    mpFlameGPUSimulation->simulate();

    // Extract results
    flamegpu::AgentVector out_pop(*mpCellAgentDescription);
    mpFlameGPUSimulation->getPopulationData(*mpCellAgentVector);

    // Apply results to chaste - TODO: Assumes no change in pop size. Should always be true for force resolution?
    i = 0;
    for (auto iter = rMesh.GetNodeIteratorBegin(); iter != rMesh.GetNodeIteratorEnd(); ++iter) {
        iter->rGetModifiableLocation()[0] = cellVector[i].getVariable<float>("x");
        iter->rGetModifiableLocation()[1] = cellVector[i].getVariable<float>("y");
        i++;
    }
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
    
    // Define the location message
    flamegpu::MessageSpatial2D::Description location_message = mpFlameGPUModel->newMessage<flamegpu::MessageSpatial2D>("location_message");
    //location_message.newVariable<float>("x"); // Implicit for spatial message
    //location_message.newVariable<float>("y"); // Implicit for spatial message
    location_message.newVariable<float>("radius");
    location_message.setMin(-500.0, -500.0);
    location_message.setMax(500.0, 500.0);
    location_message.setRadius(1.5);

    // Agent functions
    flamegpu::AgentFunctionDescription output_location_desc = mpCellAgentDescription->newFunction("output_location", output_location);
    output_location_desc.setMessageOutput("location_message");
    
    flamegpu::AgentFunctionDescription compute_force_desc = mpCellAgentDescription->newFunction("csfompute_force_meineke_spring", compute_force_meineke_spring);
    compute_force_desc.setMessageInput("location_message");

    compute_force_desc.dependsOn(output_location_desc);
    
    // Set execution root
    mpFlameGPUModel->addExecutionRoot(output_location_desc);
    
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

