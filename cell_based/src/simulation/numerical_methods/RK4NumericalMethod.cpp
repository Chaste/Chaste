/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "RK4NumericalMethod.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RK4NumericalMethod<ELEMENT_DIM,SPACE_DIM>::RK4NumericalMethod()
    : AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RK4NumericalMethod<ELEMENT_DIM,SPACE_DIM>::~RK4NumericalMethod()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RK4NumericalMethod<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt)
{
    if (!this->mUseUpdateNodeLocation)
    {
        // Store initial node locations
        auto old_node_locations = this->SaveCurrentNodeLocations();

        // Apply boundary conditions to the initial positions, then resave
        this->ImposeBoundaryConditions(old_node_locations);
        old_node_locations = this->SaveCurrentNodeLocations();

        // Compute k1 = F(r^t) / nu
        auto k1 = this->ComputeForcesIncludingDamping();

        // Update to r^t + dt*k1/2 and compute k2
        unsigned index = 0;
        for (auto node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> new_location = node_iter->rGetLocation() + (dt/2.0) * k1[index];
            this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
        }
        this->ImposeBoundaryConditions(old_node_locations);

        // Compute k2 = F(r^t + dt*k1/2) / nu
        auto k2 = this->ComputeForcesIncludingDamping();

        // Update to r^t + dt*k2/2 (from old locations) and compute k3
        index = 0;
        for (auto node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> old_location = old_node_locations.find(&(*node_iter))->second;
            c_vector<double, SPACE_DIM> new_location = old_location + (dt/2.0) * k2[index];
            this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
        }
        this->ImposeBoundaryConditions(old_node_locations);

        // Compute k3 = F(r^t + dt*k2/2) / nu
        auto k3 = this->ComputeForcesIncludingDamping();

        // Update to r^t + dt*k3 (from old locations) and compute k4
        index = 0;
        for (auto node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> old_location = old_node_locations.find(&(*node_iter))->second;
            c_vector<double, SPACE_DIM> new_location = old_location + dt * k3[index];
            this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
        }
        this->ImposeBoundaryConditions(old_node_locations);

        // Compute k4 = F(r^t + dt*k3) / nu
        auto k4 = this->ComputeForcesIncludingDamping();

        // Final RK4 position update: r^(t+1) = r^t + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        index = 0;
        for (auto node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> effective_force = (k1[index] + 2.0*k2[index] + 2.0*k3[index] + k4[index]) / 6.0;
            c_vector<double, SPACE_DIM> old_location = old_node_locations.find(&(*node_iter))->second;
            c_vector<double, SPACE_DIM> displacement = dt * effective_force;

            this->DetectStepSizeExceptions(node_iter->GetIndex(), displacement, dt);

            c_vector<double, SPACE_DIM> new_location = old_location + displacement;
            this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
        }
    }
    else
    {
        /*
         * If this type of cell population does not support the new numerical methods, delegate
         * updating node positions to the population itself.
         *
         * This only applies to NodeBasedCellPopulationWithBuskeUpdates.
         */
        NEVER_REACHED;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RK4NumericalMethod<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(rParamsFile);
}

// Explicit instantiation
template class RK4NumericalMethod<1,1>;
template class RK4NumericalMethod<1,2>;
template class RK4NumericalMethod<2,2>;
template class RK4NumericalMethod<1,3>;
template class RK4NumericalMethod<2,3>;
template class RK4NumericalMethod<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RK4NumericalMethod)
