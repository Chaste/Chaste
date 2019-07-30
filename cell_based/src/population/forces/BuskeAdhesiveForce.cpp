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

#include "BuskeAdhesiveForce.hpp"

template<unsigned DIM>
BuskeAdhesiveForce<DIM>::BuskeAdhesiveForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mAdhesionEnergyParameter(0.2)        // Denoted by epsilon in Buske et al (2011) (doi:10.1371/journal.pcbi.1001045).
{
}

template<unsigned DIM>
double BuskeAdhesiveForce<DIM>::GetAdhesionEnergyParameter()
{
    return mAdhesionEnergyParameter;
}

template<unsigned DIM>
void BuskeAdhesiveForce<DIM>::SetAdhesionEnergyParameter(double adhesionEnergyParameter)
{
    mAdhesionEnergyParameter = adhesionEnergyParameter;
}

template<unsigned DIM>
c_vector<double, DIM> BuskeAdhesiveForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                          unsigned nodeBGlobalIndex,
                                                                          AbstractCellPopulation<DIM>& rCellPopulation)
{
    // This force class is defined for NodeBasedCellPopulations only
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr);

    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    const c_vector<double, DIM>& r_node_a_location = p_node_a->rGetLocation();
    const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
    c_vector<double, DIM> unit_vector = r_node_b_location - r_node_a_location;

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_vector);

    // Account for any cutoff in the force law
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(DIM);
        }
    }

    // Assert that the nodes are a finite, non-zero distance apart
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    // Normalize the unit vector
    unit_vector /= distance_between_nodes;

    double radius_of_cell_one = p_node_a->GetRadius();
    double radius_of_cell_two = p_node_b->GetRadius();

    // Compute the force vector
    c_vector<double, DIM> force_between_nodes = GetMagnitudeOfForce(distance_between_nodes,radius_of_cell_one,radius_of_cell_two) * unit_vector;

    return force_between_nodes;
}

template<unsigned DIM>
double BuskeAdhesiveForce<DIM>::GetMagnitudeOfForce(double distanceBetweenNodes, double radiusOfCellOne, double radiusOfCellTwo)
{
    double dWAdd = 0.0;

    // If the cells are close enough...
    if (distanceBetweenNodes < radiusOfCellOne + radiusOfCellTwo)
    {
        // ...calculate the force contribution from their adhesive interaction energy
        double xij = 0.5*(radiusOfCellOne*radiusOfCellOne - radiusOfCellTwo*radiusOfCellTwo + distanceBetweenNodes*distanceBetweenNodes)/distanceBetweenNodes;
        double dxijdd = 1.0 - xij/distanceBetweenNodes;
        dWAdd = 2.0*mAdhesionEnergyParameter*M_PI*xij*dxijdd;
    }

    return dWAdd;
}

template<unsigned DIM>
void BuskeAdhesiveForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AdhesionEnergyParameter>" << mAdhesionEnergyParameter << "</AdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class BuskeAdhesiveForce<1>;
template class BuskeAdhesiveForce<2>;
template class BuskeAdhesiveForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeAdhesiveForce)
