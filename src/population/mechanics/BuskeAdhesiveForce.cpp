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

#include "BuskeAdhesiveForce.hpp"

template<unsigned DIM>
BuskeAdhesiveForce<DIM>::BuskeAdhesiveForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mAdhesionEnergyParameter(200)        // Denoted by epsilon in Buske et al (2011) (doi:10.1371/journal.pcbi.1001045).
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
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL);

    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = rCellPopulation.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = rCellPopulation.GetNode(nodeBGlobalIndex)->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
    c_vector<double, DIM> unit_vector = node_b_location - node_a_location;

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

    // Get the cell radii
    NodesOnlyMesh<DIM>& r_mesh = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh();
    double radius_of_cell_one = r_mesh.GetCellRadius(nodeAGlobalIndex);
    double radius_of_cell_two = r_mesh.GetCellRadius(nodeBGlobalIndex);

    // Compute the force vector
    c_vector<double, DIM> force_between_nodes = GetMagnitudeOfForce(distance_between_nodes,radius_of_cell_one,radius_of_cell_two) * unit_vector;

    return force_between_nodes;
}

template<unsigned DIM>
double BuskeAdhesiveForce<DIM>::GetMagnitudeOfForce(double distanceBetweenNodes, double radiusOfCellOne, double radiusOfCellTwo)
{
    double xij = 0.5*(radiusOfCellOne*radiusOfCellOne - radiusOfCellTwo*radiusOfCellTwo + distanceBetweenNodes*distanceBetweenNodes)/distanceBetweenNodes;

    double dxijdd = 1.0 - xij/distanceBetweenNodes;

    // Calculate contribution from adhesive interaction energy
    double dWAdd = 2.0*mAdhesionEnergyParameter*M_PI*xij*dxijdd;

    return dWAdd; //
}

template<unsigned DIM>
void BuskeAdhesiveForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AdhesionEnergyParameter>" << mAdhesionEnergyParameter << "</AdhesionEnergyParameter> \n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BuskeAdhesiveForce<1>;
template class BuskeAdhesiveForce<2>;
template class BuskeAdhesiveForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeAdhesiveForce)
