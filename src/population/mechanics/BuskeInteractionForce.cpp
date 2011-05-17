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

#include "BuskeInteractionForce.hpp"

template<unsigned DIM>
BuskeInteractionForce<DIM>::BuskeInteractionForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mAdhesionEnergyParameter(200),        // Denoted by epsilon in Buske et al (2011) (doi:10.1371/journal.pcbi.1001045).
     mDeformationEnergyParameter(4.0/3.0), // Denoted by D in Buske et al (2011) (doi:10.1371/journal.pcbi.1001045).
     mCompressionEnergyParameter(1.0)      // Denoted by K in Buske et al (2011) (doi:10.1371/journal.pcbi.1001045).

{
}

template<unsigned DIM>
c_vector<double, DIM> BuskeInteractionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                             unsigned nodeBGlobalIndex,
                                                                             AbstractCellPopulation<DIM>& rCellPopulation)
{
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

    // Determine cell radii
//    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
//    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double radius_of_cell_one = 1.0;
    double radius_of_cell_two = 1.0;

    // Compute the force vector
    c_vector<double, DIM> force_between_nodes = GetMagnitudeOfForce(distance_between_nodes,radius_of_cell_one,radius_of_cell_two) * unit_vector;

    return force_between_nodes;
}

template<unsigned DIM>
double BuskeInteractionForce<DIM>::GetMagnitudeOfForce(double distanceBetweenNodes, double radiusOfCellOne, double radiusOfCellTwo)
{
	double xij = 0.5*(radiusOfCellOne*radiusOfCellOne - radiusOfCellTwo*radiusOfCellTwo + distanceBetweenNodes*distanceBetweenNodes)/distanceBetweenNodes;

	double dxijdd = 1.0 - 2.0*xij/distanceBetweenNodes;

	// Calculate contribution from adhesive interaction energy
	double dWAdd = -2.0*mAdhesionEnergyParameter*M_PI*xij*dxijdd;

	// Calculate contribution from deformation interaction energy
	double dWDdd;

	if (distanceBetweenNodes < radiusOfCellOne + radiusOfCellTwo)
	{
		dWDdd = -pow(radiusOfCellOne + radiusOfCellTwo - distanceBetweenNodes,1.5)
                *pow(radiusOfCellOne*radiusOfCellTwo/(radiusOfCellOne+radiusOfCellTwo),0.5)
			    /mDeformationEnergyParameter;
	}
	else  // no deformation energy contribution as too far apart
	{
		dWDdd = 0.0;
	}

	// Calculate contribution from compression interaction energy Still todo
	double dWKdd = 0.0;

	return dWAdd + dWDdd + dWKdd;
}

template<unsigned DIM>
void BuskeInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	// Call method on direct parent class
	AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BuskeInteractionForce<1>;
template class BuskeInteractionForce<2>;
template class BuskeInteractionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeInteractionForce)
