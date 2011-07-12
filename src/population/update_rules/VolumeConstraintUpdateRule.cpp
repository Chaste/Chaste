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


#include "VolumeConstraintUpdateRule.hpp"

template<unsigned DIM>
VolumeConstraintUpdateRule<DIM>::VolumeConstraintUpdateRule()
    : AbstractPottsUpdateRule<DIM>()
{
}

template<unsigned DIM>
VolumeConstraintUpdateRule<DIM>::~VolumeConstraintUpdateRule()
{
}
template<unsigned DIM>
void VolumeConstraintUpdateRule<DIM>::EvaluateHamiltonianContribution(double& DeltaH, unsigned CurrentNodeIndex, unsigned TargetNodeIndex,
																	  AbstractCellPopulation<2>& rCellPopulation)
{
	// Make sure that we are in the correct dimension - this code will be eliminated at compile time
	assert(DIM == 2); // this method only works in 2D at present

	// Throw an exception message if not using a PottsBasedCellPopulation
	/** TODO this is probably not the best way of doing this it will slow things down probably should pass in a PottsBasedCellPopulation #1665 */
	if (dynamic_cast<PottsBasedCellPopulation*>(&rCellPopulation) == NULL)
	{
		EXCEPTION("VolumeConstraintUpdateRule is to be used with a PottsBasedCellPopulation only");
	}

	// Helper variable that is a static cast of the cell population
	PottsBasedCellPopulation* p_cell_population = static_cast<PottsBasedCellPopulation*>(&rCellPopulation);

	std::set<unsigned> containing_elements = p_cell_population->GetNode(CurrentNodeIndex)->rGetContainingElementIndices();
	std::set<unsigned> new_location_containing_elements = p_cell_population->GetNode(TargetNodeIndex)->rGetContainingElementIndices();

	// All nodes should be in at most one element.
	assert(new_location_containing_elements.size() <= 1);

	// Both elements should be different to use this method
	assert(new_location_containing_elements.begin() != containing_elements.begin());
	assert((containing_elements.size()>0) || (new_location_containing_elements.size()>0));

	double delta_H = 0.0; // H_1-H_0 differnece in Hamiltonian after and before swap.

	//// VOLUME CONSTRAINT UPDATE RULE ////
	double lambda_volume = 0.5;
	double target_volume = 16.0;

	if (containing_elements.size() == 1) // current node is in an element
	{
		unsigned current_element = (*containing_elements.begin());
		delta_H += lambda_volume*(pow(p_cell_population->rGetMesh().GetVolumeOfElement(current_element) + 1.0 - target_volume, 2.0) - pow(p_cell_population->rGetMesh().GetVolumeOfElement(current_element) - target_volume, 2.0));
	}
	if (new_location_containing_elements.size() == 1) // target node is in an element
	{
		unsigned target_element = (*new_location_containing_elements.begin());
		delta_H += lambda_volume*(pow(p_cell_population->rGetMesh().GetVolumeOfElement(target_element) - 1.0 - target_volume, 2.0) - pow(p_cell_population->rGetMesh().GetVolumeOfElement(target_element) - target_volume, 2.0));
	}


	//// CELL CELL ADHESION UPDATE RULE ////
	double lambda_contact = 0.1;


	// Iterate over nodes neighbouring the target node to work out the contact energy contribution
	std::set<unsigned> target_neighboring_node_indices = p_cell_population->rGetMesh().GetNeighbouringNodeIndices(TargetNodeIndex);

	for (std::set<unsigned>::iterator iter = target_neighboring_node_indices.begin();
		 iter != target_neighboring_node_indices.end();
		 ++iter)
	{
		std::set<unsigned> neighboring_node_containing_elements = p_cell_population->rGetMesh().GetNode(*iter)->rGetContainingElementIndices();


		if ( neighboring_node_containing_elements.size() == 1u )
		{
			unsigned neighbour_element = (*neighboring_node_containing_elements.begin());

			if (new_location_containing_elements.size() == 1u) // target node is in an element
			{
				unsigned target_element = (*new_location_containing_elements.begin());
				// If the nodes are currently from different elements
				if ( target_element != neighbour_element )
				{
					delta_H -= lambda_contact;
				}
			}
			else
			{
				delta_H -= lambda_contact;
			}
			if (containing_elements.size() == 1u) // current node is in an element
			{
				unsigned current_element = (*containing_elements.begin());
				// If the nodes will be in different elements after swap
				if ( current_element != neighbour_element )
				{
					delta_H += lambda_contact;
				}
			}
			else
			{
				delta_H += lambda_contact;
			}
		}
	}
}


template<unsigned DIM>
void VolumeConstraintUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VolumeConstraintUpdateRule<1>;
template class VolumeConstraintUpdateRule<2>;
template class VolumeConstraintUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeConstraintUpdateRule)
