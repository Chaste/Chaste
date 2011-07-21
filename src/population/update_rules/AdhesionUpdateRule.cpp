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

#include "AdhesionUpdateRule.hpp"

template<unsigned DIM>
AdhesionUpdateRule<DIM>::AdhesionUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mCellCellAdhesionEnergyParameter(0.1), // Educated guess
      mCellBoundaryAdhesionEnergyParameter(0.2) // Educated guess
{
}

template<unsigned DIM>
AdhesionUpdateRule<DIM>::~AdhesionUpdateRule()
{
}

template<unsigned DIM>
double AdhesionUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                unsigned targetNodeIndex,
                                                                PottsBasedCellPopulation& rCellPopulation)
{
    // This method only works in 2D at present
	assert(DIM == 2);

	std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
	std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();

	bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    // At least one of the current node and target node must be in an element
	assert(current_node_contained || target_node_contained);

	// Every node must each be in at most one element
	assert(new_location_containing_elements.size() < 2);

	// The current node and target node must not be in the same element
	assert(new_location_containing_elements.begin() != containing_elements.begin());

	// Iterate over nodes neighbouring the target node to work out the contact energy contribution
    double delta_H = 0.0;
	std::set<unsigned> target_neighbouring_node_indices = rCellPopulation.rGetMesh().GetNeighbouringNodeIndices(targetNodeIndex);
	for (std::set<unsigned>::iterator iter = target_neighbouring_node_indices.begin();
		 iter != target_neighbouring_node_indices.end();
		 ++iter)
	{
		std::set<unsigned> neighbouring_node_containing_elements = rCellPopulation.rGetMesh().GetNode(*iter)->rGetContainingElementIndices();

        // Every node must each be in at most one element
		assert(neighbouring_node_containing_elements.size() < 2);
        bool neighbouring_node_contained = !neighbouring_node_containing_elements.empty();

		// Before move (H_0)
		if (neighbouring_node_contained && current_node_contained)
		{
			unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
			unsigned current_element = (*containing_elements.begin());

			// The nodes are currently contained in different elements
			if (current_element != neighbour_element)
			{
				delta_H += mCellCellAdhesionEnergyParameter;
			}
		}
		else if ( (neighbouring_node_contained && !current_node_contained) || (!neighbouring_node_contained && current_node_contained) )
		{
			// One node is in an element and the other is in the medium
			delta_H += mCellBoundaryAdhesionEnergyParameter;
		}

		// After move (H_1)
		if (neighbouring_node_contained && target_node_contained)
		{
			unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
			unsigned target_element = (*new_location_containing_elements.begin());

            // The nodes are currently contained in different elements
			if ( target_element != neighbour_element )
			{
				delta_H -= mCellCellAdhesionEnergyParameter;
			}
		}
		if ( (neighbouring_node_contained && !target_node_contained) || (!neighbouring_node_contained && target_node_contained) )
		{
			// One node is in an element and the other is in the medium
			delta_H -= mCellBoundaryAdhesionEnergyParameter;
		}
	}

	return delta_H;
}

template<unsigned DIM>
double AdhesionUpdateRule<DIM>::GetCellCellAdhesionEnergyParameter()
{
    return mCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double AdhesionUpdateRule<DIM>::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void AdhesionUpdateRule<DIM>::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void AdhesionUpdateRule<DIM>::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void AdhesionUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<CellCellAdhesionEnergyParameter>" << mCellCellAdhesionEnergyParameter << "</CellCellAdhesionEnergyParameter> \n";
    *rParamsFile << "\t\t\t<CellBoundaryAdhesionEnergyParameter>" << mCellBoundaryAdhesionEnergyParameter << "</CellBoundaryAdhesionEnergyParameter> \n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AdhesionUpdateRule<1>;
template class AdhesionUpdateRule<2>;
template class AdhesionUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AdhesionUpdateRule)
