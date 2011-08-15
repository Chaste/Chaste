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

#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

template<unsigned DIM>
SurfaceAreaConstraintPottsUpdateRule<DIM>::SurfaceAreaConstraintPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mDeformationEnergyParameter(0.5), // Educated guess
      mMatureCellTargetSurfaceArea(16.0) // Defaults to a 4*4 cell size
{
    /// \todo Default values don't apply in 3D. 
}

template<unsigned DIM>
SurfaceAreaConstraintPottsUpdateRule<DIM>::~SurfaceAreaConstraintPottsUpdateRule()
{
}

template<unsigned DIM>
double SurfaceAreaConstraintPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        PottsBasedCellPopulation& rCellPopulation)
{
	double delta_H = 0.0;

	// This method only works in 2D at present
	assert(DIM == 2);

	std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
	std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();

    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    // Every node must each be in at most one element.
    assert(new_location_containing_elements.size() < 2);

    if(!current_node_contained && !target_node_contained)
    {
    	EXCEPTION("At least one of the current node or target node must be in an element.");
    }

    if (current_node_contained && target_node_contained)
    {
        if(*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }

    // Iterate over nodes neighbouring the target node to work out the change in surface area
	unsigned neighbours_in_same_element_as_current_node = 0;
	unsigned neighbours_in_same_element_as_target_node = 0;
	std::set<unsigned> target_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(targetNodeIndex);
	for (std::set<unsigned>::iterator iter = target_neighbouring_node_indices.begin();
		 iter != target_neighbouring_node_indices.end();
		 ++iter)
	{
		std::set<unsigned> neighbouring_node_containing_elements = rCellPopulation.rGetMesh().GetNode(*iter)->rGetContainingElementIndices();

		// Every node must each be in at most one element
		assert(neighbouring_node_containing_elements.size() < 2);

		bool neighbouring_node_contained = !neighbouring_node_containing_elements.empty();

		if (neighbouring_node_contained && target_node_contained)
		{
			unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
			unsigned target_element = (*new_location_containing_elements.begin());
			if (target_element == neighbour_element)
			{
				neighbours_in_same_element_as_target_node++;
			}
		}
		if (neighbouring_node_contained && current_node_contained)
		{
			unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
			unsigned current_element = (*containing_elements.begin());
			if (current_element == neighbour_element)
			{
				neighbours_in_same_element_as_current_node++;
			}
		}
	}
	assert(neighbours_in_same_element_as_current_node<5);
	assert(neighbours_in_same_element_as_target_node<5);

	double change_in_surface_area[5] = {4,2,0,-2,-4};

	if (current_node_contained) // current node is in an element
	{
		unsigned current_element = (*containing_elements.begin());
		double current_surface_area = rCellPopulation.rGetMesh().GetSurfaceAreaOfElement(current_element);
		double current_surface_area_difference = current_surface_area - mMatureCellTargetSurfaceArea;
		double current_surface_area_difference_after_switch = current_surface_area_difference + change_in_surface_area[neighbours_in_same_element_as_current_node];

		delta_H += mDeformationEnergyParameter*(current_surface_area_difference_after_switch*current_surface_area_difference_after_switch - current_surface_area_difference*current_surface_area_difference);
	}
	if (target_node_contained) // target node is in an element
	{
		unsigned target_element = (*new_location_containing_elements.begin());
        double target_surface_area = rCellPopulation.rGetMesh().GetSurfaceAreaOfElement(target_element);
        double target_surface_area_difference = target_surface_area - mMatureCellTargetSurfaceArea;
        double target_surface_area_difference_after_switch = target_surface_area_difference - change_in_surface_area[neighbours_in_same_element_as_target_node];

		delta_H += mDeformationEnergyParameter*(target_surface_area_difference_after_switch*target_surface_area_difference_after_switch - target_surface_area_difference*target_surface_area_difference);
	}

	return delta_H;
}

template<unsigned DIM>
double SurfaceAreaConstraintPottsUpdateRule<DIM>::GetDeformationEnergyParameter()
{
    return mDeformationEnergyParameter;
}

template<unsigned DIM>
void SurfaceAreaConstraintPottsUpdateRule<DIM>::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
double SurfaceAreaConstraintPottsUpdateRule<DIM>::GetMatureCellTargetSurfaceArea() const
{
    return mMatureCellTargetSurfaceArea;
}

template<unsigned DIM>
void SurfaceAreaConstraintPottsUpdateRule<DIM>::SetMatureCellTargetSurfaceArea(double matureCellTargetSurfaceArea)
{
    assert(matureCellTargetSurfaceArea >= 0.0);
    mMatureCellTargetSurfaceArea = matureCellTargetSurfaceArea;
}

template<unsigned DIM>
void SurfaceAreaConstraintPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DeformationEnergyParameter>" << mDeformationEnergyParameter << "</DeformationEnergyParameter> \n";
	*rParamsFile << "\t\t\t<MatureCellTargetSurfaceArea>" << mMatureCellTargetSurfaceArea << "</MatureCellTargetSurfaceArea> \n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SurfaceAreaConstraintPottsUpdateRule<1>;
template class SurfaceAreaConstraintPottsUpdateRule<2>;
template class SurfaceAreaConstraintPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceAreaConstraintPottsUpdateRule)
