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
    : AbstractPottsUpdateRule<DIM>(),
      mDeformationEnergyParameter(0.5), // Educated guess
      mMatureCellTargetVolume(16.0) // Defaults to a 4*4 cell size
{
}

template<unsigned DIM>
VolumeConstraintUpdateRule<DIM>::~VolumeConstraintUpdateRule()
{
}

template<unsigned DIM>
double VolumeConstraintUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        AbstractCellPopulation<2>& rCellPopulation)
{
	double delta_H = 0.0;

	// Make sure that we are in the correct dimension - this code will be eliminated at compile time
	assert(DIM == 2); // this method only works in 2D at present

	// Throw an exception message if not using a PottsBasedCellPopulation
	///\todo this is probably not the best way of doing this it will slow things down probably should pass in a PottsBasedCellPopulation (#1665)
	if (dynamic_cast<PottsBasedCellPopulation*>(&rCellPopulation) == NULL)
	{
		EXCEPTION("VolumeConstraintUpdateRule is to be used with a PottsBasedCellPopulation only");
	}

	// Helper variable that is a static cast of the cell population
	PottsBasedCellPopulation* p_cell_population = static_cast<PottsBasedCellPopulation*>(&rCellPopulation);

	std::set<unsigned> containing_elements = p_cell_population->GetNode(currentNodeIndex)->rGetContainingElementIndices();
	std::set<unsigned> new_location_containing_elements = p_cell_population->GetNode(targetNodeIndex)->rGetContainingElementIndices();

	// All nodes should be in at most one element.
	assert(new_location_containing_elements.size() <= 1);

	// Both elements should be different to use this method
	assert(new_location_containing_elements.begin() != containing_elements.begin());
	assert((containing_elements.size()>0) || (new_location_containing_elements.size()>0));

	if (containing_elements.size() == 1) // current node is in an element
	{
		unsigned current_element = (*containing_elements.begin());
		delta_H += mDeformationEnergyParameter*(pow(p_cell_population->rGetMesh().GetVolumeOfElement(current_element) + 1.0 - mMatureCellTargetVolume, 2.0) - pow(p_cell_population->rGetMesh().GetVolumeOfElement(current_element) - mMatureCellTargetVolume, 2.0));
	}
	if (new_location_containing_elements.size() == 1) // target node is in an element
	{
		unsigned target_element = (*new_location_containing_elements.begin());
		delta_H += mDeformationEnergyParameter*(pow(p_cell_population->rGetMesh().GetVolumeOfElement(target_element) - 1.0 - mMatureCellTargetVolume, 2.0) - pow(p_cell_population->rGetMesh().GetVolumeOfElement(target_element) - mMatureCellTargetVolume, 2.0));
	}

	return delta_H;
}

template<unsigned DIM>
double VolumeConstraintUpdateRule<DIM>::GetDeformationEnergyParameter()
{
    return mDeformationEnergyParameter;
}

template<unsigned DIM>
void VolumeConstraintUpdateRule<DIM>::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
double VolumeConstraintUpdateRule<DIM>::GetMatureCellTargetVolume() const
{
    return mMatureCellTargetVolume;
}

template<unsigned DIM>
void VolumeConstraintUpdateRule<DIM>::SetMatureCellTargetVolume(double matureCellTargetVolume)
{
    assert(matureCellTargetVolume >= 0.0);
    mMatureCellTargetVolume = matureCellTargetVolume;
}

template<unsigned DIM>
void VolumeConstraintUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DeformationEnergyParameter>" << mDeformationEnergyParameter << "</DeformationEnergyParameter> \n";
	*rParamsFile << "\t\t\t<MatureCellTargetVolume>" << mMatureCellTargetVolume << "</MatureCellTargetVolume> \n";

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
