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

#include "VolumeConstraintPottsUpdateRule.hpp"

template<unsigned DIM>
VolumeConstraintPottsUpdateRule<DIM>::VolumeConstraintPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mDeformationEnergyParameter(0.5), // Educated guess
      mMatureCellTargetVolume(16.0) // Defaults to a 4*4 cell size
{
        /// \todo Default values don't apply in 3D.
}

template<unsigned DIM>
VolumeConstraintPottsUpdateRule<DIM>::~VolumeConstraintPottsUpdateRule()
{
}

template<unsigned DIM>
double VolumeConstraintPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    double delta_H = 0.0;

    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();

    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    // Every node must each be in at most one element
    assert(new_location_containing_elements.size() < 2);

    if (!current_node_contained && !target_node_contained)
    {
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }

    if (current_node_contained && target_node_contained)
    {
        if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }

    if (current_node_contained) // current node is in an element
    {
        unsigned current_element = (*containing_elements.begin());
        double current_volume = rCellPopulation.rGetMesh().GetVolumeOfElement(current_element);
        double current_volume_difference = current_volume - mMatureCellTargetVolume;

        delta_H += mDeformationEnergyParameter*((current_volume_difference + 1.0)*(current_volume_difference + 1.0) - current_volume_difference*current_volume_difference);
    }
    if (target_node_contained) // target node is in an element
    {
        unsigned target_element = (*new_location_containing_elements.begin());
        double target_volume = rCellPopulation.rGetMesh().GetVolumeOfElement(target_element);
        double target_volume_difference = target_volume - mMatureCellTargetVolume;

        delta_H += mDeformationEnergyParameter*((target_volume_difference - 1.0)*(target_volume_difference - 1.0) - target_volume_difference*target_volume_difference);
    }

    return delta_H;
}

template<unsigned DIM>
double VolumeConstraintPottsUpdateRule<DIM>::GetDeformationEnergyParameter()
{
    return mDeformationEnergyParameter;
}

template<unsigned DIM>
void VolumeConstraintPottsUpdateRule<DIM>::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
double VolumeConstraintPottsUpdateRule<DIM>::GetMatureCellTargetVolume() const
{
    return mMatureCellTargetVolume;
}

template<unsigned DIM>
void VolumeConstraintPottsUpdateRule<DIM>::SetMatureCellTargetVolume(double matureCellTargetVolume)
{
    assert(matureCellTargetVolume >= 0.0);
    mMatureCellTargetVolume = matureCellTargetVolume;
}

template<unsigned DIM>
void VolumeConstraintPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DeformationEnergyParameter>" << mDeformationEnergyParameter << "</DeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<MatureCellTargetVolume>" << mMatureCellTargetVolume << "</MatureCellTargetVolume>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class VolumeConstraintPottsUpdateRule<1>;
template class VolumeConstraintPottsUpdateRule<2>;
template class VolumeConstraintPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeConstraintPottsUpdateRule)
