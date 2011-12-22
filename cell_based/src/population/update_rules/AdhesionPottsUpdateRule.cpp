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

#include "AdhesionPottsUpdateRule.hpp"

template<unsigned DIM>
AdhesionPottsUpdateRule<DIM>::AdhesionPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mCellCellAdhesionEnergyParameter(0.1), // Educated guess
      mCellBoundaryAdhesionEnergyParameter(0.2) // Educated guess
{
}

template<unsigned DIM>
AdhesionPottsUpdateRule<DIM>::~AdhesionPottsUpdateRule()
{
}

template<unsigned DIM>
double AdhesionPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                unsigned targetNodeIndex,
                                                                PottsBasedCellPopulation<DIM>& rCellPopulation)
{
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

    // Iterate over nodes neighbouring the target node to work out the contact energy contribution
    double delta_H = 0.0;
    std::set<unsigned> target_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(targetNodeIndex);
    for (std::set<unsigned>::iterator iter = target_neighbouring_node_indices.begin();
         iter != target_neighbouring_node_indices.end();
         ++iter)
    {
        std::set<unsigned> neighbouring_node_containing_elements = rCellPopulation.rGetMesh().GetNode(*iter)->rGetContainingElementIndices();

        // Every node must each be in at most one element
        assert(neighbouring_node_containing_elements.size() < 2);

        bool neighbouring_node_contained = !neighbouring_node_containing_elements.empty();

        /**
         * Before the move, we have a negative contribution (H_0) to the Hamiltonian if:
         * the target node and neighbouring node are NOT contained in the same Potts element;
         * the neighbouring node is contained in a Potts element, but the target node is not; or
         * the target node is contained in a Potts element, but the neighbouring node is not.
         */
        if (neighbouring_node_contained && target_node_contained)
        {
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            unsigned target_element = (*new_location_containing_elements.begin());
            if (target_element != neighbour_element)
            {
                // The nodes are currently contained in different elements
                delta_H -= GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(target_element), rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
            }
        }
        else if (neighbouring_node_contained && !target_node_contained)
        {
            // The neighbouring node is contained in a Potts element, but the target node is not
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            delta_H -= GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
        }
        else if (!neighbouring_node_contained && target_node_contained)
        {
            // The target node is contained in a Potts element, but the neighbouring node is not
            unsigned target_element = (*new_location_containing_elements.begin());
            delta_H -= GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(target_element));
        }

        /**
         * After the move, we have a positive contribution (H_1) to the Hamiltonian if:
         * the current node and neighbouring node are contained in different Potts elements;
         * the neighbouring node is contained in a Potts element, but the current node is not; or
         * the current node is contained in a Potts element, but the neighbouring node is not.
         */
        if (neighbouring_node_contained && current_node_contained)
        {
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            unsigned current_element = (*containing_elements.begin());
            if (current_element != neighbour_element)
            {
                // The nodes are currently contained in different elements
                delta_H += GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(current_element),rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
            }
        }
        else if (neighbouring_node_contained && !current_node_contained)
        {
            // The neighbouring node is contained in a Potts element, but the current node is not
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            delta_H += GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
        }
        else if (!neighbouring_node_contained && current_node_contained)
        {
            // The current node is contained in a Potts element, but the neighbouring node is not
            unsigned current_element = (*containing_elements.begin());
            delta_H += GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(current_element));
        }
    }

    return delta_H;
}

template<unsigned DIM>
double AdhesionPottsUpdateRule<DIM>::GetCellCellAdhesionEnergy(CellPtr pCellA, CellPtr pCellB)
{
    return GetCellCellAdhesionEnergyParameter();
}

template<unsigned DIM>
double AdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergy(CellPtr pCell)
{
    return GetCellBoundaryAdhesionEnergyParameter();
}

template<unsigned DIM>
double AdhesionPottsUpdateRule<DIM>::GetCellCellAdhesionEnergyParameter()
{
    return mCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double AdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void AdhesionPottsUpdateRule<DIM>::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void AdhesionPottsUpdateRule<DIM>::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void AdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCellAdhesionEnergyParameter>" << mCellCellAdhesionEnergyParameter << "</CellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<CellBoundaryAdhesionEnergyParameter>" << mCellBoundaryAdhesionEnergyParameter << "</CellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AdhesionPottsUpdateRule<1>;
template class AdhesionPottsUpdateRule<2>;
template class AdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AdhesionPottsUpdateRule)
