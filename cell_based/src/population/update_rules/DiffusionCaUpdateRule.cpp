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

#include "DiffusionCaUpdateRule.hpp"

#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
DiffusionCaUpdateRule<DIM>::DiffusionCaUpdateRule(double diffusionConstant)
    : AbstractCaUpdateRule<DIM>(),
      mDiffusionConstant(diffusionConstant)
{
}

template<unsigned DIM>
DiffusionCaUpdateRule<DIM>::~DiffusionCaUpdateRule()
{
}

template<unsigned DIM>
unsigned DiffusionCaUpdateRule<DIM>::GetNewLocationOfCell(unsigned currentLocationIndex,
                                                        CaBasedCellPopulation<DIM>& rCellPopulation,
                                                        double dt)
{
    assert(DIM == 2); // this method only works in 2D at present

    // Make sure we have a cell at this node
    if (rCellPopulation.IsEmptySite(currentLocationIndex))
    {
        EXCEPTION("There is no cell at the current location.");
    }

    unsigned new_location_index = currentLocationIndex;

    double probability_of_moving = dt*mDiffusionConstant;

    if (RandomNumberGenerator::Instance()->ranf() < probability_of_moving)
    {
        // Find a random site from all of the available neighbouring nodes to move into
        std::set<unsigned> free_neighbour_indices = rCellPopulation.GetFreeNeighbouringNodeIndices(currentLocationIndex);

        if (!free_neighbour_indices.empty()) // only move if there is any free space
        {
            unsigned num_free_neighbours = free_neighbour_indices.size();
            unsigned chosen_neighbour = RandomNumberGenerator::Instance()->randMod(num_free_neighbours);

            std::set<unsigned>::iterator neighbour_iter = free_neighbour_indices.begin();
            for (unsigned i=0; i<chosen_neighbour; i++)
            {
                neighbour_iter++;
            }
            new_location_index = *neighbour_iter;
        }
    }

    return new_location_index;
}

template<unsigned DIM>
double DiffusionCaUpdateRule<DIM>::GetDiffusionConstant()
{
    return mDiffusionConstant;
}

template<unsigned DIM>
void DiffusionCaUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DiffusionConstant>" << mDiffusionConstant << "</DiffusionConstant>\n";

    // Call method on direct parent class
    AbstractCaUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DiffusionCaUpdateRule<1>;
template class DiffusionCaUpdateRule<2>;
template class DiffusionCaUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionCaUpdateRule)
