/*

Copyright (C) University of Oxford, 2005-2010

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

#include "DiffusionUpdateRule.hpp"

template<unsigned DIM>
DiffusionUpdateRule<DIM>::DiffusionUpdateRule(double diffusionConstant)
    : AbstractUpdateRule<DIM>(),
      mDiffusionConstant(diffusionConstant)
{
}

template<unsigned DIM>
DiffusionUpdateRule<DIM>::~DiffusionUpdateRule()
{
}

template<unsigned DIM>
unsigned DiffusionUpdateRule<DIM>::GetNewLocationOfCell(unsigned currentLocationIndex,
                                                        LatticeBasedTissue<DIM>& rTissue,
                                                        double dt)
{
    assert(DIM == 2); // this method only works in 2D at present

    // Make sure we have a cell at this node
    if (rTissue.IsEmptySite(currentLocationIndex))
    {
        EXCEPTION("There is no cell at the current location.");
    }

    unsigned new_location_index = currentLocationIndex;

    double probability_of_moving = dt*mDiffusionConstant;

    if (RandomNumberGenerator::Instance()->ranf() < probability_of_moving)
    {
        // Find a random site from all of the available neighbouring nodes to move into
        std::set<unsigned> free_neighbour_indices = rTissue.GetFreeNeighbouringNodeIndices(currentLocationIndex);

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
double DiffusionUpdateRule<DIM>::GetDiffusionConstant()
{
    return mDiffusionConstant;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DiffusionUpdateRule<1>;
template class DiffusionUpdateRule<2>;
template class DiffusionUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionUpdateRule)
