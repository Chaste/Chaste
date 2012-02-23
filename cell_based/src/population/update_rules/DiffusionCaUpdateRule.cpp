/*

Copyright (c) 2005-2012, University of Oxford.
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
