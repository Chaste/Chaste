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

#include "ExclusionCaBasedDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned SPACE_DIM>
bool ExclusionCaBasedDivisionRule<SPACE_DIM>::IsRoomToDivide(CellPtr pParentCell, CaBasedCellPopulation<SPACE_DIM>& rCellPopulation)
{
    bool is_room = false;

    // Get node index corresponding to this cell
    unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);

    // Get the set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<SPACE_DIM>*>(&(rCellPopulation.rGetMesh()))->GetMooreNeighbouringNodeIndices(node_index);

    // Iterate through the neighbours to see if there are any available sites
    for (std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
         neighbour_iter != neighbouring_node_indices.end();
         ++neighbour_iter)
    {
        if (rCellPopulation.IsSiteAvailable(*neighbour_iter, pParentCell))
        {
            is_room = true;
            break;
        }
    }

    return is_room;
}

template<unsigned SPACE_DIM>
unsigned ExclusionCaBasedDivisionRule<SPACE_DIM>::CalculateDaughterNodeIndex(CellPtr pNewCell, CellPtr pParentCell, CaBasedCellPopulation<SPACE_DIM>& rCellPopulation)
{
    if (!IsRoomToDivide(pParentCell,rCellPopulation))
    {
        EXCEPTION("Trying to divide when there is no room to divide, check your division rule");
    }

    // Get node index corresponding to the parent cell
    unsigned parent_node_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);

    PottsMesh<SPACE_DIM>* static_cast_mesh = static_cast<PottsMesh<SPACE_DIM>*>(&(rCellPopulation.rGetMesh()));

    // Get the set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices = static_cast_mesh->GetMooreNeighbouringNodeIndices(parent_node_index);
    unsigned num_neighbours = neighbouring_node_indices.size();

    // Each node must have at least one neighbour
    assert(!neighbouring_node_indices.empty());

    std::vector<double> neighbouring_node_propensities;
    std::vector<unsigned> neighbouring_node_indices_vector;

    double total_propensity = 0.0;

    // Select neighbour at random
    for (std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
         neighbour_iter != neighbouring_node_indices.end();
         ++neighbour_iter)
    {
        neighbouring_node_indices_vector.push_back(*neighbour_iter);

        double propensity_dividing_into_neighbour = rCellPopulation.EvaluateDivisionPropensity(parent_node_index,*neighbour_iter,pParentCell);

        if (!(rCellPopulation.IsSiteAvailable(*neighbour_iter, pParentCell)))
        {
            propensity_dividing_into_neighbour = 0.0;
        }

        neighbouring_node_propensities.push_back(propensity_dividing_into_neighbour);
        total_propensity += propensity_dividing_into_neighbour;
    }
    assert(total_propensity > 0); // If this trips the cell can't divide, so we need to include this in the IsSiteAvailable() method

    for (unsigned i=0; i<num_neighbours; i++)
    {
        neighbouring_node_propensities[i] /= total_propensity;
    }

    // Sample random number to specify which move to make
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double random_number = p_gen->ranf();

    double total_probability = 0.0;
    unsigned daughter_node_index = UNSIGNED_UNSET;

    unsigned counter;
    for (counter=0; counter < num_neighbours; counter++)
    {
        total_probability += neighbouring_node_propensities[counter];
        if (total_probability >= random_number)
        {
            // Divide the parent cell to this neighbour location
            daughter_node_index = neighbouring_node_indices_vector[counter];
            break;
        }
    }
    // This loop should always break as sum(neighbouring_node_propensities) = 1

    assert(daughter_node_index != UNSIGNED_UNSET);
    assert(daughter_node_index < static_cast_mesh->GetNumNodes());

    return daughter_node_index;
}

// Explicit instantiation
template class ExclusionCaBasedDivisionRule<1>;
template class ExclusionCaBasedDivisionRule<2>;
template class ExclusionCaBasedDivisionRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExclusionCaBasedDivisionRule)
