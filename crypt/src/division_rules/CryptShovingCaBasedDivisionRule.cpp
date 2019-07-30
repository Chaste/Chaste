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

#include "CryptShovingCaBasedDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"

bool CryptShovingCaBasedDivisionRule::IsNodeOnBase(unsigned NodeIndex, PottsMesh<2>* pPottsMesh)
{
    std::set<unsigned> neighbouring_node_indices = pPottsMesh->GetVonNeumannNeighbouringNodeIndices(NodeIndex);
    unsigned num_neighbours = neighbouring_node_indices.size();

    // No strange neighbourhoods and in 2D so need 3 or 4 neighbours
    if (num_neighbours == 4)
    {
        return true;
    }
    else if (num_neighbours == 3)
    {
        // Quick and dirty check to see if cells are in bottom or top half of the domain
        if (NodeIndex < 0.5*pPottsMesh->GetNumNodes())
        {
            return false;
        }
        else
        {
            EXCEPTION("Cells reaching the top of the crypt need to increase length to at least double the sloughing height.");
        }
    }
    else
    {
        // If here then have <2 or >4 neighbours and not possible for 2d periodic crypt
        NEVER_REACHED;
    }
}

bool CryptShovingCaBasedDivisionRule::IsRoomToDivide(CellPtr pParentCell, CaBasedCellPopulation<2>& rCellPopulation)
{
    return true;
}

unsigned CryptShovingCaBasedDivisionRule::CalculateDaughterNodeIndex(CellPtr pNewCell,
    CellPtr pParentCell,
    CaBasedCellPopulation<2>& rCellPopulation)
{
    // Get node index corresponding to the parent cell
    unsigned parent_node_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);

    PottsMesh<2>* static_cast_mesh = static_cast<PottsMesh<2>*>(&(rCellPopulation.rGetMesh()));

    // This tracks if the cell is on the base of the crypt, and offsets the neighbours accordingly
    bool is_not_on_base = IsNodeOnBase(parent_node_index,static_cast_mesh);

    /*
     * Select a neighbour at random.
     * Sample random number to specify which move to make either 1 (E) 2 (W) or 3 (N)
     * This is as they are ordered in node index and that moves from south west to north east.
     */
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    unsigned direction = p_gen->randMod(3)+ (unsigned) is_not_on_base;

    // Stem Cells only divide vertically
    if (pParentCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        direction = 2;
    }


    std::set<unsigned> neighbouring_node_indices = static_cast_mesh->GetVonNeumannNeighbouringNodeIndices(parent_node_index);

    std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
    for (unsigned  i=0; i<direction; i++)
    {
        ++neighbour_iter;
    }
    assert(neighbour_iter != neighbouring_node_indices.end());

    unsigned daughter_node_index = *neighbour_iter;

    assert(daughter_node_index < static_cast_mesh->GetNumNodes());

    // If daughter node is occupied then move the cell north (which is always the last one in the set of neighbours)
    if (!(rCellPopulation.IsSiteAvailable(daughter_node_index, pNewCell)))
    {
        std::list<std::pair<unsigned,unsigned> > cell_moves;

        bool is_neighbour_occupied = true;

        unsigned current_node_index = parent_node_index;
        unsigned target_node_index = daughter_node_index;
        while (is_neighbour_occupied)
        {
            current_node_index = target_node_index;

            std::set<unsigned> neighbouring_node_indices = static_cast_mesh->GetVonNeumannNeighbouringNodeIndices(current_node_index);
            unsigned num_neighbours = neighbouring_node_indices.size();

            // Check to see if the current node is on the boundary
            IsNodeOnBase(current_node_index, static_cast_mesh);

            // Select the appropriate neighbour
            std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
            for (unsigned i=0; i<num_neighbours-1; i++)
            {
                ++neighbour_iter;
            }
            assert(neighbour_iter != neighbouring_node_indices.end());

            target_node_index = *neighbour_iter;

            std::pair<unsigned, unsigned> new_move(current_node_index, target_node_index);

            cell_moves.push_back(new_move);

            // If target node is unoccupied move the cell on the current node to the target node and stop shoving cells
            if (rCellPopulation.IsSiteAvailable(target_node_index, pNewCell))
            {
                is_neighbour_occupied = false;
            }

            // If target node is occupied then keep shoving the cells out of the way
            current_node_index = target_node_index;
        }

        // Do moves to free up the daughter node index
        for (std::list<std::pair<unsigned, unsigned> >::reverse_iterator reverse_iter = cell_moves.rbegin();
             reverse_iter != cell_moves.rend();
             ++reverse_iter)
        {
            assert(rCellPopulation.IsSiteAvailable(reverse_iter->second, pNewCell));
            assert(!(rCellPopulation.IsSiteAvailable(reverse_iter->first, pNewCell)));

            // Move cell from first() to second()
            rCellPopulation.MoveCellInLocationMap(rCellPopulation.GetCellUsingLocationIndex(reverse_iter->first), reverse_iter->first, reverse_iter->second);
        }

        // Check daughter site is now free
        assert(rCellPopulation.IsSiteAvailable(daughter_node_index, pNewCell));

    }
    return daughter_node_index;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptShovingCaBasedDivisionRule)
