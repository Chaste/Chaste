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

#include "CellPopulationAdjacencyMatrixWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellLabel.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPopulationAdjacencyMatrixWriter<ELEMENT_DIM, SPACE_DIM>::CellPopulationAdjacencyMatrixWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("cellpopulationadjacency.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAdjacencyMatrixWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    ///\todo #2645 - if efficiency is an issue, check if this is really needed
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();

    // Store a map between cells numbered 1 to n and location indices
    std::map<unsigned,unsigned> local_cell_id_location_index_map;

    unsigned local_cell_id = 0;
    for (typename AbstractCellPopulation<SPACE_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        local_cell_id_location_index_map[pCellPopulation->GetLocationIndexUsingCell(*cell_iter)] = local_cell_id;
        local_cell_id++;
    }
    assert(local_cell_id = num_cells+1);

    // Iterate over cells and calculate the adjacency matrix (stored as a long vector)
    std::vector<double> adjacency_matrix(num_cells*num_cells);
    for (unsigned i=0; i<num_cells*num_cells; i++)
    {
        adjacency_matrix[i] = 0;
    }

    for (typename AbstractCellPopulation<SPACE_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        // Determine whether this cell is labelled
        bool cell_is_labelled = cell_iter->template HasCellProperty<CellLabel>();

        // Get the location index corresponding to this cell
        unsigned index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = pCellPopulation->GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
        if (!neighbour_indices.empty())
        {
            unsigned local_cell_index = local_cell_id_location_index_map[index];

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
                // If both cell_iter and p_neighbour_cell are not labelled, then set type_of_link to 1
                unsigned type_of_link = 1;

                // Determine whether this neighbour is labelled
                CellPtr p_neighbour_cell = pCellPopulation->GetCellUsingLocationIndex(*neighbour_iter);
                bool neighbour_is_labelled = p_neighbour_cell->template HasCellProperty<CellLabel>();

                if (cell_is_labelled != neighbour_is_labelled)
                {
                    // Here cell_iter is labelled but p_neighbour_cell is not, or vice versa, so set type_of_link to 3
                    type_of_link = 3;
                }
                else if (cell_is_labelled)
                {
                    // Here both cell_iter and p_neighbour_cell are labelled, so set type_of_link to 2
                    type_of_link = 2;
                }

                unsigned local_neighbour_index = local_cell_id_location_index_map[*neighbour_iter];
                adjacency_matrix[local_cell_index + num_cells*local_neighbour_index] = type_of_link;
                adjacency_matrix[num_cells*local_cell_index + local_neighbour_index] = type_of_link;
            }
        }
    }

    // Output the number of cells and the elements of the adjacency matrix
    *this->mpOutStream << num_cells << "\t";
    for (unsigned i=0; i<num_cells*num_cells; i++)
    {
        *this->mpOutStream << adjacency_matrix[i] << "\t";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAdjacencyMatrixWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    ///\todo #2645 - if efficiency is an issue, check if this is really needed
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();

    // Store a map between cells numbered 1 to n and location indices
    std::map<unsigned,unsigned> local_cell_id_location_index_map;

    unsigned local_cell_id = 0;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        local_cell_id_location_index_map[pCellPopulation->GetLocationIndexUsingCell(*cell_iter)] = local_cell_id;
        local_cell_id++;
    }
    assert(local_cell_id = num_cells+1);

    // Iterate over cells and calculate the adjacency matrix (stored as a long vector)
    std::vector<double> adjacency_matrix(num_cells*num_cells);
    for (unsigned i=0; i<num_cells*num_cells; i++)
    {
        adjacency_matrix[i] = 0;
    }

    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        // Determine whether this cell is labelled
        bool cell_is_labelled = cell_iter->template HasCellProperty<CellLabel>();

        // Get the location index corresponding to this cell
        unsigned index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = pCellPopulation->GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
        if (!neighbour_indices.empty())
        {
            unsigned local_cell_index = local_cell_id_location_index_map[index];

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
                // If both cell_iter and p_neighbour_cell are not labelled, then set type_of_link to 1
                unsigned type_of_link = 1;

                // Determine whether this neighbour is labelled
                CellPtr p_neighbour_cell = pCellPopulation->GetCellUsingLocationIndex(*neighbour_iter);
                bool neighbour_is_labelled = p_neighbour_cell->template HasCellProperty<CellLabel>();

                if (cell_is_labelled != neighbour_is_labelled)
                {
                    // Here cell_iter is labelled but p_neighbour_cell is not, or vice versa, so set type_of_link to 3
                    type_of_link = 3;
                }
                else if (cell_is_labelled)
                {
                    // Here both cell_iter and p_neighbour_cell are labelled, so set type_of_link to 2
                    type_of_link = 2;
                }

                unsigned local_neighbour_index = local_cell_id_location_index_map[*neighbour_iter];
                adjacency_matrix[local_cell_index + num_cells*local_neighbour_index] = type_of_link;
                adjacency_matrix[num_cells*local_cell_index + local_neighbour_index] = type_of_link;
            }
        }
    }

    // Output the number of cells and the elements of the adjacency matrix
    *this->mpOutStream << num_cells << "\t";
    for (unsigned i=0; i<num_cells*num_cells; i++)
    {
        *this->mpOutStream << adjacency_matrix[i] << "\t";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAdjacencyMatrixWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAdjacencyMatrixWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAdjacencyMatrixWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAdjacencyMatrixWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

// Explicit instantiation
template class CellPopulationAdjacencyMatrixWriter<1,1>;
template class CellPopulationAdjacencyMatrixWriter<1,2>;
template class CellPopulationAdjacencyMatrixWriter<2,2>;
template class CellPopulationAdjacencyMatrixWriter<1,3>;
template class CellPopulationAdjacencyMatrixWriter<2,3>;
template class CellPopulationAdjacencyMatrixWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPopulationAdjacencyMatrixWriter)
