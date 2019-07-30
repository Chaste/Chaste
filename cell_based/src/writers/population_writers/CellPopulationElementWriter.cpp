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

#include "CellPopulationElementWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPopulationElementWriter<ELEMENT_DIM, SPACE_DIM>::CellPopulationElementWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("results.vizelements")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationElementWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    for (typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetElementIteratorBegin();
         elem_iter != static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetElementIteratorEnd();
         ++elem_iter)
    {
        bool element_contains_dead_cells_or_deleted_nodes = false;

        // Hack that covers the case where the element contains a node that is associated with a cell that has just been killed (#1129)
        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
        {
            unsigned node_index = elem_iter->GetNodeGlobalIndex(i);

            if (pCellPopulation->GetNode(node_index)->IsDeleted())
            {
                element_contains_dead_cells_or_deleted_nodes = true;
                break;
            }
            else if (pCellPopulation->IsCellAttachedToLocationIndex(node_index))
            {
                if (pCellPopulation->GetCellUsingLocationIndex(node_index)->IsDead())
                {
                    element_contains_dead_cells_or_deleted_nodes = true;
                    break;
                }
            }
        }
        if (!element_contains_dead_cells_or_deleted_nodes)
        {
            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                *this->mpOutStream << elem_iter->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationElementWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationElementWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationElementWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationElementWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationElementWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // Loop over cells and find associated elements so in the same order as the cells in output files
    for (typename AbstractCellPopulation<SPACE_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Hack that covers the case where the element is associated with a cell that has just been killed (#1129)
        bool elem_corresponds_to_dead_cell = false;

        if (pCellPopulation->IsCellAttachedToLocationIndex(elem_index))
        {
            elem_corresponds_to_dead_cell = pCellPopulation->GetCellUsingLocationIndex(elem_index)->IsDead();
        }

        // Write element data to file
        if (!(pCellPopulation->GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
            PottsElement<SPACE_DIM>* p_element = pCellPopulation->rGetMesh().GetElement(elem_index);
            unsigned num_nodes_in_element = p_element->GetNumNodes();

            // First write the number of Nodes belonging to this VertexElement
            *this->mpOutStream << num_nodes_in_element << " ";

            // Then write the global index of each Node in this element
            for (unsigned i=0; i<num_nodes_in_element; i++)
            {
                *this->mpOutStream << p_element->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationElementWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // Loop over cells and find associated elements so in the same order as the cells in output files
    for (typename AbstractCellPopulation<SPACE_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Hack that covers the case where the element is associated with a cell that has just been killed (#1129)
        bool elem_corresponds_to_dead_cell = false;

        if (pCellPopulation->IsCellAttachedToLocationIndex(elem_index))
        {
            elem_corresponds_to_dead_cell = pCellPopulation->GetCellUsingLocationIndex(elem_index)->IsDead();
        }

        // Write element data to file
        if (!(pCellPopulation->GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
            VertexElement<SPACE_DIM, SPACE_DIM>* p_element = pCellPopulation->rGetMesh().GetElement(elem_index);
            unsigned num_nodes_in_element = p_element->GetNumNodes();

            // First write the number of Nodes belonging to this VertexElement
            *this->mpOutStream << num_nodes_in_element << " ";

            // Then write the global index of each Node in this element
            for (unsigned i=0; i<num_nodes_in_element; i++)
            {
                *this->mpOutStream << p_element->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
}

// Explicit instantiation
template class CellPopulationElementWriter<1,1>;
template class CellPopulationElementWriter<1,2>;
template class CellPopulationElementWriter<2,2>;
template class CellPopulationElementWriter<1,3>;
template class CellPopulationElementWriter<2,3>;
template class CellPopulationElementWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPopulationElementWriter)
