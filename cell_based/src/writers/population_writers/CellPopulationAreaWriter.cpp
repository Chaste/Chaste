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

#include "CellPopulationAreaWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "ApoptoticCellProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPopulationAreaWriter<ELEMENT_DIM, SPACE_DIM>::CellPopulationAreaWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("cellpopulationareas.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM==2 || SPACE_DIM==3); // LCOV_EXCL_LINE

    VertexMesh<ELEMENT_DIM, SPACE_DIM>* voronoi_tessellation = pCellPopulation->GetVoronoiTessellation();

    assert (voronoi_tessellation != nullptr);

    // Don't use the Voronoi tessellation to calculate the total area of the mesh because it gives huge areas for boundary cells
    double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetVolume();
    double apoptotic_area = 0.0;

    // Loop over elements of mpVoronoiTessellation
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = voronoi_tessellation->GetElementIteratorBegin();
         elem_iter != voronoi_tessellation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in mpVoronoiTessellation
        unsigned elem_index = elem_iter->GetIndex();

        // Get the index of the corresponding node in mrMesh
        unsigned node_index = voronoi_tessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        // Discount ghost nodes
        if (!pCellPopulation->IsGhostNode(node_index))
        {
            // Get the cell corresponding to this node
            CellPtr p_cell =  pCellPopulation->GetCellUsingLocationIndex(node_index);

            // Only bother calculating the area/volume of apoptotic cells
            bool cell_is_apoptotic = p_cell->HasCellProperty<ApoptoticCellProperty>();
            if (cell_is_apoptotic)
            {
                double cell_volume = voronoi_tessellation->GetVolumeOfElement(elem_index);
                apoptotic_area += cell_volume;
            }
        }
    }
    *this->mpOutStream << total_area << " " << apoptotic_area;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationAreaWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationAreaWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationAreaWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationAreaWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationAreaWriter cannot be used with a VertexBasedCellPopulation");
}

// Explicit instantiation
template class CellPopulationAreaWriter<1,1>;
template class CellPopulationAreaWriter<1,2>;
template class CellPopulationAreaWriter<2,2>;
template class CellPopulationAreaWriter<1,3>;
template class CellPopulationAreaWriter<2,3>;
template class CellPopulationAreaWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPopulationAreaWriter)
