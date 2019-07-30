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

#include "VoronoiDataWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VoronoiDataWriter<ELEMENT_DIM, SPACE_DIM>::VoronoiDataWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("voronoi.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VoronoiDataWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM==2 || SPACE_DIM==3); // LCOV_EXCL_LINE
    VertexMesh<ELEMENT_DIM,SPACE_DIM>* voronoi_tesselation = pCellPopulation->GetVoronoiTessellation();

    // Loop over elements of voronoi_tesselation
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = voronoi_tesselation->GetElementIteratorBegin();
         elem_iter != voronoi_tesselation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in voronoi_tesselation
        unsigned elem_index = elem_iter->GetIndex();

        // Get the index of the corresponding node in mrMesh
        unsigned node_index = voronoi_tesselation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        // Write node index and location to file
        *this->mpOutStream << node_index << " ";
        c_vector<double, SPACE_DIM> node_location = pCellPopulation->GetNode(node_index)->rGetLocation();
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << node_location[i] << " ";
        }

        double cell_volume = voronoi_tesselation->GetVolumeOfElement(elem_index);
        double cell_surface_area = voronoi_tesselation->GetSurfaceAreaOfElement(elem_index);
        *this->mpOutStream << cell_volume << " " << cell_surface_area << " ";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VoronoiDataWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VoronoiDataWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VoronoiDataWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VoronoiDataWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VoronoiDataWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VoronoiDataWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VoronoiDataWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VoronoiDataWriter cannot be used with a VertexBasedCellPopulation");
}

// Explicit instantiation
template class VoronoiDataWriter<1,1>;
template class VoronoiDataWriter<1,2>;
template class VoronoiDataWriter<2,2>;
template class VoronoiDataWriter<1,3>;
template class VoronoiDataWriter<2,3>;
template class VoronoiDataWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VoronoiDataWriter)
