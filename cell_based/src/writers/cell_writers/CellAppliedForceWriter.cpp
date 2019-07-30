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

#include "CellAppliedForceWriter.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "UblasVectorInclude.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellAppliedForceWriter<ELEMENT_DIM, SPACE_DIM>::CellAppliedForceWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellappliedforce.dat")
{
    this->mVtkVectorCellDataName = "Cell applied force";
    this->mOutputScalarData = false;
    this->mOutputVectorData = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> CellAppliedForceWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    c_vector<double, SPACE_DIM> applied_force = scalar_vector<double>(SPACE_DIM, DOUBLE_UNSET);

    if (dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(pCellPopulation))
    {
        unsigned node_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
        applied_force = pCellPopulation->GetNode(node_index)->rGetAppliedForce();
    }

    return applied_force;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellAppliedForceWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    c_vector<double, SPACE_DIM> applied_force = GetVectorCellDataForVtkOutput(pCell, pCellPopulation);

    *this->mpOutStream << location_index << " " << cell_id << " ";
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << applied_force[i] << " ";
    }
}

// Explicit instantiation
template class CellAppliedForceWriter<1,1>;
template class CellAppliedForceWriter<1,2>;
template class CellAppliedForceWriter<2,2>;
template class CellAppliedForceWriter<1,3>;
template class CellAppliedForceWriter<2,3>;
template class CellAppliedForceWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellAppliedForceWriter)
