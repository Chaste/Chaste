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

#include "CellCycleModelProteinConcentrationsWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellCycleModelOdeHandler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellCycleModelProteinConcentrationsWriter<ELEMENT_DIM, SPACE_DIM>::CellCycleModelProteinConcentrationsWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("proteinconcentrations.dat")
{
    this->mVtkCellDataName = "Protein concentrations";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellCycleModelProteinConcentrationsWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    /*
     * At present it is not possible to output all cell variables via this method, we just return zero.
     * If the user requires cell variables to be output to VTK, the easiest way to do this is to store
     * them in the cell's CellData using a modifier object, since any CellData is automatically output
     * to VTK (if CHASTE_VTK is turned on).
     */
    return 0.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellCycleModelProteinConcentrationsWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    CellCycleModelOdeHandler* p_model = dynamic_cast<CellCycleModelOdeHandler*>(pCell->GetCellCycleModel());

    if (p_model)
    {
        // Write location index corresponding to cell
        *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

        // Write cell variables
        std::vector<double> proteins = p_model->GetProteinConcentrations();
        for (unsigned i=0; i<proteins.size(); i++)
        {
            *this->mpOutStream << proteins[i] << " ";
        }
    }

    else
    {
        EXCEPTION("CellCycleModelProteinConcentrationsWriter cannot be used with a cell-cycle model that does not inherit from CellCycleModelOdeHandler");
    }
}

// Explicit instantiation
template class CellCycleModelProteinConcentrationsWriter<1,1>;
template class CellCycleModelProteinConcentrationsWriter<1,2>;
template class CellCycleModelProteinConcentrationsWriter<2,2>;
template class CellCycleModelProteinConcentrationsWriter<1,3>;
template class CellCycleModelProteinConcentrationsWriter<2,3>;
template class CellCycleModelProteinConcentrationsWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellCycleModelProteinConcentrationsWriter)
