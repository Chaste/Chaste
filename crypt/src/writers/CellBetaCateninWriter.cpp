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

#include "CellBetaCateninWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellBetaCateninWriter<ELEMENT_DIM, SPACE_DIM>::CellBetaCateninWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizbetacatenin")
{
    this->mVtkCellDataName = "Beta catenin";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellBetaCateninWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(pCell->GetCellCycleModel());
    double b_cat_cytoplasm = p_model->GetCytoplasmicBetaCateninLevel();

    return b_cat_cytoplasm;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellBetaCateninWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE

    unsigned global_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    double x = pCellPopulation->GetLocationOfCellCentre(pCell)[0];
    double y = pCellPopulation->GetLocationOfCellCentre(pCell)[1];

    AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(pCell->GetCellCycleModel());
    double b_cat_membrane = p_model->GetMembraneBoundBetaCateninLevel();
    double b_cat_cytoplasm = p_model->GetCytoplasmicBetaCateninLevel();
    double b_cat_nuclear = p_model->GetNuclearBetaCateninLevel();

    *this->mpOutStream << global_index << " " << x << " " << y << " " << b_cat_membrane << " " << b_cat_cytoplasm << " " << b_cat_nuclear << " ";
}

// Explicit instantiation
template class CellBetaCateninWriter<1,1>;
template class CellBetaCateninWriter<1,2>;
template class CellBetaCateninWriter<2,2>;
template class CellBetaCateninWriter<1,3>;
template class CellBetaCateninWriter<2,3>;
template class CellBetaCateninWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellBetaCateninWriter)
