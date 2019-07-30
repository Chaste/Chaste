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

#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "CellLabel.hpp"

template<unsigned DIM>
DifferentialAdhesionPottsUpdateRule<DIM>::DifferentialAdhesionPottsUpdateRule()
    : AdhesionPottsUpdateRule<DIM>(),
      mLabelledCellLabelledCellAdhesionEnergyParameter(0.1), // Educated guess
      mLabelledCellCellAdhesionEnergyParameter(0.1), // Educated guess
      mLabelledCellBoundaryAdhesionEnergyParameter(0.2) // Educated guess
{
}

template<unsigned DIM>
DifferentialAdhesionPottsUpdateRule<DIM>::~DifferentialAdhesionPottsUpdateRule()
{
}

template<unsigned DIM>
double DifferentialAdhesionPottsUpdateRule<DIM>::GetCellCellAdhesionEnergy(CellPtr pCellA, CellPtr pCellB)
{
    if (pCellA->HasCellProperty<CellLabel>() && pCellB->HasCellProperty<CellLabel>())
    {
        return GetLabelledCellLabelledCellAdhesionEnergyParameter();
    }
    else if (pCellA->HasCellProperty<CellLabel>() || pCellB->HasCellProperty<CellLabel>())
    {
        return GetLabelledCellCellAdhesionEnergyParameter();
    }
    else
    {
        return this->GetCellCellAdhesionEnergyParameter();
    }
}

template<unsigned DIM>
double DifferentialAdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergy(CellPtr pCell)
{
    if (pCell->HasCellProperty<CellLabel>())
    {
        return GetLabelledCellBoundaryAdhesionEnergyParameter();
    }
    else
    {
        return this->GetCellBoundaryAdhesionEnergyParameter();
    }
}

template<unsigned DIM>
double DifferentialAdhesionPottsUpdateRule<DIM>::GetLabelledCellLabelledCellAdhesionEnergyParameter()
{
    return mLabelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionPottsUpdateRule<DIM>::GetLabelledCellCellAdhesionEnergyParameter()
{
    return mLabelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionPottsUpdateRule<DIM>::GetLabelledCellBoundaryAdhesionEnergyParameter()
{
    return mLabelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionPottsUpdateRule<DIM>::SetLabelledCellLabelledCellAdhesionEnergyParameter(double labelledCellLabelledCellAdhesionEnergyParameter)
{
    mLabelledCellLabelledCellAdhesionEnergyParameter = labelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionPottsUpdateRule<DIM>::SetLabelledCellCellAdhesionEnergyParameter(double labelledCellCellAdhesionEnergyParameter)
{
    mLabelledCellCellAdhesionEnergyParameter = labelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionPottsUpdateRule<DIM>::SetLabelledCellBoundaryAdhesionEnergyParameter(double labelledCellBoundaryAdhesionEnergyParameter)
{
    mLabelledCellBoundaryAdhesionEnergyParameter = labelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<LabelledCellLabelledCellAdhesionEnergyParameter>" << mLabelledCellLabelledCellAdhesionEnergyParameter << "</LabelledCellLabelledCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellCellAdhesionEnergyParameter>" << mLabelledCellCellAdhesionEnergyParameter << "</LabelledCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellBoundaryAdhesionEnergyParameter>" << mLabelledCellBoundaryAdhesionEnergyParameter << "</LabelledCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialAdhesionPottsUpdateRule<1>;
template class DifferentialAdhesionPottsUpdateRule<2>;
template class DifferentialAdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionPottsUpdateRule)
