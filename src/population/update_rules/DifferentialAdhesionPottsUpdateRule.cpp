/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "DifferentialAdhesionPottsUpdateRule.hpp"

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
    if( pCellA->HasCellProperty<CellLabel>() && pCellB->HasCellProperty<CellLabel>() )
    {
        return GetLabelledCellLabelledCellAdhesionEnergyParameter();
    }
    else if( pCellA->HasCellProperty<CellLabel>() || pCellB->HasCellProperty<CellLabel>() )
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
    if( pCell->HasCellProperty<CellLabel>() )
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
	*rParamsFile << "\t\t\t<LabelledCellLabelledCellAdhesionEnergyParameter>" << mLabelledCellLabelledCellAdhesionEnergyParameter << "</LabelledCellLabelledCellAdhesionEnergyParameter> \n";
    *rParamsFile << "\t\t\t<LabelledCellCellAdhesionEnergyParameter>" << mLabelledCellCellAdhesionEnergyParameter << "</LabelledCellCellAdhesionEnergyParameter> \n";
    *rParamsFile << "\t\t\t<LabelledCellBoundaryAdhesionEnergyParameter>" << mLabelledCellBoundaryAdhesionEnergyParameter << "</LabelledCellBoundaryAdhesionEnergyParameter> \n";

    // Call method on direct parent class
    AdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DifferentialAdhesionPottsUpdateRule<1>;
template class DifferentialAdhesionPottsUpdateRule<2>;
template class DifferentialAdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionPottsUpdateRule)
