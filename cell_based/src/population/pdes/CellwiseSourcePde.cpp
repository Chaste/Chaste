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

#include "CellwiseSourcePde.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Exception.hpp"

template<unsigned DIM>
CellwiseSourcePde<DIM>::CellwiseSourcePde(MeshBasedCellPopulation<DIM>& rCellPopulation, double coefficient)
    : mrCellPopulation(rCellPopulation),
      mCoefficient(coefficient)
{
}

template<unsigned DIM>
const MeshBasedCellPopulation<DIM>& CellwiseSourcePde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double CellwiseSourcePde<DIM>::GetCoefficient() const
{
    return mCoefficient;
}

template<unsigned DIM>
double CellwiseSourcePde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return 0.0;
}

template<unsigned DIM>
double CellwiseSourcePde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}

template<unsigned DIM>
double CellwiseSourcePde<DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode)
{
    double coefficient = 0.0;

    CellPtr p_cell = mrCellPopulation.GetCellUsingLocationIndex(rNode.GetIndex());

    bool cell_is_apoptotic = p_cell->HasCellProperty<ApoptoticCellProperty>();

    if (!cell_is_apoptotic)
    {
        coefficient = mCoefficient;
    }

    return coefficient;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> CellwiseSourcePde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX)
{
    return identity_matrix<double>(DIM);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellwiseSourcePde<1>;
template class CellwiseSourcePde<2>;
template class CellwiseSourcePde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellwiseSourcePde)
