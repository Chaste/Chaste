/*

Copyright (C) University of Oxford, 2005-2010

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


#include "CellwiseNutrientSinkPde.hpp"

template<unsigned DIM>
CellwiseNutrientSinkPde<DIM>::CellwiseNutrientSinkPde(MeshBasedTissue<DIM>& rTissue, double coefficient)
    : mrTissue(rTissue),
      mCoefficient(coefficient)
{
}

template<unsigned DIM>
double CellwiseNutrientSinkPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX)
{
    return 0.0;
}

template<unsigned DIM>
double CellwiseNutrientSinkPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}

template<unsigned DIM>
double CellwiseNutrientSinkPde<DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode)
{
    TissueCell& r_cell = mrTissue.rGetCellUsingLocationIndex(rNode.GetIndex());
    if (r_cell.GetCellProliferativeType() != APOPTOTIC)
    {
        return -mCoefficient;
    }
    else
    {
        return 0.0;
    }
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> CellwiseNutrientSinkPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX)
{
    return identity_matrix<double>(DIM);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class CellwiseNutrientSinkPde<1>;
template class CellwiseNutrientSinkPde<2>;
template class CellwiseNutrientSinkPde<3>;
