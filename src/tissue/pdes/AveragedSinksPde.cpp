/*

Copyright (C) University of Oxford, 2005-2009

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
#include "AveragedSinksPde.hpp"

template<unsigned DIM>
AveragedSinksPde<DIM>::AveragedSinksPde(MeshBasedTissue<DIM>& rTissue, double coefficient)
    : mrTissue(rTissue),
      mCoefficient(coefficient)
{
}

template<unsigned DIM>
void AveragedSinksPde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh) // must be called before solve
{
    // Allocate memory
    mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        mCellDensityOnCoarseElements[elem_index] = 0.0;
    }

    // Loop over cells, find which coarse element it is in, and add 1 to the mSourceTermOnCoarseElements[elem_index];
    for (typename AbstractTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
        cell_iter != mrTissue.End();
        ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = mrTissue.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);

        if (cell_iter->GetCellProliferativeType() != APOPTOTIC)
        {
            mCellDensityOnCoarseElements[elem_index] += 1.0;
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}

template<unsigned DIM>
double AveragedSinksPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX)
{
    return 0.0;
}

template<unsigned DIM>
double AveragedSinksPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement) // now takes in element
{
    assert(mCellDensityOnCoarseElements.size() > 0);
    return mCoefficient*mCellDensityOnCoarseElements[pElement->GetIndex()];
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> AveragedSinksPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX)
{
    return identity_matrix<double>(DIM);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AveragedSinksPde<1>;
template class AveragedSinksPde<2>;
template class AveragedSinksPde<3>;
