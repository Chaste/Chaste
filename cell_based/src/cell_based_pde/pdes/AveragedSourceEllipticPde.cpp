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

#include "AveragedSourceEllipticPde.hpp"
#include "ApoptoticCellProperty.hpp"

template<unsigned DIM>
AveragedSourceEllipticPde<DIM>::AveragedSourceEllipticPde(AbstractCellPopulation<DIM>& rCellPopulation,
                                                          double sourceCoefficient,
                                                          double diffusionCoefficient)
    : mrCellPopulation(rCellPopulation),
      mSourceCoefficient(sourceCoefficient),
      mDiffusionCoefficient(diffusionCoefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM>& AveragedSourceEllipticPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double AveragedSourceEllipticPde<DIM>::GetCoefficient() const
{
    return mSourceCoefficient;
}

template<unsigned DIM>
void AveragedSourceEllipticPde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map< CellPtr, unsigned >* pCellPdeElementMap) // must be called before solve
{
    // Allocate memory
    mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        mCellDensityOnCoarseElements[elem_index] = 0.0;
    }

    // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        unsigned elem_index = 0;
        const ChastePoint<DIM>& r_position_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

        if (pCellPdeElementMap != nullptr)
        {
            elem_index = (*pCellPdeElementMap)[*cell_iter];
        }
        else
        {
            elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
        }

        // Update element map if cell has moved
        bool cell_is_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();

        if (!cell_is_apoptotic)
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
double AveragedSourceEllipticPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return 0.0;
}

template<unsigned DIM>
double AveragedSourceEllipticPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement) // now takes in element
{
    assert(!mCellDensityOnCoarseElements.empty());
    return mSourceCoefficient * mCellDensityOnCoarseElements[pElement->GetIndex()];
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> AveragedSourceEllipticPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX)
{
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

template<unsigned DIM>
double AveragedSourceEllipticPde<DIM>::GetUptakeRateForElement(unsigned elementIndex)
{
    return this->mCellDensityOnCoarseElements[elementIndex];
}

// Explicit instantiation
template class AveragedSourceEllipticPde<1>;
template class AveragedSourceEllipticPde<2>;
template class AveragedSourceEllipticPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AveragedSourceEllipticPde)
