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

#include "CellwiseSourceEllipticPde.hpp"

template<unsigned DIM>
CellwiseSourceEllipticPde<DIM>::CellwiseSourceEllipticPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation, double sourceCoefficient)
    : mrCellPopulation(rCellPopulation),
      mSourceCoefficient(sourceCoefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& CellwiseSourceEllipticPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::GetCoefficient() const
{
    return mSourceCoefficient;
}

template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return 0.0;
}

// LCOV_EXCL_START
template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
double CellwiseSourceEllipticPde<DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode)
{
    double source_coefficient = 0.0;

    if (mrCellPopulation.IsPdeNodeAssociatedWithNonApoptoticCell(rNode.GetIndex()))
    {
        source_coefficient = mSourceCoefficient;
    }

    return source_coefficient;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> CellwiseSourceEllipticPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX)
{
    return identity_matrix<double>(DIM);
}

// Explicit instantiation
template class CellwiseSourceEllipticPde<1>;
template class CellwiseSourceEllipticPde<2>;
template class CellwiseSourceEllipticPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellwiseSourceEllipticPde)
