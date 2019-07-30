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

#include "TargetedCellKiller.hpp"

template<unsigned DIM>
TargetedCellKiller<DIM>::TargetedCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, unsigned targetedIndex, bool bloodLust)
: AbstractCellKiller<DIM>(pCellPopulation),
  mTargetIndex(targetedIndex),
  mBloodLust(bloodLust)
{
}

template<unsigned DIM>
unsigned TargetedCellKiller<DIM>::GetTargetIndex() const
{
    return mTargetIndex;
}

template<unsigned DIM>
unsigned TargetedCellKiller<DIM>::GetBloodLust() const
{
    return mBloodLust;
}

template<unsigned DIM>
void TargetedCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    if (!mBloodLust || this->mpCellPopulation->GetNumRealCells()==0 || this->mpCellPopulation->GetNumRealCells()<mTargetIndex+1)
    {
        return;
    }
    this->mpCellPopulation->GetCellUsingLocationIndex(mTargetIndex)->Kill();
    mBloodLust = false;
}

template<unsigned DIM>
void TargetedCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<TargetIndex>" << mTargetIndex << "</TargetIndex>\n";
    *rParamsFile << "\t\t\t<BloodLust>" << mBloodLust << "</BloodLust>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class TargetedCellKiller<1>;
template class TargetedCellKiller<2>;
template class TargetedCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TargetedCellKiller)
