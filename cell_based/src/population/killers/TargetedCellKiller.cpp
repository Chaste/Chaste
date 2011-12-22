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
void TargetedCellKiller<DIM>::TestAndLabelCellsForApoptosisOrDeath()
{
    if ( !mBloodLust || this->mpCellPopulation->GetNumRealCells()==0 || this->mpCellPopulation->GetNumRealCells()<mTargetIndex+1)
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

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class TargetedCellKiller<1>;
template class TargetedCellKiller<2>;
template class TargetedCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TargetedCellKiller)
