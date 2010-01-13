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
#include "OxygenBasedCellKiller.hpp"

template <unsigned SPACE_DIM>
OxygenBasedCellKiller<SPACE_DIM>::OxygenBasedCellKiller(AbstractTissue<SPACE_DIM>* pTissue,
                                                        double concentration)
    : AbstractCellKiller<SPACE_DIM>(pTissue),
      mHypoxicConcentration(concentration)
{
}

template <unsigned SPACE_DIM>
void OxygenBasedCellKiller<SPACE_DIM>::SetHypoxicConcentration(double hypoxicConcentration)
{
    mHypoxicConcentration = hypoxicConcentration;
}

template <unsigned SPACE_DIM>
double OxygenBasedCellKiller<SPACE_DIM>::GetHypoxicConcentration() const
{
    return mHypoxicConcentration;
}

template <unsigned SPACE_DIM>
void OxygenBasedCellKiller<SPACE_DIM>::TestAndLabelSingleCellForApoptosis(TissueCell& rCell)
{
    if (rCell.GetCellProliferativeType()==APOPTOTIC && !(rCell.HasApoptosisBegun()))
    {
        rCell.StartApoptosis();
    }
}

template <unsigned SPACE_DIM>
void OxygenBasedCellKiller<SPACE_DIM>::TestAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractTissue<SPACE_DIM>::Iterator cell_iter = this->mpTissue->Begin();
        cell_iter != this->mpTissue->End();
        ++cell_iter)
    {
        TestAndLabelSingleCellForApoptosis(*cell_iter);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OxygenBasedCellKiller<1>;
template class OxygenBasedCellKiller<2>;
template class OxygenBasedCellKiller<3>;
