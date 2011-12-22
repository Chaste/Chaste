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
#include "SloughingCellKiller.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "Exception.hpp"

template<unsigned DIM>
SloughingCellKiller<DIM>::SloughingCellKiller(AbstractCellPopulation<DIM>* pCrypt, double sloughHeight, bool sloughSides, double sloughWidth)
    : AbstractCellKiller<DIM>(pCrypt),
      mSloughSides(sloughSides)
{
    assert(sloughHeight > 0.0);
    mSloughHeight = sloughHeight;

    assert(sloughWidth > 0.0);
    mSloughWidth = sloughWidth;
}

template<unsigned DIM>
bool SloughingCellKiller<DIM>::GetSloughSides() const
{
    return mSloughSides;
}

template<unsigned DIM>
double SloughingCellKiller<DIM>::GetSloughHeight() const
{
    return mSloughHeight;
}

template<unsigned DIM>
double SloughingCellKiller<DIM>::GetSloughWidth() const
{
    return mSloughWidth;
}


template<unsigned DIM>
void SloughingCellKiller<DIM>::TestAndLabelCellsForApoptosisOrDeath()
{
    switch (DIM)
    {
        case 1:
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                 cell_iter != this->mpCellPopulation->End();
                 ++cell_iter)
            {
                double x = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];

                if (x > mSloughHeight)
                {
                    cell_iter->Kill();
                }
            }
            break;
        }
        case 2:
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                 cell_iter != this->mpCellPopulation->End();
                 ++cell_iter)
            {
                c_vector<double, 2> location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
                double x = location[0];
                double y = location[1];

                if ((y>mSloughHeight) || (mSloughSides && ((x<0.0) || (x>mSloughWidth))))
                {
                    cell_iter->Kill();
                }
            }
            break;
        }
        case 3:
        {
            EXCEPTION("SloughingCellKiller is not yet implemented in 3D");
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;
    }
}

template<unsigned DIM>
void SloughingCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SloughLength>" << mSloughHeight << "</SloughLength>\n";
    *rParamsFile << "\t\t\t<SloughSides>" << mSloughSides << "</SloughSides>\n";
    *rParamsFile << "\t\t\t<SloughWidth>" << mSloughWidth << "</SloughWidth>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SloughingCellKiller<1>;
template class SloughingCellKiller<2>;
template class SloughingCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SloughingCellKiller)
