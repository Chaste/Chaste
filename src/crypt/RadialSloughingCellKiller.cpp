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
#include "RadialSloughingCellKiller.hpp"
#include "AbstractCellCentreBasedTissue.hpp"

RadialSloughingCellKiller::RadialSloughingCellKiller(AbstractTissue<2>* pTissue, c_vector<double,2> centre, double radius)
    : AbstractCellKiller<2>(pTissue),
      mCentre(centre),
      mRadius(radius)
{
}

c_vector<double,2> RadialSloughingCellKiller::GetCentre() const
{
    return mCentre;
}

double RadialSloughingCellKiller::GetRadius() const
{
    return mRadius;
}

void RadialSloughingCellKiller::TestAndLabelCellsForApoptosisOrDeath()
{
    for (AbstractTissue<2>::Iterator cell_iter = this->mpTissue->Begin();
         cell_iter != this->mpTissue->End();
         ++cell_iter)
    {
        // Get distance from centre of tissue
        double r = norm_2(this->mpTissue->GetLocationOfCellCentre(*cell_iter) - mCentre);

        if (r > mRadius)
        {
            cell_iter->Kill();
        }
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialSloughingCellKiller)
