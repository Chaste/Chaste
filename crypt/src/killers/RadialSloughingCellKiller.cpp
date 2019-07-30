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
#include "RadialSloughingCellKiller.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

RadialSloughingCellKiller::RadialSloughingCellKiller(AbstractCellPopulation<2>* pCellPopulation, c_vector<double,2> centre, double radius)
    : AbstractCellKiller<2>(pCellPopulation),
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

void RadialSloughingCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Get distance from centre of cell population
        double r = norm_2(this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter) - mCentre);

        if (r > mRadius)
        {
            cell_iter->Kill();
        }
    }
}

void RadialSloughingCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<xCentre>" << mCentre[0] << "</xCentre>\n";
    *rParamsFile << "\t\t\t<yCentre>" << mCentre[1] << "</yCentre>\n";
    *rParamsFile << "\t\t\t<mRadius>" << mRadius << "</mRadius>\n";

    // Call method on direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialSloughingCellKiller)
