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
void SloughingCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
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
                c_vector<double, 2> location;
                location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
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

// Explicit instantiation
template class SloughingCellKiller<1>;
template class SloughingCellKiller<2>;
template class SloughingCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SloughingCellKiller)
