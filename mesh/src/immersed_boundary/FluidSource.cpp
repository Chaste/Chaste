/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "FluidSource.hpp"

template<unsigned SPACE_DIM>
FluidSource<SPACE_DIM>::FluidSource(unsigned index, ChastePoint<SPACE_DIM> point)
    : mIndex(index),
      mLocation(point.rGetLocation()),
      mStrength(0.0),
      mIsSourceAssociatedWithElement(false),
      mAssociatedElementIndex(UINT_MAX)
{
}

template<unsigned SPACE_DIM>
FluidSource<SPACE_DIM>::FluidSource(
    unsigned index,
    c_vector<double, SPACE_DIM> location)
    : mIndex(index),
      mLocation(location),
      mStrength(0.0),
      mIsSourceAssociatedWithElement(false),
      mAssociatedElementIndex(UINT_MAX)
{
}

template<unsigned SPACE_DIM>
FluidSource<SPACE_DIM>::FluidSource(unsigned index, double v1, double v2, double v3)
    : mIndex(index),
      mStrength(0.0),
      mIsSourceAssociatedWithElement(false),
      mAssociatedElementIndex(UINT_MAX)
{
    mLocation[0] = v1;
    if (SPACE_DIM > 1)
    {
        mLocation[1] = v2;
        if (SPACE_DIM > 2)
        {
            mLocation[2] = v3;
        }
    }
}

template<unsigned SPACE_DIM>
FluidSource<SPACE_DIM>::~FluidSource()
{
}

template<unsigned SPACE_DIM>
unsigned FluidSource<SPACE_DIM>::GetIndex() const
{
    return mIndex;
}

template<unsigned SPACE_DIM>
void FluidSource<SPACE_DIM>::SetIndex(unsigned index)
{
    mIndex = index;
}

template<unsigned SPACE_DIM>
ChastePoint<SPACE_DIM> FluidSource<SPACE_DIM>::GetPoint() const
{
    return ChastePoint<SPACE_DIM>(mLocation);
}

template<unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& FluidSource<SPACE_DIM>::rGetLocation() const
{
    return mLocation;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM>& FluidSource<SPACE_DIM>::rGetModifiableLocation()
{
    return mLocation;
}

template<unsigned SPACE_DIM>
double FluidSource<SPACE_DIM>::GetStrength() const
{
    return mStrength;
}

template<unsigned SPACE_DIM>
void FluidSource<SPACE_DIM>::SetStrength(double strength)
{
    mStrength = strength;
}

template<unsigned SPACE_DIM>
void FluidSource<SPACE_DIM>::SetIfSourceIsAssociatedWithElement(bool associated)
{
    mIsSourceAssociatedWithElement = associated;
}

template<unsigned SPACE_DIM>
bool FluidSource<SPACE_DIM>::IsSourceAssociatedWithElement()
{
    return mIsSourceAssociatedWithElement;
}

template<unsigned SPACE_DIM>
unsigned FluidSource<SPACE_DIM>::GetAssociatedElementIndex() const
{
    return mIsSourceAssociatedWithElement ? mAssociatedElementIndex : UINT_MAX;
}

template<unsigned SPACE_DIM>
void FluidSource<SPACE_DIM>::SetAssociatedElementIndex(unsigned associatedElementIndex)
{
    mIsSourceAssociatedWithElement = true;
    mAssociatedElementIndex = associatedElementIndex;
}

// Explicit instantiation
template class FluidSource<1>;
template class FluidSource<2>;
template class FluidSource<3>;
