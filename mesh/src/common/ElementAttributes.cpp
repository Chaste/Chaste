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

#include "ElementAttributes.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementAttributes<ELEMENT_DIM, SPACE_DIM>::ElementAttributes()
    :   mAttributes(std::vector<double>())
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double>& ElementAttributes<ELEMENT_DIM, SPACE_DIM>::rGetAttributes()
{
    return mAttributes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementAttributes<ELEMENT_DIM, SPACE_DIM>::AddAttribute(double attribute)
{
    mAttributes.push_back(attribute);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementAttributes<ELEMENT_DIM, SPACE_DIM>::SetFirstAttribute(double attribute)
{
    // Make sure that the first entry exists
    if (mAttributes.empty())
    {
        mAttributes.resize(1u);
    }

    mAttributes[0] = attribute;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ElementAttributes<ELEMENT_DIM, SPACE_DIM>::GetFirstAttribute()
{
    if (mAttributes.empty())
    {
        ///\todo #2739 This should throw: EXCEPTION("Attempting to get element attribute when there are none defined");
        return 0.0;
    }
    // Otherwise
    return(mAttributes[0]);
}

// Explicit instantiation
template class ElementAttributes<0,1>;
template class ElementAttributes<1,1>;
template class ElementAttributes<0,2>;
template class ElementAttributes<1,2>;
template class ElementAttributes<2,2>;
template class ElementAttributes<0,3>;
template class ElementAttributes<1,3>;
template class ElementAttributes<2,3>;
template class ElementAttributes<3,3>;
