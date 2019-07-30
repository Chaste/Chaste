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

#include "ChastePoint.hpp"
#include "Exception.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
ChastePoint<DIM>::ChastePoint(double v1, double v2, double v3)
{
    if (DIM > 0)
    {
        mLocation[0] = v1;
    }
    if (DIM > 1)
    {
        mLocation[1] = v2;
    }
    if (DIM > 2)
    {
        mLocation[2] = v3;
    }
}

template<unsigned DIM>
ChastePoint<DIM>::ChastePoint(std::vector<double> coords)
{
    for (unsigned i=0; i<DIM; i++)
    {
        mLocation(i) = coords.at(i);
    }
}

template<unsigned DIM>
ChastePoint<DIM>::ChastePoint(c_vector<double, DIM> location)
    : mLocation(location)
{
}

template<unsigned DIM>
c_vector<double, DIM>& ChastePoint<DIM>::rGetLocation()
{
    return mLocation;
}

template<unsigned DIM>
const c_vector<double, DIM>& ChastePoint<DIM>::rGetLocation() const
{
    return mLocation;
}

template<unsigned DIM>
double ChastePoint<DIM>::operator[] (unsigned i) const
{
    assert(i<DIM);
    return mLocation(i);
}

template<unsigned DIM>
double ChastePoint<DIM>::GetWithDefault(unsigned i, double def) const
{
    if (i<DIM)
    {
        return mLocation(i);
    }
    else
    {
        return def;
    }
}

template<unsigned DIM>
void ChastePoint<DIM>::SetCoordinate(unsigned i, double value)
{
    assert(i < DIM);
    mLocation(i) = value;
}

template<unsigned DIM>
bool ChastePoint<DIM>::IsSamePoint(const ChastePoint<DIM>& rPoint) const
{
    bool returned_value = true;
    for (unsigned dim=0; dim<DIM; dim++)
    {
        if (rPoint[dim] != mLocation[dim])
        {
            returned_value = false;
            break;
        }
    }
    return returned_value;
}

// Methods of the 0d version

ChastePoint<0>::ChastePoint(double v1, double v2, double v3)
{
}

double ChastePoint<0>::operator[] (unsigned i) const
{
    EXCEPTION("Zero-dimensional point has no data");
}

// Explicit instantiation
template class ChastePoint<1>;
template class ChastePoint<2>;
template class ChastePoint<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChastePoint)
