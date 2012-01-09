/*

Copyright (C) University of Oxford, 2005-2012

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

//////////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////////

template class ChastePoint<1>;
template class ChastePoint<2>;
template class ChastePoint<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChastePoint)

