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

#include "ChasteEllipsoid.hpp"
#include "Exception.hpp"

template <unsigned SPACE_DIM>
ChasteEllipsoid<SPACE_DIM>::ChasteEllipsoid(ChastePoint<SPACE_DIM>& rCentre, ChastePoint<SPACE_DIM>& rRadii)
    : mCentre(rCentre),
      mRadii(rRadii)
{
    for (unsigned dim=0; dim<SPACE_DIM; dim++)
    {
        if (mRadii[dim] < 0.0)
        {
            EXCEPTION("Attempted to create an ellipsoid with a negative radius");
        }
    }
}

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
bool ChasteEllipsoid<1>::DoesContain(const ChastePoint<1>& rPointToCheck) const
{
    if (rPointToCheck[0] < mCentre[0] - mRadii[0] - 100.0 * DBL_EPSILON ||
        rPointToCheck[0] > mCentre[0] + mRadii[0] + 100.0 * DBL_EPSILON )
    {
        return false;
    }
    else
    {
        return true;
    }
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
bool ChasteEllipsoid<2>::DoesContain(const ChastePoint<2>& rPointToCheck) const
{
    double x_distance = (rPointToCheck[0]-mCentre[0])/mRadii[0];
    double y_distance = (rPointToCheck[1]-mCentre[1])/mRadii[1];

    if ( (x_distance*x_distance + y_distance*y_distance) > (1.0 + 100.0 * DBL_EPSILON) )
    {
        return false;
    }
    else
    {
        return true;
    }
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
bool ChasteEllipsoid<3>::DoesContain(const ChastePoint<3>& rPointToCheck) const
{
    double x_distance = (rPointToCheck[0]-mCentre[0])/mRadii[0];
    double y_distance = (rPointToCheck[1]-mCentre[1])/mRadii[1];
    double z_distance = (rPointToCheck[2]-mCentre[2])/mRadii[2];

    if ( (x_distance*x_distance + y_distance*y_distance + z_distance*z_distance) > (1.0 + 100.0 * DBL_EPSILON) )
    {
        return false;
    }
    else
    {
        return true;
    }
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

template <unsigned SPACE_DIM>
const ChastePoint<SPACE_DIM>& ChasteEllipsoid<SPACE_DIM>::rGetCentre() const
{
    return mCentre;
}

template <unsigned SPACE_DIM>
const ChastePoint<SPACE_DIM>& ChasteEllipsoid<SPACE_DIM>::rGetRadii() const
{
    return mRadii;
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ChasteEllipsoid<1>;
template class ChasteEllipsoid<2>;
template class ChasteEllipsoid<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChasteEllipsoid)
