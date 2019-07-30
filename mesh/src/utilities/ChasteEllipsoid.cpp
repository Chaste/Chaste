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

    if ((x_distance*x_distance + y_distance*y_distance) > (1.0 + 100.0 * DBL_EPSILON))
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

    if ((x_distance*x_distance + y_distance*y_distance + z_distance*z_distance) > (1.0 + 100.0 * DBL_EPSILON))
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


///////// Explicit instantiation///////

template class ChasteEllipsoid<1>;
template class ChasteEllipsoid<2>;
template class ChasteEllipsoid<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChasteEllipsoid)
