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

#include "ChasteCuboid.hpp"
#include "Exception.hpp"

template <unsigned SPACE_DIM>
    ChasteCuboid<SPACE_DIM>::ChasteCuboid(ChastePoint<SPACE_DIM>& rLowerPoint, ChastePoint<SPACE_DIM>& rUpperPoint)
        : mLowerCorner(rLowerPoint),
          mUpperCorner(rUpperPoint)
    {
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            if (mLowerCorner[dim] > mUpperCorner[dim])
            {
                EXCEPTION("Attempt to create a cuboid with MinCorner greater than MaxCorner in some dimension");
            }
        }
    }

template <unsigned SPACE_DIM>
    bool ChasteCuboid<SPACE_DIM>::DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const
    {
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            if (rPointToCheck[dim] < mLowerCorner[dim] - 100*DBL_EPSILON
                || mUpperCorner[dim] + 100* DBL_EPSILON < rPointToCheck[dim])
            {
                return false;
            }
        }
        return true;
    }

template <unsigned SPACE_DIM>
    const ChastePoint<SPACE_DIM>& ChasteCuboid<SPACE_DIM>::rGetUpperCorner() const
    {
        return mUpperCorner;
    }

template <unsigned SPACE_DIM>
    const ChastePoint<SPACE_DIM>& ChasteCuboid<SPACE_DIM>::rGetLowerCorner() const
    {
        return mLowerCorner;
    }

template <unsigned SPACE_DIM>
    double ChasteCuboid<SPACE_DIM>::GetWidth(unsigned rDimension) const
    {
        assert(rDimension<SPACE_DIM);
        return mUpperCorner[rDimension] - mLowerCorner[rDimension];
    }

template <unsigned SPACE_DIM>
     unsigned ChasteCuboid<SPACE_DIM>::GetLongestAxis() const
     {
        unsigned axis=0;
        double max_dimension = 0.0;
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            double dimension =  mUpperCorner[i] - mLowerCorner[i];
            if ( dimension > max_dimension)
            {
                axis=i;
                max_dimension = dimension;
            }
        }
        return axis;
     }

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class ChasteCuboid<1>;
template class ChasteCuboid<2>;
template class ChasteCuboid<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChasteCuboid)

