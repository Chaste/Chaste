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
#ifndef _GAUSSIANQUADRATURERULE_HPP_
#define _GAUSSIANQUADRATURERULE_HPP_

#include <vector>
#include "ChastePoint.hpp"

/**
 * This class encapsulates tables of gaussian quadrature points and the
 * associated weights.
 *
 * Data is available for 1d, 2d and 3d quadrature over (canonical) triangles,
 * with between 1 and 3 (inclusive) Gauss points in each dimension.
 * The values are computed when an object is instantiated.
 */
template<unsigned ELEMENT_DIM>
class GaussianQuadratureRule
{
    /** The total number of Gauss points. */
    unsigned mNumQuadPoints;

    /** The gaussian quadrature points. */
    std::vector<ChastePoint<ELEMENT_DIM> > mPoints;

    /** The associated weights. */
    std::vector<double> mWeights;

public:

    /**
     * The constructor builds the appropriate table for the dimension (given
     * by the template argument) and number of points in each dimension (given
     * as a constructor argument).
     *
     * An exception is thrown if data is not available for the requested
     * parameters.
     *
     * @param numPointsInEachDimension  number of Gauss points in each dimension
     */
    GaussianQuadratureRule(unsigned numPointsInEachDimension);

    /**
     * Get a quadrature point.
     *
     * @param index The index of the point to return.
     * @return A gaussian quadrature point.
     */
    const ChastePoint<ELEMENT_DIM>& rGetQuadPoint(unsigned index) const;

    /**
     * Get the weight associated with a quadrature point.
     *
     * @param index The index of the point to return.
     */
    double GetWeight(unsigned index) const;

    /**
     * Get the number of quadrature points. This is the number of points in
     * each dimension, raised to the power of the number of dimensions.
     */
    unsigned GetNumQuadPoints() const;
};

#endif //_GAUSSIANQUADRATURERULE_HPP_
