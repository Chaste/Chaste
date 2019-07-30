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
#ifndef _GAUSSIANQUADRATURERULE_HPP_
#define _GAUSSIANQUADRATURERULE_HPP_

#include <vector>
#include "ChastePoint.hpp"

/**
 * This class encapsulates tables of Gaussian quadrature points and the
 * associated weights.
 *
 * Data is available for 1d, 2d and 3d quadrature over (canonical) triangles,
 * with appropriate numbers of Gauss points.  Weights sum to 1 and are non-negative.
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
     * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
     */
    GaussianQuadratureRule(unsigned quadratureOrder);

    /**
     * Get a quadrature point.
     *
     * @param index The index of the point to return.
     * @return A gaussian quadrature point.
     */
    const ChastePoint<ELEMENT_DIM>& rGetQuadPoint(unsigned index) const;

    /**
     * @return the weight associated with a quadrature point.
     *
     * @param index The index of the point to return.
     */
    double GetWeight(unsigned index) const;

    /**
     * @return the number of quadrature points. This is the number of points in
     * each dimension, raised to the power of the number of dimensions.
     */
    unsigned GetNumQuadPoints() const;
};

#endif //_GAUSSIANQUADRATURERULE_HPP_
