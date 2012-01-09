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
#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "ChastePoint.hpp"

/**
 * Linear basis functions for the finite element method,
 * computed on a canonical element.
 */
template <unsigned ELEMENT_DIM>
class LinearBasisFunction
{
public:

    /**
     * Compute a basis function at a point within an element.
     *
     * @param rPoint The point at which to compute the basis function. The results
     *     are undefined if this is not within the canonical element.
     * @param basisIndex Which basis function to compute. This is a local index
     *     within a canonical element.
     * @return The value of the basis function.
     */
    static double ComputeBasisFunction(const ChastePoint<ELEMENT_DIM>& rPoint, unsigned basisIndex);

    /**
     * Compute the derivative of a basis function at a point within a
     * canonical element.
     *
     * @param rPoint (unused) The point at which to compute the basis function.
     *     The results are undefined if this is not within the canonical element.
     * @param basisIndex Which basis function to compute. This is a local index
     *     within a canonical element.
     * @return The derivative of the basis function. This is a vector
     *     (c_vector<double, ELEMENT_DIM> instance) giving the derivative
     *     along each axis.
     */
    static c_vector<double, ELEMENT_DIM> ComputeBasisFunctionDerivative(const ChastePoint<ELEMENT_DIM>& rPoint, unsigned basisIndex);

    static void ComputeBasisFunctions(const ChastePoint<ELEMENT_DIM>& rPoint, c_vector<double, ELEMENT_DIM+1>& rReturnValue);
    static void ComputeBasisFunctionDerivatives(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rReturnValue);

    static void ComputeTransformedBasisFunctionDerivatives(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                           const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM>& rInverseJacobian,
                                                           c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rReturnValue);
};

/**
 * We need to specialise for the 0d case, because 0x0 matrices don't work.
 */
template <>
class LinearBasisFunction<0>
{
public:
    static double ComputeBasisFunction(const ChastePoint<0>& rPoint, unsigned basisIndex);
    static void ComputeBasisFunctions(const ChastePoint<0>& rPoint,c_vector<double,1>& rReturnValue);
};

#endif //_LINEARBASISFUNCTION_HPP_
