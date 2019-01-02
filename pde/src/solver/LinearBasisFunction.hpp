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
#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "ChastePoint.hpp"
#include "UblasMatrixInclude.hpp"

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
