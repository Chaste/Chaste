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

#include "UblasCustomFunctions.hpp"
#include "LinearBasisFunction.hpp"
#include "ChastePoint.hpp"
#include <cassert>

/**
 * Compute a basis function at a point within an element (3d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
template <>
double LinearBasisFunction<3>::ComputeBasisFunction(
    const ChastePoint<3>& rPoint,
    unsigned basisIndex)
{
    assert(basisIndex <= 3);

    switch (basisIndex)
    {
        case 0:
            return 1.0 - rPoint[0] - rPoint[1] - rPoint[2];
            break;
        case 1:
            return rPoint[0];
            break;
        case 2:
            return rPoint[1];
            break;
        case 3:
            return rPoint[2];
            break;
        default:
            NEVER_REACHED; //not possible to get here because of assertions above
    }

    return 0.0; // LCOV_EXCL_LINE // Avoid compiler warning
}

/**
 * Compute a basis function at a point within an element (2d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
template <>
double LinearBasisFunction<2>::ComputeBasisFunction(
    const ChastePoint<2>& rPoint,
    unsigned basisIndex)
{
    assert(basisIndex <= 2);

    switch (basisIndex)
    {
        case 0:
            return 1.0 - rPoint[0] - rPoint[1];
            break;
        case 1:
            return rPoint[0];
            break;
        case 2:
            return rPoint[1];
            break;
        default:
           NEVER_REACHED; //not possible to get here because of assertions above
    }
    return 0.0; // LCOV_EXCL_LINE // Avoid compiler warning
}

/**
 * Compute a basis function at a point within an element (1d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
template <>
double LinearBasisFunction<1>::ComputeBasisFunction(
    const ChastePoint<1>& rPoint,
    unsigned basisIndex)
{
    assert(basisIndex <= 1);

    switch (basisIndex)
    {
        case 0:
            return 1.0 - rPoint[0];
            break;
        case 1:
            return rPoint[0];
            break;
        default:
           NEVER_REACHED; //not possible to get here because of assertions above
    }
    return 0.0; // LCOV_EXCL_LINE  // Avoid compiler warning
}

/**
 * Compute a basis function at a point within an element (0d case).
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
double LinearBasisFunction<0>::ComputeBasisFunction(const ChastePoint<0>& rPoint, unsigned basisIndex)
{
    assert(basisIndex == 0);
    return 1.0;
}

/**
 * Compute the derivative of a basis function at a point within a
 * canonical element (3d case).
 *
 * @param rPoint (unused) The point at which to compute the basis function.
 *     The results are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector
 *     (c_vector<double, ELEMENT_DIM> instance) giving the derivative
 *     along each axis.
 */
template <>
c_vector<double, 3> LinearBasisFunction<3>::ComputeBasisFunctionDerivative(
    const ChastePoint<3>& rPoint,
    unsigned basisIndex)
{
    assert(basisIndex <= 3);

    c_vector<double, 3> gradN;
    switch (basisIndex)
    {
        case 0:
            gradN(0) = -1;
            gradN(1) = -1;
            gradN(2) = -1;
            break;
        case 1:
            gradN(0) =  1;
            gradN(1) =  0;
            gradN(2) =  0;
            break;
        case 2:
            gradN(0) =  0;
            gradN(1) =  1;
            gradN(2) =  0;
            break;
        case 3:
            gradN(0) =  0;
            gradN(1) =  0;
            gradN(2) =  1;
            break;
        default:
           ; //not possible to get here because of assertions above
    }
    return gradN;
}

/**
 * Compute the derivative of a basis function at a point within a
 * canonical element (2d case).
 *
 * @param rPoint (unused) The point at which to compute the basis function.
 *     The results are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector
 *     (c_vector<double, ELEMENT_DIM> instance) giving the derivative
 *     along each axis.
 */
template <>
c_vector<double, 2> LinearBasisFunction<2>::ComputeBasisFunctionDerivative(
    const ChastePoint<2>& rPoint,
    unsigned basisIndex)
{
    assert(basisIndex <= 2);

    c_vector<double, 2> gradN;
    switch (basisIndex)
    {
        case 0:
            gradN(0) = -1;
            gradN(1) = -1;
            break;
        case 1:
            gradN(0) =  1;
            gradN(1) =  0;
            break;
        case 2:
            gradN(0) =  0;
            gradN(1) =  1;
            break;
        default:
           ; //not possible to get here because of assertions above
    }
    return gradN;
}

/**
 * Compute the derivative of a basis function at a point within a
 * canonical element (1d case).
 *
 * @param rPoint (unused) The point at which to compute the basis function.
 *     The results are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector
 *     (c_vector<double, ELEMENT_DIM> instance) giving the derivative
 *     along each axis.
 */
template <>
c_vector<double,1> LinearBasisFunction<1>::ComputeBasisFunctionDerivative(
    const ChastePoint<1>& rPoint,
    unsigned basisIndex)
{
    assert(basisIndex <= 1);

    c_vector<double,1> gradN;
    switch (basisIndex)
    {
        case 0:
            gradN(0) = -1;
            break;
        case 1:
            gradN(0) =  1;
            break;
        default:
           ; //not possible to get here because of assertions above
    }
    return gradN;
}


/**
 * Compute all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @param rReturnValue A reference to a vector, to be filled in
 */
template <unsigned ELEMENT_DIM>
void LinearBasisFunction<ELEMENT_DIM>::ComputeBasisFunctions(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                          c_vector<double, ELEMENT_DIM+1>& rReturnValue)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM > 0);
    for (unsigned i=0; i<ELEMENT_DIM+1; i++)
    {
        rReturnValue(i) = ComputeBasisFunction(rPoint, i);
    }
}

/**
 * Compute all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @param rReturnValue A reference to a vector, to be filled in
 *
 */
void LinearBasisFunction<0>::ComputeBasisFunctions(const ChastePoint<0>& rPoint,
                                                   c_vector<double,1>& rReturnValue)
{
    rReturnValue(0) = ComputeBasisFunction(rPoint, 0);
}

/**
 * Compute the derivatives of all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @param rReturnValue A reference to a vector, to be filled in
 * @return The derivatives of the basis functions as the column vectors of
 *     a matrix in local index order.
 */
template <unsigned ELEMENT_DIM>
void LinearBasisFunction<ELEMENT_DIM>::ComputeBasisFunctionDerivatives(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                                    c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rReturnValue)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM > 0);

    for (unsigned j=0; j<ELEMENT_DIM+1; j++)
    {
        matrix_column<c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> > column(rReturnValue, j);
        column = ComputeBasisFunctionDerivative(rPoint, j);
    }
}

/**
 * Compute the derivatives of all basis functions at a point within an element.
 * This method will transform the results, for use within gaussian quadrature
 * for example.
 *
 * @param rPoint The point at which to compute the basis functions. The
 *     results are undefined if this is not within the canonical element.
 * @param rInverseJacobian The inverse of the Jacobian matrix mapping the real
 *     element into the canonical element.
 * @param rReturnValue A reference to a vector, to be filled in
 * @return The derivatives of the basis functions, in local index order. Each
 *     entry is a vector (c_vector<double, SPACE_DIM> instance) giving the
 *     derivative along each axis.
 */
template <unsigned ELEMENT_DIM>
void LinearBasisFunction<ELEMENT_DIM>::ComputeTransformedBasisFunctionDerivatives(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                                               const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM>& rInverseJacobian,
                                                                               c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>& rReturnValue)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM > 0);

    ComputeBasisFunctionDerivatives(rPoint, rReturnValue);
    rReturnValue = prod(trans(rInverseJacobian), rReturnValue);
}

// Explicit instantiation
template class LinearBasisFunction<1>;
template class LinearBasisFunction<2>;
template class LinearBasisFunction<3>;
