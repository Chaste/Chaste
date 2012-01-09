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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "QuadraticBasisFunction.hpp"
#include "Exception.hpp"

/**
 * Specialization for 0d.
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
double QuadraticBasisFunction<0>::ComputeBasisFunction(const ChastePoint<0>& rPoint, unsigned basisIndex)
{
    assert(basisIndex == 0);
    return 1.0;
}

/**
 * Specialization for 0d.
 *
 * @param rPoint The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @param rReturnValue Filled in with the values of the basis functions.
 */
void QuadraticBasisFunction<0>::ComputeBasisFunctions(const ChastePoint<0>& rPoint,
                                                      c_vector<double, 1>& rReturnValue)
{
    rReturnValue(0) = ComputeBasisFunction(rPoint, 0);
}

/**
 * Compute a basis function at a point within an element.
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
template <unsigned ELEMENT_DIM>
double QuadraticBasisFunction<ELEMENT_DIM>::ComputeBasisFunction(const ChastePoint<ELEMENT_DIM>& rPoint, unsigned basisIndex)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM >= 0);
    double x, y, z;
    switch (ELEMENT_DIM)
    {
    case 0:
        assert(basisIndex == 0);
        return 1.0;
        break;

    case 1:
        x = rPoint[0];
        switch (basisIndex)
        {
            case 0:
                return 2.0*(x-1.0)*(x-0.5);
                break;
            case 1:
                return 2.0*x*(x-0.5);
                break;
            case 2:
                return 4.0*x*(1.0-x);
                break;
            default:
                NEVER_REACHED;
        }
        break;

    case 2:
        x = rPoint[0];
        y = rPoint[1];
        switch (basisIndex)
        {
            case 0: // the node at (0,0)
                return 2.0 * (1.0 - x - y) * (0.5 - x - y);
                break;
            case 1: // the node at (1,0)
                return 2.0*x*(x-0.5);
                break;
            case 2: // the node at (0,1)
                return 2.0*y*(y-0.5);
                break;
            case 3: // the node opposite 0, which is (1/2,1/2)
                return 4.0 * y * x;
                break;
            case 4: // the node opposite 1, which is (0,1/2)
                return 4.0 * (1.0 - x - y) * y;
                break;
            case 5: // the node opposite 2, which is (1/2,0)
                return 4.0 * (1.0 - x - y) * x;
                break;
            default:
                NEVER_REACHED;
        }
        break;

    case 3:
        x = rPoint[0];
        y = rPoint[1];
        z = rPoint[2];
        switch (basisIndex)
        {
            case 0: // the node at (0,0,0)
                return 2.0 * (1.0 - x - y - z) * (0.5 - x - y - z);
                break;
            case 1: // the node at (1,0,0)
                return 2.0*x*(x-0.5);
                break;
            case 2: // the node at (0,1,0)
                return 2.0*y*(y-0.5);
                break;
            case 3: // the node at (0,0,1)
                return 2.0*z*(z-0.5);
                break;
            case 4: // our (tetgen convention), node4 is between nodes 0 and 1, (1/2,0,0)
                return 4.0 * (1.0 - x - y - z) * x;
                break;
            case 5: // our (tetgen convention), node5 is between nodes 1 and 2, (1/2,1/2,0)
                return 4 * x * y;
                break;
            case 6: // our (tetgen convention), node6 is between nodes 0 and 2, (0,1/2,0)
                return 4.0 * (1.0 - x - y - z) * y;
                break;
            case 7: // our (tetgen convention), node7 is between nodes 0 and 3, (0,0,1/2)
                return 4.0 * (1.0 - x - y - z) * z;
                break;
            case 8: // our (tetgen convention), node8 is between nodes 1 and 3, (1/2,0,1/2)
                return 4.0 * x * z;
                break;
            case 9: // our (tetgen convention), node9 is between nodes 2 and 3, (0,1/2,1/2)
                return 4.0 * y * z;
                break;
            default:
                NEVER_REACHED;
        }
        break;
    }
    return 0.0; // Avoid compiler warning
}

/**
 * Compute the derivative of a basis function at a point within an canonical element.
 *
 * @param rPoint The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector (c_vector
 *     instance) giving the derivative along each axis.
 */
template <unsigned ELEMENT_DIM>
c_vector<double, ELEMENT_DIM> QuadraticBasisFunction<ELEMENT_DIM>::ComputeBasisFunctionDerivative(const ChastePoint<ELEMENT_DIM>& rPoint, unsigned basisIndex)
{
    c_vector<double, ELEMENT_DIM> gradN;
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM > 0);

    double x, y, z;
    switch (ELEMENT_DIM)
    {
    case 1:
        x = rPoint[0];
        switch (basisIndex)
        {
            case 0:
                gradN(0) = 4.0*x-3.0;
                break;
            case 1:
                gradN(0) = 4.0*x-1.0;
                break;
            case 2:
                gradN(0) = 4.0-8.0*x;
                break;
            default:
                NEVER_REACHED;
        }
        break;

    case 2:
        x = rPoint[0];
        y = rPoint[1];
        switch (basisIndex)
        {
            case 0:
                gradN(0) = -3.0 + 4.0*x + 4.0*y;
                gradN(1) = -3.0 + 4.0*x + 4.0*y;
                break;
            case 1:
                gradN(0) = 4.0*x - 1.0;
                gradN(1) = 0.0;
                break;
            case 2:
                gradN(0) = 0.0;
                gradN(1) = 4.0*y - 1.0;
                break;
            case 3:
                gradN(0) = 4.0*y;
                gradN(1) = 4.0*x;
                break;
            case 4:
                gradN(0) = -4.0*y;
                gradN(1) = 4.0-4.0*x-8.0*y;
                break;
            case 5:
                gradN(0) = 4.0-8.0*x-4.0*y;
                gradN(1) = -4.0*x;
                break;
            default:
                NEVER_REACHED;
        }
        break;

    case 3:
        x = rPoint[0];
        y = rPoint[1];
        z = rPoint[2];
        switch (basisIndex)
        {
            case 0:
                gradN(0) = -3.0 + 4.0*(x+y+z);
                gradN(1) = -3.0 + 4.0*(x+y+z);
                gradN(2) = -3.0 + 4.0*(x+y+z);
                break;
            case 1:
                gradN(0) =  4.0*x-1.0;
                gradN(1) =  0;
                gradN(2) =  0;
                break;
            case 2:
                gradN(0) =  0;
                gradN(1) =  4.0*y-1.0;
                gradN(2) =  0;
                break;
            case 3:
                gradN(0) =  0;
                gradN(1) =  0;
                gradN(2) =  4.0*z-1.0;
                break;
            case 4:
                gradN(0) =  4.0-8.0*x-4.0*y-4.0*z;
                gradN(1) =  -4.0*x;
                gradN(2) =  -4.0*x;
                break;
            case 5:
                gradN(0) =  4.0*y;
                gradN(1) =  4.0*x;
                gradN(2) =  0.0;
                break;
            case 6:
                gradN(0) =  -4.0*y;
                gradN(1) =  4.0-4.0*x-8.0*y-4.0*z;
                gradN(2) =  -4.0*y;
                break;
            case 7:
                gradN(0) =  -4.0*z;
                gradN(1) =  -4.0*z;
                gradN(2) =  4.0-4.0*x-4.0*y-8.0*z;
                break;
            case 8:
                gradN(0) =  4.0*z;
                gradN(1) =  0;
                gradN(2) =  4.0*x;
                break;
            case 9:
                gradN(0) =  0;
                gradN(1) =  4.0*z;
                gradN(2) =  4.0*y;
                break;
            default:
                NEVER_REACHED;
        }
        break;
    }
    return gradN;
}

/**
 * Compute all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @param rReturnValue The values of the basis functions, in local index order.
 */
template <unsigned ELEMENT_DIM>
void QuadraticBasisFunction<ELEMENT_DIM>::ComputeBasisFunctions(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                             c_vector<double, (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2>& rReturnValue)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM >= 0);

    for (unsigned i=0; i<(ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2; i++)
    {
        rReturnValue(i) = ComputeBasisFunction(rPoint, i);
    }
}

/**
 * Compute the derivatives of all basis functions at a point within an element.
 *
 * @param rPoint The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @param rReturnValue The derivatives of the basis functions, in local index order. Each
 *     column of the matrix gives the derivative along
 *     each axis.
 */
template <unsigned ELEMENT_DIM>
void QuadraticBasisFunction<ELEMENT_DIM>::ComputeBasisFunctionDerivatives(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                                    c_matrix<double, ELEMENT_DIM, (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2>& rReturnValue)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM > 0);

    for (unsigned j=0; j<(ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2; j++)
    {
        matrix_column<c_matrix<double, ELEMENT_DIM, (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2> > column(rReturnValue, j);
        column = ComputeBasisFunctionDerivative(rPoint, j);
    }
}

/**
 * Compute the derivatives of all basis functions at a point within an element.
 * This method will transform the results, for use within gaussian quadrature
 * for example.
 *
 * @param rPoint The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @param rInverseJacobian The inverse of the Jacobian matrix mapping the real
 *     element into the canonical element.
 * @param rReturnValue The derivatives of the basis functions, in local index order. Each
 *     entry is a vector (VectorDouble instance) giving the derivative along
 *     each axis.
 */
template <unsigned ELEMENT_DIM>
void QuadraticBasisFunction<ELEMENT_DIM>::ComputeTransformedBasisFunctionDerivatives(const ChastePoint<ELEMENT_DIM>& rPoint,
                                                                                  const c_matrix<double, ELEMENT_DIM, ELEMENT_DIM>& rInverseJacobian,
                                                                                  c_matrix<double, ELEMENT_DIM, (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2>& rReturnValue)
{
    assert(ELEMENT_DIM < 4 && ELEMENT_DIM > 0);

    ComputeBasisFunctionDerivatives(rPoint, rReturnValue);
    rReturnValue = prod(trans(rInverseJacobian), rReturnValue);
}

//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class QuadraticBasisFunction<1>;
template class QuadraticBasisFunction<2>;
template class QuadraticBasisFunction<3>;
