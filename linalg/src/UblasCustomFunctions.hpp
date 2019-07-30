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

#ifndef UBLASCUSTOMFUNCTIONS_HPP_
#define UBLASCUSTOMFUNCTIONS_HPP_

/**
 * @file
 * A collection of useful functions extending the functionality of the
 * Boost Ublas library.
 */

#include "UblasIncludes.hpp"

#include "Exception.hpp"
#include "MathsCustomFunctions.hpp"

// COMMON DETERMINANTS - SQUARE

/**
 * 1x1 Determinant.
 * @return the determinant of a ublas matrix.
 *
 * @param rM The matrix of which to find the determinant.
 */
template<class T>
inline T Determinant(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM)
{
    using namespace boost::numeric::ublas;

    return rM(0,0);
}

/**
 * 2x2 Determinant.
 * @return the determinant of a ublas matrix.
 *
 * @param rM The matrix of which to find the determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T,2,2>& rM)
{
    using namespace boost::numeric::ublas;

    return rM(0,0)*rM(1,1) - rM(1,0)*rM(0,1);
}

/**
 * 3x3 Determinant.
 * @return the determinant of a ublas matrix.
 *
 * @param rM The matrix of which to find the determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM)
{
    using namespace boost::numeric::ublas;

    return    rM(0,0) * (rM(1,1)*rM(2,2) - rM(1,2)*rM(2,1))
            - rM(0,1) * (rM(1,0)*rM(2,2) - rM(1,2)*rM(2,0))
            + rM(0,2) * (rM(1,0)*rM(2,1) - rM(1,1)*rM(2,0));
}

// COMMON GENERALIZED DETERMINANTS - NOT SQUARE

/**
 * 3x2 (Generalized determinant).
 * @return calculated generalized determinant of a 3x2 matrix.
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM)
{
    using namespace boost::numeric::ublas;
    c_matrix<T,2,2> product = prod(trans(rM), rM);
    return std::sqrt(Determinant(product));
}

/**
 * 3x1 (Generalized determinant).
 * @return calculated generalized determinant of a 3x1 matrix.
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM)
{
    using namespace boost::numeric::ublas;
    return std::sqrt(rM(0,0)*rM(0,0) + rM(1,0)*rM(1,0) + rM(2,0)*rM(2,0));
}

/**
 * 2x1 (Generalized determinant).
 * @return calculated generalized determinant of a 2x1 matrix.
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM)
{
    using namespace boost::numeric::ublas;
    return   std::sqrt(rM(0,0) * rM(0,0) + rM(1,0) * rM(1,0));
}

/**
 * @return (nothing)
 * 3x0 (Generalized determinant) - not implement, but needed by some compilers.
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 0>& rM)
{
    NEVER_REACHED;
}

/**
 * @return (nothing)
 * 2x0 (Generalized determinant) - not implement, but needed by some compilers.
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 2, 0>& rM)
{
    NEVER_REACHED;
}

/**
 * @return (nothing)
 * 1x0 (Generalized determinant) - not implement, but needed by some compilers.
 *
 * @param rM The matrix of which to find the generalized determinant.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 1, 0>& rM)
{
    NEVER_REACHED;
}

// COMMON SUBDETERMINANTS - SQUARE

/**
 * 1x1 SubDeterminant.
 * @return the determinant of a submatrix after removing a particular row and column
 * For a 1x1 matrix this should always remove the only row and column (0,0).
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow == 0);
    assert(misscol == 0);
    return 1.0;
}

/**
 * 2x2 SubDeterminant.
 * @return the determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 2>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 2);
    assert(misscol < 2);

    unsigned row = (missrow==1) ? 0 : 1;
    unsigned col = (misscol==1) ? 0 : 1;
    return rM(row,col);
}

/**
 * SubDeterminant 3x3.
 * @return determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    assert(misscol < 3);

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

// COMMON SUBDETERMINANTS - NOT SQUARE

/**
 * SubDeterminant 3x2.
 * @return determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    //assert(misscol < 2); //Don't assert this since it is used

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

/**
 * SubDeterminant 3x1.
 * @return determinant of a submatrix after removing a particular row and column.
 * @param rM The matrix of which to find the subdeterminant.
 *
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    assert(misscol < 1);

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

/**
 * SubDeterminant 2x1.
 * @return determinant of a submatrix after removing a particular row and column.
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 2);
    assert(misscol < 1);

    unsigned row = (missrow==1) ? 0 : 1;
    unsigned col = (misscol==1) ? 0 : 1;
    return rM(row,col);
}

#if defined(__xlC__)
/* IBM compiler doesn't support zero-sized arrays*/
#else //#if defined(__xlC__)
/**
 * SubDeterminant 3x0 - Not implemented, but needed by some compilers for recursive template calls.
 * @return (nothing)
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}

/**
 * SubDeterminant 2x0 - Not implemented, but needed by some compilers for recursive template calls.
 * @return (nothing)
 *
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}

/**
 * SubDeterminant 1x0 - Not implemented, but needed by some compilers for recursive template calls.
 * Determinant of a submatrix after removing a particular row and column.
 * @return (nothing)
 * @param rM The matrix of which to find the subdeterminant.
 * @param missrow The index to the row to remove
 * @param misscol The index to the column to remove
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 1, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}
#endif //#if defined(__xlC__)

// COMMON INVERSES - SQUARE

/**
 * 1x1 Inverse.
 * Get the inverse of a ublas matrix.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 1> Inverse(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T,1,1> inverse;
    T det = Determinant(rM);
    assert(fabs(det) > DBL_EPSILON); // else it is a singular matrix
    inverse(0,0) =  1.0/det;
    return inverse;
}

/**
 * 2x2 Inverse.
 * Get the inverse of a ublas matrix.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 2, 2> Inverse(const boost::numeric::ublas::c_matrix<T, 2, 2>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 2, 2> inverse;
    T det = Determinant(rM);

    assert( fabs(det) > DBL_EPSILON ); // else it is a singular matrix
    inverse(0,0)  =  rM(1,1)/det;
    inverse(0,1)  = -rM(0,1)/det;
    inverse(1,0)  = -rM(1,0)/det;
    inverse(1,1)  =  rM(0,0)/det;
    return inverse;
}

/**
 * 3x3 Inverse.
 * Get the inverse of a ublas matrix.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 3, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 3, 3> inverse;
    T det = Determinant(rM);
    assert(fabs(det) > DBL_EPSILON); // else it is a singular matrix

    inverse(0,0) =  (rM(1,1)*rM(2,2) - rM(1,2)*rM(2,1))/det;
    inverse(1,0) = -(rM(1,0)*rM(2,2) - rM(1,2)*rM(2,0))/det;
    inverse(2,0) =  (rM(1,0)*rM(2,1) - rM(1,1)*rM(2,0))/det;
    inverse(0,1) = -(rM(0,1)*rM(2,2) - rM(0,2)*rM(2,1))/det;
    inverse(1,1) =  (rM(0,0)*rM(2,2) - rM(0,2)*rM(2,0))/det;
    inverse(2,1) = -(rM(0,0)*rM(2,1) - rM(0,1)*rM(2,0))/det;
    inverse(0,2) =  (rM(0,1)*rM(1,2) - rM(0,2)*rM(1,1))/det;
    inverse(1,2) = -(rM(0,0)*rM(1,2) - rM(0,2)*rM(1,0))/det;
    inverse(2,2) =  (rM(0,0)*rM(1,1) - rM(0,1)*rM(1,0))/det;

    return inverse;
}

// COMMON PSEUDO-INVERSES - NOT SQUARE

/**
 * 2x3 pseudo-inverse of a  matrix.
 * The pseudo-inverse is given by pinv(T) = (T'T)^(-1)*T'.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 2, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 2, 3> inverse;

    //
    // calculate (T'T)^-1, where T'T = (a b)
    //                                 (c d)

    T a = rM(0,0)*rM(0,0) + rM(1,0)*rM(1,0) + rM(2,0)*rM(2,0);
    T b = rM(0,0)*rM(0,1) + rM(1,0)*rM(1,1) + rM(2,0)*rM(2,1);
    T c = b;
    T d = rM(0,1)*rM(0,1) + rM(1,1)*rM(1,1) + rM(2,1)*rM(2,1);

    T det = a*d - b*c;

    T a_inv =  d/det;
    T b_inv = -b/det;
    T c_inv = -c/det;
    T d_inv =  a/det;

    inverse(0,0) = a_inv*rM(0,0) + b_inv*rM(0,1);
    inverse(1,0) = c_inv*rM(0,0) + d_inv*rM(0,1);
    inverse(0,1) = a_inv*rM(1,0) + b_inv*rM(1,1);
    inverse(1,1) = c_inv*rM(1,0) + d_inv*rM(1,1);
    inverse(0,2) = a_inv*rM(2,0) + b_inv*rM(2,1);
    inverse(1,2) = c_inv*rM(2,0) + d_inv*rM(2,1);

    return inverse;
}

/**
 * 2x1 pseudo-inverse of a  matrix.
 * The pseudo-inverse is given by pinv(T) = (T'T)^(-1)*T'.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 2> Inverse(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 1, 2> inverse;
    T det = Determinant(rM);

    inverse(0,0) = rM(0,0)/det/det;
    inverse(0,1) = rM(1,0)/det/det;

    return inverse;
}

/**
 * 3x1 pseudo-inverse of a  matrix.
 * The pseudo-inverse is given by pinv(T) = (T'T)^(-1)*T'.
 *
 * @param rM The matrix of which to find the inverse.
 * @return The inverse
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 1, 3> inverse;
    T det = Determinant(rM);

    inverse(0,0) = rM(0,0)/det/det;
    inverse(0,1) = rM(1,0)/det/det;
    inverse(0,2) = rM(2,0)/det/det;

    return inverse;
}

// COMMON MATRIX TRACES

/**
 * 1x1 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
inline T Trace(const c_matrix<T, 1, 1>& rM)
{
    return rM(0,0);
}

/**
 * 2x2 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
inline T Trace(const c_matrix<T, 2, 2>& rM)
{
    return rM(0,0) + rM(1,1);
}

/**
 * 3x3 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
inline T Trace(const c_matrix<T, 3, 3>& rM)
{
    return rM(0,0) + rM(1,1) + rM(2,2);
}

/**
 * 4x4 matrix trace (sum of diagonal elements).
 *
 * @param rM The matrix of which to find the trace.
 * @return The trace
 */
template<class T>
inline T Trace(const c_matrix<T, 4, 4>& rM)
{
    return rM(0,0) + rM(1,1) + rM(2,2) + rM(3,3);
}

// OTHER MATRIX FUNCTIONS (INVARIANTS, EIGENVECTORS)

/**
 * @return 3x3 second invariant.
 * @note Implementation only correct for
 * a SYMMETRIC matrix though. It is up to the user to check the
 * input matrix is symmetric.
 *
 * @param rM The matrix
 */
template<class T>
T SecondInvariant(const c_matrix<T, 3, 3>& rM)
{
    return    rM(0,0)*rM(1,1) + rM(1,1)*rM(2,2) + rM(2,2)*rM(0,0)
            - rM(1,0)*rM(1,0) - rM(2,1)*rM(2,1) - rM(2,0)*rM(2,0);
}

/**
 * @return 2x2 second invariant.
 * Second invariant of a 2d matrix, i.e. the determinant. This function
 * is mainly here just so that the same code can be used in 2d and 3d.
 *
 * @param rM The matrix
 */
template<class T>
T SecondInvariant(const c_matrix<T, 2, 2>& rM)
{
    return Determinant(rM);
}

/**
 * Find the eigenvector corresponding
 * real eigenvalue which is smallest in magnitude.
 * Caveat: if there are zero eigenvalues they are ignored.
 * It's the smallest magnitude non-zero real eigenvalue which is used.
 *
 * @param rA 3x3 matrix is question.  This should be symmetric and positive definite.
 * @return 3-vector corresponding to right-eigenvector in question
 */
c_vector<double,3> CalculateEigenvectorForSmallestNonzeroEigenvalue(c_matrix<double, 3, 3>& rA);

/**
 * Helper function to get maximum eigenpair from a 3x3 matrix by the power method
 * @param rA 3x3 matrix is question.
 * @param rEigenvector a guess eigenvector which will be refined
 * @return the maximum eigenvalue
 */
double CalculateMaxEigenpair(c_matrix<double, 3, 3>& rA, c_vector<double, 3>& rEigenvector);

                             //COMMON VECTOR FUNCTIONS


/**
 * This is a fake cross-product aka vector-product.  Fake because it's only implemented for 3-vectors.
 * This version is to satisfy template compilation.
 *
 * @param rA first vector
 * @param rB second vector
 * @return Does not return rA x rB
 */
template<class T>
c_vector<T, 1> VectorProduct(const c_vector<T, 1>& rA, const c_vector<T, 1>& rB)
{
    NEVER_REACHED;
}
/**
 * This is a fake cross-product aka vector-product.  Fake because it's only implemented for 3-vectors.
 * This version is to satisfy template compilation.
 *
 * @param rA first vector
 * @param rB second vector
 * @return Does not return rA x rB
 */
template<class T>
c_vector<T, 2> VectorProduct(const c_vector<T, 2>& rA, const c_vector<T, 2>& rB)
{
    NEVER_REACHED;
}

/**
 * This is a cross-product aka vector-product, only implemented for 3-vectors.
 *
 * @param rA first vector
 * @param rB second vector
 * @return rA x rB
 */
template<class T>
c_vector<T, 3> VectorProduct(const c_vector<T, 3>& rA, const c_vector<T, 3>& rB)
{

    c_vector<T, 3> result;

    result(0) = rA(1)*rB(2) - rA(2)*rB(1);
    result(1) = rA(2)*rB(0) - rA(0)*rB(2);
    result(2) = rA(0)*rB(1) - rA(1)*rB(0);

    return result;
}

/**
 * Convenience function for quickly creating test vectors (1D).
 *
 * @param x entry in vector
 * @return vector=(x)
 */
c_vector<double, 1> Create_c_vector(double x);

/**
 * Convenience function for quickly creating test vectors (2D).
 *
 * @param x entry in vector
 * @param y entry in vector
 * @return vector=(x,y)
 */
c_vector<double, 2> Create_c_vector(double x, double y);

/**
 * Convenience function for quickly creating test vectors (3D).
 *
 * @param x entry in vector
 * @param y entry in vector
 * @param z entry in vector
 * @return vector=(x,y,z)
 */
c_vector<double, 3> Create_c_vector(double x, double y, double z);

#endif /*UBLASCUSTOMFUNCTIONS_HPP_*/
