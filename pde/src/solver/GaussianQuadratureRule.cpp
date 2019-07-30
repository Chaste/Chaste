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

#include <cmath>

#include "GaussianQuadratureRule.hpp"
#include "Exception.hpp"
#include "UblasCustomFunctions.hpp"

template<unsigned ELEMENT_DIM>
const ChastePoint<ELEMENT_DIM>& GaussianQuadratureRule<ELEMENT_DIM>::rGetQuadPoint(unsigned index) const
{
    assert(index < mNumQuadPoints);
    return mPoints[index];
}

template<unsigned ELEMENT_DIM>
double GaussianQuadratureRule<ELEMENT_DIM>::GetWeight(unsigned index) const
{
    assert(index < mNumQuadPoints);
    return mWeights[index];
}

template<unsigned ELEMENT_DIM>
unsigned GaussianQuadratureRule<ELEMENT_DIM>::GetNumQuadPoints() const
{
    return mNumQuadPoints;
}

/**
 * Constructor specialization for 0d.
 *
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly (ignored in 0-d case)
 */
template<>
GaussianQuadratureRule<0>::GaussianQuadratureRule(unsigned quadratureOrder)
{
    mNumQuadPoints = 1;
    mWeights.push_back(1);
    mPoints.push_back(ChastePoint<0>());
}

/**
 * Constructor specialization for 1d.
 *
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 */
template<>
GaussianQuadratureRule<1>::GaussianQuadratureRule(unsigned quadratureOrder)
{
    switch (quadratureOrder)
    {
        case 0:
        case 1: // 1d, 1st order
                // 1 point rule
            mWeights.push_back(1);
            mPoints.push_back(ChastePoint<1>(0.5));
            break;
        case 2:
        case 3: // 1d, 3rd order
                // 2 point rule
            mWeights.push_back(0.5);
            mWeights.push_back(0.5);
            {
                double sqrt_one_third = sqrt(1.0/3.0);
                mPoints.push_back(ChastePoint<1>((-sqrt_one_third+1.0)/2.0));
                mPoints.push_back(ChastePoint<1>((sqrt_one_third+1.0)/2.0));
            }
            break;
        case 4:
        case 5: // 1d, 5th order
                // 3 point rule
            mWeights.push_back(5.0/18.0);
            mWeights.push_back(4.0/9.0);
            mWeights.push_back(5.0/18.0);

            {
                double sqrt_three_fifths = sqrt(3.0/5.0);
                mPoints.push_back(ChastePoint<1>((-sqrt_three_fifths+1.0)/2.0));
                mPoints.push_back(ChastePoint<1>(0.5));
                mPoints.push_back(ChastePoint<1>((sqrt_three_fifths+1.0)/2.0));
            }
            break;
         default:
            EXCEPTION("Gauss quadrature order not supported.");
    }
    assert(mPoints.size() == mWeights.size());
    mNumQuadPoints = mPoints.size();
}

/**
 * Constructor specialization for 2d.
 *
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 *
 */
template<>
GaussianQuadratureRule<2>::GaussianQuadratureRule(unsigned quadratureOrder)
{
    double one_third = 1.0/3.0;
    double one_sixth = 1.0/6.0;
    double two_thirds = 2.0/3.0;
    switch (quadratureOrder)
    {
        case 0: // 2d, 0th order
        case 1: // 2d, 1st order
                // 1 point rule
            mWeights.push_back(0.5);
            mPoints.push_back(ChastePoint<2>(one_third, one_third));
            break;

        case 2: // 2d, 2nd order
                // 3 point rule
            mWeights.push_back(one_sixth);
            mWeights.push_back(one_sixth);
            mWeights.push_back(one_sixth);

            mPoints.push_back(ChastePoint<2>(two_thirds, one_sixth));
            mPoints.push_back(ChastePoint<2>(one_sixth,  one_sixth));
            mPoints.push_back(ChastePoint<2>(one_sixth,  two_thirds));
            break;

        case 3: // 2d, 3rd order - derived by hand and using a Macsyma script to solve the cubic
                //                60*x^3  - 60*x^2  + 15*x - 1;
                // 6 point rule
            {
                double w = 1.0/12.0;
                mWeights.push_back(w);
                mWeights.push_back(w);
                mWeights.push_back(w);
                mWeights.push_back(w);
                mWeights.push_back(w);
                mWeights.push_back(w);

                double inverse_tan = atan(0.75);
                double cos_third = cos(inverse_tan/3.0);
                double sin_third = sin(inverse_tan/3.0);
                // a = 0.23193336461755
                double a =   sin_third/(2*sqrt(3.0)) - cos_third/6.0 + 1.0/3.0;
                // b = 0.10903901046622
                double b = - sin_third/(2*sqrt(3.0)) - cos_third/6.0 + 1.0/3.0;
                // c = 0.659028
                double c = cos_third/3.0 + 1.0/3.0;

                mPoints.push_back(ChastePoint<2>(a, b));
                mPoints.push_back(ChastePoint<2>(a, c));
                mPoints.push_back(ChastePoint<2>(b, a));
                mPoints.push_back(ChastePoint<2>(b, c));
                mPoints.push_back(ChastePoint<2>(c, a));
                mPoints.push_back(ChastePoint<2>(c, b));
            }
            break;
        default:
            EXCEPTION("Gauss quadrature order not supported.");
    }
    assert(mPoints.size() == mWeights.size());
    mNumQuadPoints = mPoints.size();
}

/**
 * Constructor specialization for 3d.
 *
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 */
template<>
GaussianQuadratureRule<3>::GaussianQuadratureRule(unsigned quadratureOrder)
{
    switch (quadratureOrder)
    {
        case 0: // 3d, 0th order
        case 1: // 3d, 1st order
                // 1 point rule
            mWeights.push_back(1.0/6.0);
            mPoints.push_back(ChastePoint<3>(0.25, 0.25, 0.25));
            break;

        case 2: //2nd order
                // 4 point rule
            {
                double sqrt_fifth = 1.0/sqrt(5.0);
                double a = (1.0 + 3.0*sqrt_fifth)/4.0; //0.585410196624969;
                double b = (1.0 - sqrt_fifth)/4.0; //0.138196601125011;
                double w = 1.0/24.0;

                mWeights.push_back(w);
                mWeights.push_back(w);
                mWeights.push_back(w);
                mWeights.push_back(w);

                mPoints.push_back(ChastePoint<3>(a,b,b));
                mPoints.push_back(ChastePoint<3>(b,a,b));
                mPoints.push_back(ChastePoint<3>(b,b,a));
                mPoints.push_back(ChastePoint<3>(b,b,b));
            }
            break;

        case 3: // 3d, 3rd order
            // 8 point rule
            /* The main options were
             *  5-point rule.  Commonly published rule has four symmetric points and
             *                 a negative weight in the centre.  We would like to avoid
             *                 negative weight (certainly for interpolation.
             *  8-point rule.  Uses two sets of symmetric points (as 4 point rule with a,b and then with c,d).
             *                 This one is hard to derive a closed form solution to.
             */
            {
                double root_seventeen = sqrt(17.0);
                double root_term = sqrt(1022.0-134.0*root_seventeen);
                double b = (55.0 - 3.0*root_seventeen + root_term)/196; //b = 0.328055
                double d = (55.0 - 3.0*root_seventeen - root_term)/196; //d = 0.106952

                double a = 1.0 - 3.0*b; // a = 0.0158359
                double c = 1.0 - 3.0*d; // c = 0.679143

                // w1 = 0.023088 (= 0.138528/6)
                double w1 = (20.0*d*d - 10.0*d + 1.0)/(240.0*(2.0*d*d - d - 2.0*b*b + b)); // w1 = 0.0362942
                double w2 = 1.0/24.0 - w1; // w2 = 0.0185787 (=0.111472/6)

                mWeights.push_back(w1);
                mWeights.push_back(w1);
                mWeights.push_back(w1);
                mWeights.push_back(w1);

                mWeights.push_back(w2);
                mWeights.push_back(w2);
                mWeights.push_back(w2);
                mWeights.push_back(w2);

                mPoints.push_back(ChastePoint<3>(a, b, b));
                mPoints.push_back(ChastePoint<3>(b, a, b));
                mPoints.push_back(ChastePoint<3>(b, b, a));
                mPoints.push_back(ChastePoint<3>(b, b, b));

                mPoints.push_back(ChastePoint<3>(c, d, d));
                mPoints.push_back(ChastePoint<3>(d, c, d));
                mPoints.push_back(ChastePoint<3>(d, d, c));
                mPoints.push_back(ChastePoint<3>(d, d, d));

            }
break;

        default:
            EXCEPTION("Gauss quadrature order not supported.");
    }
    assert(mPoints.size() == mWeights.size());
    mNumQuadPoints = mPoints.size();
}

template<unsigned ELEMENT_DIM>
GaussianQuadratureRule<ELEMENT_DIM>::GaussianQuadratureRule(unsigned quadratureOrder)
{
    EXCEPTION("Gauss quadrature rule not available for this dimension.");
}

// Explicit instantiation
template class GaussianQuadratureRule<0>;
template class GaussianQuadratureRule<1>;
template class GaussianQuadratureRule<2>;
template class GaussianQuadratureRule<3>;
template class GaussianQuadratureRule<4>;
