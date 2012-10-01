/*

Copyright (c) 2005-2012, University of Oxford.
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
 * @param deprecated  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly (ignored in 0-d case)
 */
template<>
GaussianQuadratureRule<0>::GaussianQuadratureRule(unsigned deprecated, unsigned quadratureOrder)
{
    mNumQuadPoints = 1;
    mWeights.push_back(1);
    mPoints.push_back(ChastePoint<0>());
}

/**
 * Constructor specialization for 1d.
 *
 * @param deprecated  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 */
template<>
GaussianQuadratureRule<1>::GaussianQuadratureRule(unsigned deprecated, unsigned quadratureOrder)
{
    switch (quadratureOrder)
    {
		case 0:
        case 1: // 1d, 1st order
            mWeights.push_back(1);
            mPoints.push_back(ChastePoint<1>(0.5));
            break;
        case 2:
        case 3: // 1d, 3rd order
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
 * @param deprecated  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 *
 */
template<>
GaussianQuadratureRule<2>::GaussianQuadratureRule(unsigned deprecated, unsigned quadratureOrder)
{
    double one_third = 1.0/3.0;
    double one_sixth = 1.0/6.0;
    double two_thirds = 2.0/3.0;
    switch (quadratureOrder)
    {
        case 0: // 2d, 0th order
        case 1: // 2d, 1st order
        	//This is now order 1
            mNumQuadPoints = 1;
            mWeights.push_back(0.5);
            mPoints.push_back(ChastePoint<2>(one_third, one_third));
            break;

        case 2: // 2d, 2nd order
        	mNumQuadPoints = 3;
            mWeights.push_back(one_sixth);
            mWeights.push_back(one_sixth);
            mWeights.push_back(one_sixth);
            
            mPoints.push_back(ChastePoint<2>(two_thirds, one_sixth));
            mPoints.push_back(ChastePoint<2>(one_sixth,  one_sixth));
            mPoints.push_back(ChastePoint<2>(one_sixth,  two_thirds));
            break;

        case 3: // 2d, 3rd order - derived by hand and using a Macsyma script to solve the cubic
                //                60*x^3  - 60*x^2  + 15*x - 1; .
            mNumQuadPoints = 6;
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
        //case 3:
        case 4: // 2d, 4th order
        	mNumQuadPoints = 9;
            mWeights.push_back(0.06846437766975);
            mWeights.push_back(0.10954300427160);
            mWeights.push_back(0.06846437766975);
            mWeights.push_back(0.06172839506173);
            mWeights.push_back(0.09876543209877);
            mWeights.push_back(0.06172839506173);
            mWeights.push_back(0.00869611615741);
            mWeights.push_back(0.01391378585185);
            mWeights.push_back(0.00869611615741);

            mPoints.push_back(ChastePoint<2>(0.10000000001607,0.11270166540000));
            mPoints.push_back(ChastePoint<2>(0.44364916730000,0.11270166540000));
            mPoints.push_back(ChastePoint<2>(0.78729833458393,0.11270166540000));
            mPoints.push_back(ChastePoint<2>(0.05635083270000,0.50000000000000));
            mPoints.push_back(ChastePoint<2>(0.25000000000000,0.50000000000000));
            mPoints.push_back(ChastePoint<2>(0.44364916730000,0.50000000000000));
            mPoints.push_back(ChastePoint<2>(0.01270166538393,0.88729833460000));
            mPoints.push_back(ChastePoint<2>(0.05635083270000,0.88729833460000));
            mPoints.push_back(ChastePoint<2>(0.10000000001607,0.88729833460000));
            break;
        default:
            EXCEPTION("Gauss quadrature order not supported.");
    }
    assert(mPoints.size() == mWeights.size());
    assert(mNumQuadPoints == mPoints.size());
}

/**
 * Constructor specialization for 3d.
 *
 * @param deprecated  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 */
template<>
GaussianQuadratureRule<3>::GaussianQuadratureRule(unsigned deprecated, unsigned quadratureOrder)
{
    switch (quadratureOrder)
    {
        case 0: // 3d, 0th order
            mWeights.push_back(1.0/6.0);
            mPoints.push_back(ChastePoint<3>(0.25000000000000,0.50000000000000,0.12500000000000));
            break;

        case 1: // 3d, 1st order
            mWeights.push_back(0.06132032652029);
            mWeights.push_back(0.01643073197073);
            mWeights.push_back(0.00440260136261);
            mWeights.push_back(0.00117967347971);
            mWeights.push_back(0.06132032652029);
            mWeights.push_back(0.01643073197073);
            mWeights.push_back(0.00440260136261);
            mWeights.push_back(0.00117967347971);

            mPoints.push_back(ChastePoint<3>(0.16666666666667,   0.21132486540519,   0.13144585576580));
            mPoints.push_back(ChastePoint<3>(0.62200846792815,   0.21132486540519,   0.03522081090086));
            mPoints.push_back(ChastePoint<3>(0.04465819873852,   0.78867513459481,   0.03522081090086));
            mPoints.push_back(ChastePoint<3>(0.16666666666667,   0.78867513459481,   0.00943738783766));
            mPoints.push_back(ChastePoint<3>(0.16666666666667,   0.21132486540519,   0.49056261216234));
            mPoints.push_back(ChastePoint<3>(0.62200846792815,   0.21132486540519,   0.13144585576580));
            mPoints.push_back(ChastePoint<3>(0.04465819873852,   0.78867513459481,   0.13144585576580));
            mPoints.push_back(ChastePoint<3>(0.16666666666667,   0.78867513459481,   0.03522081090086));
            break;

        case 2: //2nd order 4 point rule
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
            mWeights.push_back(0.01497274736603);
            mWeights.push_back(0.01349962850795);
            mWeights.push_back(0.00190178826891);
            mWeights.push_back(0.00760715307442);
            mWeights.push_back(0.00685871056241);
            mWeights.push_back(0.00096623512860);
            mWeights.push_back(0.00024155878219);
            mWeights.push_back(0.00021779261632);
            mWeights.push_back(0.00003068198821);
            mWeights.push_back(0.02395639578565);
            mWeights.push_back(0.02159940561273);
            mWeights.push_back(0.00304286123026);
            mWeights.push_back(0.01217144491907);
            mWeights.push_back(0.01097393689986);
            mWeights.push_back(0.00154597620576);
            mWeights.push_back(0.00038649405150);
            mWeights.push_back(0.00034846818612);
            mWeights.push_back(0.00004909118114);
            mWeights.push_back(0.01497274736603);
            mWeights.push_back(0.01349962850795);
            mWeights.push_back(0.00190178826891);
            mWeights.push_back(0.00760715307442);
            mWeights.push_back(0.00685871056241);
            mWeights.push_back(0.00096623512860);
            mWeights.push_back(0.00024155878219);
            mWeights.push_back(0.00021779261632);
            mWeights.push_back(0.00003068198821);

            mPoints.push_back(ChastePoint<3>(0.10000000001607,   0.11270166540000,   0.08872983347426));
            mPoints.push_back(ChastePoint<3>(0.44364916730000,   0.11270166540000,   0.05000000000803));
            mPoints.push_back(ChastePoint<3>(0.78729833458393,   0.11270166540000,   0.01127016654181));
            mPoints.push_back(ChastePoint<3>(0.05635083270000,   0.50000000000000,   0.05000000000803));
            mPoints.push_back(ChastePoint<3>(0.25000000000000,   0.50000000000000,   0.02817541635000));
            mPoints.push_back(ChastePoint<3>(0.44364916730000,   0.50000000000000,   0.00635083269197));
            mPoints.push_back(ChastePoint<3>(0.01270166538393,   0.88729833460000,   0.01127016654181));
            mPoints.push_back(ChastePoint<3>(0.05635083270000,   0.88729833460000,   0.00635083269197));
            mPoints.push_back(ChastePoint<3>(0.10000000001607,   0.88729833460000,   0.00143149884212));
            mPoints.push_back(ChastePoint<3>(0.10000000001607,   0.11270166540000,   0.39364916729197));
            mPoints.push_back(ChastePoint<3>(0.44364916730000,   0.11270166540000,   0.22182458365000));
            mPoints.push_back(ChastePoint<3>(0.78729833458393,   0.11270166540000,   0.05000000000803));
            mPoints.push_back(ChastePoint<3>(0.05635083270000,   0.50000000000000,   0.22182458365000));
            mPoints.push_back(ChastePoint<3>(0.25000000000000,   0.50000000000000,   0.12500000000000));
            mPoints.push_back(ChastePoint<3>(0.44364916730000,   0.50000000000000,   0.02817541635000));
            mPoints.push_back(ChastePoint<3>(0.01270166538393,   0.88729833460000,   0.05000000000803));
            mPoints.push_back(ChastePoint<3>(0.05635083270000,   0.88729833460000,   0.02817541635000));
            mPoints.push_back(ChastePoint<3>(0.10000000001607,   0.88729833460000,   0.00635083269197));
            mPoints.push_back(ChastePoint<3>(0.10000000001607,   0.11270166540000,   0.69856850110968));
            mPoints.push_back(ChastePoint<3>(0.44364916730000,   0.11270166540000,   0.39364916729197));
            mPoints.push_back(ChastePoint<3>(0.78729833458393,   0.11270166540000,   0.08872983347426));
            mPoints.push_back(ChastePoint<3>(0.05635083270000,   0.50000000000000,   0.39364916729197));
            mPoints.push_back(ChastePoint<3>(0.25000000000000,   0.50000000000000,   0.22182458365000));
            mPoints.push_back(ChastePoint<3>(0.44364916730000,   0.50000000000000,   0.05000000000803));
            mPoints.push_back(ChastePoint<3>(0.01270166538393,   0.88729833460000,   0.08872983347426));
            mPoints.push_back(ChastePoint<3>(0.05635083270000,   0.88729833460000,   0.05000000000803));
            mPoints.push_back(ChastePoint<3>(0.10000000001607,   0.88729833460000,   0.01127016654181));
            break;

        default:
            EXCEPTION("Gauss quadrature order not supported.");
    }
    assert(mPoints.size() == mWeights.size());
    mNumQuadPoints = mPoints.size();
}

template<unsigned ELEMENT_DIM>
GaussianQuadratureRule<ELEMENT_DIM>::GaussianQuadratureRule(unsigned numPointsInEachDimension, unsigned quadratureOrder)
{
    EXCEPTION("Gauss quadrature rule not available for this dimension.");
}

//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class GaussianQuadratureRule<0>;
template class GaussianQuadratureRule<1>;
template class GaussianQuadratureRule<2>;
template class GaussianQuadratureRule<3>;
template class GaussianQuadratureRule<4>;
