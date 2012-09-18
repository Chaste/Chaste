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
 * @param numPointsInEachDimension  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly (ignored in 0-d case)
 */
template<>
GaussianQuadratureRule<0>::GaussianQuadratureRule(unsigned numPointsInEachDimension, unsigned quadratureOrder)
{
    mNumQuadPoints = 1;
    mWeights.push_back(1);
    mPoints.push_back(ChastePoint<0>());
}

/**
 * Constructor specialization for 1d.
 *
 * @param numPointsInEachDimension  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 */
template<>
GaussianQuadratureRule<1>::GaussianQuadratureRule(unsigned numPointsInEachDimension, unsigned quadratureOrder)
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
 * @param numPointsInEachDimension  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 *
 */
template<>
GaussianQuadratureRule<2>::GaussianQuadratureRule(unsigned numPointsInEachDimension, unsigned quadratureOrder)
{
    switch (quadratureOrder)
    {
        case 0: // 2d, 0th order
        case 1: // 2d, 1st order
        	//This is now order 1
            mNumQuadPoints = 1;
            {
                double one_third = 1.0/3.0;
                mWeights.push_back(0.5);
                mPoints.push_back(ChastePoint<2>(one_third, one_third));
            }
            break;

        case 2: // 2d, 2nd order
        	mNumQuadPoints = 4;
            mWeights.push_back(0.19716878364870);
            mWeights.push_back(0.19716878364870);
            mWeights.push_back(0.05283121635130);
            mWeights.push_back(0.05283121635130);

            mPoints.push_back(ChastePoint<2>(0.16666666666667,0.21132486540519));
            mPoints.push_back(ChastePoint<2>(0.62200846792815,0.21132486540519));
            mPoints.push_back(ChastePoint<2>(0.04465819873852,0.78867513459481));
            mPoints.push_back(ChastePoint<2>(0.16666666666667,0.78867513459481));
            break;

        case 3:
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
///// these should work but aren't tested yet
//        case 4: // 2d, 4 points per dimension
//            mWeights.push_back(0.0281503507763256);
//            mWeights.push_back(0.0527752633353762);
//            mWeights.push_back(0.0527752633353762);
//            mWeights.push_back(0.0281503507763256);
//            mWeights.push_back(0.0379971374528828);
//            mWeights.push_back(0.0712356642018415);
//            mWeights.push_back(0.0712356642018415);
//            mWeights.push_back(0.0379971374528828);
//            mWeights.push_back(0.0187158102927055);
//            mWeights.push_back(0.0350877267775703);
//            mWeights.push_back(0.0350877267775703);
//            mWeights.push_back(0.0187158102927055);
//            mWeights.push_back(0.00210036275308617);
//            mWeights.push_back(0.00393768441021206);
//            mWeights.push_back(0.00393768441021206);
//            mWeights.push_back(0.00210036275308617);
//            mPoints.push_back(ChastePoint<2>(0.0646110610581461,   0.0694318417));
//            mPoints.push_back(ChastePoint<2>(0.307096312350118,    0.0694318417));
//            mPoints.push_back(ChastePoint<2>(0.623471845949882,    0.0694318417));
//            mPoints.push_back(ChastePoint<2>(0.865957097241854,    0.0694318417));
//            mPoints.push_back(ChastePoint<2>(0.046518675850118,    0.3300094782));
//            mPoints.push_back(ChastePoint<2>(0.221103222498164,    0.3300094782));
//            mPoints.push_back(ChastePoint<2>(0.448887299301836,    0.3300094782));
//            mPoints.push_back(ChastePoint<2>(0.623471845949882,    0.3300094782));
//            mPoints.push_back(ChastePoint<2>(0.022913165849882,    0.6699905218));
//            mPoints.push_back(ChastePoint<2>(0.108906255701836,    0.6699905218));
//            mPoints.push_back(ChastePoint<2>(0.221103222498164,    0.6699905218));
//            mPoints.push_back(ChastePoint<2>(0.307096312350118,    0.6699905218));
//            mPoints.push_back(ChastePoint<2>(0.00482078064185386,  0.9305681583));
//            mPoints.push_back(ChastePoint<2>(0.022913165849882,    0.9305681583));
//            mPoints.push_back(ChastePoint<2>(0.046518675850118,    0.9305681583));
//            mPoints.push_back(ChastePoint<2>(0.0646110610581462,   0.9305681583));
//            break;
//        case 5:
//            mWeights.push_back(0.013375270551691);
//            mWeights.push_back(0.0270200993092602);
//            mWeights.push_back(0.0321155735571689);
//            mWeights.push_back(0.0270200993092602);
//            mWeights.push_back(0.013375270551691);
//            mWeights.push_back(0.0218078024655245);
//            mWeights.push_back(0.0440551079739245);
//            mWeights.push_back(0.0523630592364514);
//            mWeights.push_back(0.0440551079739245);
//            mWeights.push_back(0.0218078024655245);
//            mWeights.push_back(0.0168481340447735);
//            mWeights.push_back(0.0340358165695537);
//            mWeights.push_back(0.0404543209892346);
//            mWeights.push_back(0.0340358165695537);
//            mWeights.push_back(0.0168481340447735);
//            mWeights.push_back(0.00654219752778963);
//            mWeights.push_back(0.0132162430822249);
//            mWeights.push_back(0.0157085739026559);
//            mWeights.push_back(0.0132162430822249);
//            mWeights.push_back(0.00654219752778963);
//            mWeights.push_back(0.000658316657259778);
//            mWeights.push_back(0.00132990068405387);
//            mWeights.push_back(0.00158069453237811);
//            mWeights.push_back(0.00132990068405387);
//            mWeights.push_back(0.000658316657259778);
//            mPoints.push_back(ChastePoint<2>(0.0447095217211631,  0.04691007705));
//            mPoints.push_back(ChastePoint<2>( 0.219940124837926,  0.04691007705));
//            mPoints.push_back(ChastePoint<2>(    0.476544961475,  0.04691007705));
//            mPoints.push_back(ChastePoint<2>( 0.733149798112074,  0.04691007705));
//            mPoints.push_back(ChastePoint<2>( 0.908380401228837,  0.04691007705));
//            mPoints.push_back(ChastePoint<2>(0.0360848569379257,  0.23076534495));
//            mPoints.push_back(ChastePoint<2>( 0.177512700520108,  0.23076534495));
//            mPoints.push_back(ChastePoint<2>(    0.384617327525,  0.23076534495));
//            mPoints.push_back(ChastePoint<2>( 0.591721954529893,  0.23076534495));
//            mPoints.push_back(ChastePoint<2>( 0.733149798112074,  0.23076534495));
//            mPoints.push_back(ChastePoint<2>(    0.023455038525,            0.5));
//            mPoints.push_back(ChastePoint<2>(    0.115382672475,            0.5));
//            mPoints.push_back(ChastePoint<2>(              0.25,            0.5));
//            mPoints.push_back(ChastePoint<2>(    0.384617327525,            0.5));
//            mPoints.push_back(ChastePoint<2>(    0.476544961475,            0.5));
//            mPoints.push_back(ChastePoint<2>(0.0108252201120743,  0.76923465505));
//            mPoints.push_back(ChastePoint<2>(0.0532526444298925,  0.76923465505));
//            mPoints.push_back(ChastePoint<2>(    0.115382672475,  0.76923465505));
//            mPoints.push_back(ChastePoint<2>( 0.177512700520108,  0.76923465505));
//            mPoints.push_back(ChastePoint<2>( 0.219940124837926,  0.76923465505));
//            mPoints.push_back(ChastePoint<2>(0.00220055532883694,  0.95308992295));
//            mPoints.push_back(ChastePoint<2>(0.0108252201120743,  0.95308992295));
//            mPoints.push_back(ChastePoint<2>(    0.023455038525,  0.95308992295));
//            mPoints.push_back(ChastePoint<2>(0.0360848569379257,  0.95308992295));
//            mPoints.push_back(ChastePoint<2>(0.0447095217211631,  0.95308992295));
//            break;
//        case 8:
//            mWeights.push_back(0.00251093933534381);
//            mWeights.push_back(0.00551608575378066);
//            mWeights.push_back(0.0077813864127667);
//            mWeights.push_back(0.0089962476127433);
//            mWeights.push_back(0.0089962476127433);
//            mWeights.push_back(0.0077813864127667);
//            mWeights.push_back(0.00551608575378066);
//            mWeights.push_back(0.00251093933534381);
//            mWeights.push_back(0.00505566374657279);
//            mWeights.push_back(0.0111063912918299);
//            mWeights.push_back(0.015667472579425);
//            mWeights.push_back(0.018113541124127);
//            mWeights.push_back(0.018113541124127);
//            mWeights.push_back(0.015667472579425);
//            mWeights.push_back(0.0111063912918299);
//            mWeights.push_back(0.00505566374657279);
//            mWeights.push_back(0.00605561321825424);
//            mWeights.push_back(0.0133031018843967);
//            mWeights.push_back(0.018766310182895);
//            mWeights.push_back(0.0216961816606203);
//            mWeights.push_back(0.0216961816606203);
//            mWeights.push_back(0.018766310182895);
//            mWeights.push_back(0.0133031018843967);
//            mWeights.push_back(0.00605561321825424);
//            mWeights.push_back(0.00543106981966284);
//            mWeights.push_back(0.0119310914598135);
//            mWeights.push_back(0.0168308538189853);
//            mWeights.push_back(0.0194585541004693);
//            mWeights.push_back(0.0194585541004693);
//            mWeights.push_back(0.0168308538189853);
//            mWeights.push_back(0.0119310914598135);
//            mWeights.push_back(0.00543106981966284);
//            mWeights.push_back(0.00374741731366922);
//            mWeights.push_back(0.00823240727740299);
//            mWeights.push_back(0.0116132244841987);
//            mWeights.push_back(0.0134263275848652);
//            mWeights.push_back(0.0134263275848652);
//            mWeights.push_back(0.0116132244841987);
//            mWeights.push_back(0.00823240727740299);
//            mWeights.push_back(0.00374741731366922);
//            mWeights.push_back(0.00188340292975561);
//            mWeights.push_back(0.00413750022679507);
//            mWeights.push_back(0.00583665473756204);
//            mWeights.push_back(0.00674789664256371);
//            mWeights.push_back(0.00674789664256371);
//            mWeights.push_back(0.00583665473756204);
//            mWeights.push_back(0.00413750022679507);
//            mWeights.push_back(0.00188340292975561);
//            mWeights.push_back(0.000572162909255913);
//            mWeights.push_back(0.00125693983449269);
//            mWeights.push_back(0.00177312953176681);
//            mWeights.push_back(0.00204995761308944);
//            mWeights.push_back(0.00204995761308944);
//            mWeights.push_back(0.00177312953176681);
//            mWeights.push_back(0.00125693983449269);
//            mWeights.push_back(0.000572162909255913);
//            mWeights.push_back(5.0864805016297e-05);
//            mWeights.push_back(0.000111740902048041);
//            mWeights.push_back(0.000157629735243144);
//            mWeights.push_back(0.00018223952058876);
//            mWeights.push_back(0.00018223952058876);
//            mWeights.push_back(0.000157629735243144);
//            mWeights.push_back(0.000111740902048041);
//            mWeights.push_back(5.0864805016297e-05);
//            mPoints.push_back(ChastePoint<2>(0.0194608478758024,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.0996481604597984,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.232523501027757,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.400176196869137,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.579968731380863,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.747621427222243,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.880496767790202,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.960684080374198,  0.01985507175));
//            mPoints.push_back(ChastePoint<2>(0.0178364709097984,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.0913306309467688,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.213115003436359,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.366773901106599,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.531559337593401,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.685218235263642,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.807002607753231,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.880496767790202,   0.1016667613));
//            mPoints.push_back(ChastePoint<2>(0.0151447777277575,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.0775479696863585,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.180953921536175,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.311424229416958,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.451341975533043,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.581812283413825,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.685218235263642,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.747621427222243,  0.23723379505));
//            mPoints.push_back(ChastePoint<2>(0.0117485898691365,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.0601579836565992,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.140375345716958,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.241587932982724,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.350129388267276,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.451341975533043,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.531559337593401,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.579968731380863,  0.40828267875));
//            mPoints.push_back(ChastePoint<2>(0.00810648188086345,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.0415087776434008,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.0968584493330425,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.166694745767276,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.241587932982724,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.311424229416958,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.366773901106599,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.400176196869137,  0.59171732125));
//            mPoints.push_back(ChastePoint<2>(0.00471029402224254,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.0241187916136415,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.0562798735138254,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.0968584493330425,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.140375345716957,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.180953921536175,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.213115003436359,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.232523501027757,  0.76276620495));
//            mPoints.push_back(ChastePoint<2>(0.00201860084020162,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.0103361303532312,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.0241187916136415,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.0415087776434008,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.0601579836565991,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.0775479696863585,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.0913306309467688,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.0996481604597983,   0.8983332387));
//            mPoints.push_back(ChastePoint<2>(0.000394223874197648,  0.98014492825));
//            mPoints.push_back(ChastePoint<2>(0.00201860084020162,  0.98014492825));
//            mPoints.push_back(ChastePoint<2>(0.00471029402224255,  0.98014492825));
//            mPoints.push_back(ChastePoint<2>(0.00810648188086345,  0.98014492825));
//            mPoints.push_back(ChastePoint<2>(0.0117485898691366,  0.98014492825));
//            mPoints.push_back(ChastePoint<2>(0.0151447777277575,  0.98014492825));
//            mPoints.push_back(ChastePoint<2>(0.0178364709097984,  0.98014492825));
//            mPoints.push_back(ChastePoint<2>(0.0194608478758024,  0.98014492825));
//            break;
        default:
            EXCEPTION("Gauss quadrature order not supported.");
    }
    assert(mPoints.size() == mWeights.size());
    mNumQuadPoints = mPoints.size();
}

/**
 * Constructor specialization for 3d.
 *
 * @param numPointsInEachDimension  deprecated
 * @param quadratureOrder The minimum polynomial order that the rule can integrate exactly
 */
template<>
GaussianQuadratureRule<3>::GaussianQuadratureRule(unsigned numPointsInEachDimension, unsigned quadratureOrder)
{
//    mNumQuadPoints = 4;

//    double a = 0.585410196624969;
//    double b = 0.138196601125011;
//    double w = .0416666666666666666666666666666666666666666667;
//
//
//    mWeights.push_back(w);
//    mWeights.push_back(w);
//    mWeights.push_back(w);
//    mWeights.push_back(w);
//
//    mPoints.push_back(ChastePoint<3>(a,b,b));
//    mPoints.push_back(ChastePoint<3>(b,a,b));
//    mPoints.push_back(ChastePoint<3>(b,b,a));
//    mPoints.push_back(ChastePoint<3>(b,b,b));
//
//
//return;
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

        case 2:
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

// Should be a valid quadrature rule, but not tested and the order of the rule is unknown.
//        case 4: //Doesn't appear to be order 4
//            mWeights.push_back(0.00423982561968);
//            mWeights.push_back(0.00572288385156);
//            mWeights.push_back(0.00281885467361);
//            mWeights.push_back(0.00031634320391);
//            mWeights.push_back(0.00412036229051);
//            mWeights.push_back(0.00556163317318);
//            mWeights.push_back(0.00273942929295);
//            mWeights.push_back(0.00030742976838);
//            mWeights.push_back(0.00099965677330);
//            mWeights.push_back(0.00134932898618);
//            mWeights.push_back(0.00066462336430);
//            mWeights.push_back(0.00007458670588);
//            mWeights.push_back(0.00002360309872);
//            mWeights.push_back(0.00003185928022);
//            mWeights.push_back(0.00001569255698);
//            mWeights.push_back(0.00000176108183);
//            mWeights.push_back(0.00794866986669);
//            mWeights.push_back(0.01072905315027);
//            mWeights.push_back(0.00528468555374);
//            mWeights.push_back(0.00059306865848);
//            mWeights.push_back(0.00772470439029);
//            mWeights.push_back(0.01042674628127);
//            mWeights.push_back(0.00513578175757);
//            mWeights.push_back(0.00057635807584);
//            mWeights.push_back(0.00187411992466);
//            mWeights.push_back(0.00252967258912);
//            mWeights.push_back(0.00124601155388);
//            mWeights.push_back(0.00013983242583);
//            mWeights.push_back(0.00004425022545);
//            mWeights.push_back(0.00005972861231);
//            mWeights.push_back(0.00002941983138);
//            mWeights.push_back(0.00000330161175);
//            mWeights.push_back(0.00794866986669);
//            mWeights.push_back(0.01072905315027);
//            mWeights.push_back(0.00528468555374);
//            mWeights.push_back(0.00059306865848);
//            mWeights.push_back(0.00772470439029);
//            mWeights.push_back(0.01042674628127);
//            mWeights.push_back(0.00513578175757);
//            mWeights.push_back(0.00057635807584);
//            mWeights.push_back(0.00187411992466);
//            mWeights.push_back(0.00252967258912);
//            mWeights.push_back(0.00124601155388);
//            mWeights.push_back(0.00013983242583);
//            mWeights.push_back(0.00004425022545);
//            mWeights.push_back(0.00005972861231);
//            mWeights.push_back(0.00002941983138);
//            mWeights.push_back(0.00000330161175);
//            mWeights.push_back(0.00423982561968);
//            mWeights.push_back(0.00572288385156);
//            mWeights.push_back(0.00281885467361);
//            mWeights.push_back(0.00031634320391);
//            mWeights.push_back(0.00412036229051);
//            mWeights.push_back(0.00556163317318);
//            mWeights.push_back(0.00273942929295);
//            mWeights.push_back(0.00030742976838);
//            mWeights.push_back(0.00099965677330);
//            mWeights.push_back(0.00134932898618);
//            mWeights.push_back(0.00066462336430);
//            mWeights.push_back(0.00007458670588);
//            mWeights.push_back(0.00002360309872);
//            mWeights.push_back(0.00003185928022);
//            mWeights.push_back(0.00001569255698);
//            mWeights.push_back(0.00000176108183);
//
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.06012499793653));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.04328879995478));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.02132226325621));
//            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.00448606527446));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.04328879995478));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.03116707302848));
//            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.01535160449661));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.00322987757031));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.02132226325621));
//            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.01535160449661));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.00756156217830));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.00159090341870));
//            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.00448606527446));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.00322987757031));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.00159090341870));
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00033471571455));
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.28577404826889));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.20575161800155));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.10134469352354));
//            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.02132226325621));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.20575161800155));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.14813706341321));
//            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.07296615908496));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.01535160449661));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.10134469352354));
//            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.07296615908496));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.03594009661688));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.00756156217830));
//            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.02132226325621));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.01535160449661));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.00756156217830));
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00159090341870));
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.58018304432012));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.41772022627335));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.20575161800155));
//            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.04328879995478));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.41772022627335));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.30075023588863));
//            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.14813706341321));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.03116707302848));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.20575161800155));
//            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.14813706341321));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.07296615908496));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.01535160449661));
//            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.04328879995478));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.03116707302848));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.01535160449661));
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00322987757031));
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.80583209465249));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.58018304432012));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.28577404826889));
//            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.06012499793653));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.58018304432012));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.41772022627335));
//            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.20575161800155));
//            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.04328879995478));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.28577404826889));
//            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.20575161800155));
//            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.10134469352354));
//            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.02132226325621));
//            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.06012499793653));
//            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.04328879995478));
//            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.02132226325621));
//            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00448606527446));
//            break;

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
