/*

Copyright (C) University of Oxford, 2005-2011

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
 * @param numPointsInEachDimension  number of Gauss points in each dimension
 */
template<>
GaussianQuadratureRule<0>::GaussianQuadratureRule(unsigned numPointsInEachDimension)
{
    mNumQuadPoints = 1; // numPointsInEachDimension^0
    mWeights.reserve(mNumQuadPoints);
    mPoints.reserve(mNumQuadPoints);
    mWeights.push_back(1);
    mPoints.push_back(ChastePoint<0>());
}

/**
 * Constructor specialization for 1d.
 *
 * @param numPointsInEachDimension  number of Gauss points in each dimension
 */
template<>
GaussianQuadratureRule<1>::GaussianQuadratureRule(unsigned numPointsInEachDimension)
{
    mNumQuadPoints = numPointsInEachDimension;

    mWeights.reserve(mNumQuadPoints);
    mPoints.reserve(mNumQuadPoints);
    switch (numPointsInEachDimension)
    {
        case 1: // 1d, 1 point
            mWeights.push_back(1);
            mPoints.push_back(ChastePoint<1>(0.5));
            break;
        case 2: // 1d, 2 points
            mWeights.push_back(0.5);
            mWeights.push_back(0.5);

            mPoints.push_back(ChastePoint<1>(0.21132486540519));
            mPoints.push_back(ChastePoint<1>(0.78867513459481));
            break;
        case 3: // 1d, 3 points
            mWeights.push_back(5.0/18.0);
            mWeights.push_back(4.0/9.0);
            mWeights.push_back(5.0/18.0);
            mPoints.push_back(ChastePoint<1>(0.1127016654));
            mPoints.push_back(ChastePoint<1>(0.5));
            mPoints.push_back(ChastePoint<1>(0.8872983346));
            break;
///// these should work but aren't tested yet
//        case 4: // 1d, 4 points
//            mWeights.push_back(  0.17392732255);
//            mWeights.push_back(  0.32607267745);
//            mWeights.push_back(  0.32607267745);
//            mWeights.push_back(  0.17392732255);
//            mPoints.push_back(ChastePoint<1>(   0.0694318417));
//            mPoints.push_back(ChastePoint<1>(   0.3300094782));
//            mPoints.push_back(ChastePoint<1>(   0.6699905218));
//            mPoints.push_back(ChastePoint<1>(   0.9305681583));
//            break;
//        case 5: // 1d, 5 points
//            mWeights.push_back(   0.1184634425);
//            mWeights.push_back(  0.23931433525);
//            mWeights.push_back(  0.28444444445);
//            mWeights.push_back(  0.23931433525);
//            mWeights.push_back(   0.1184634425);
//            mPoints.push_back(ChastePoint<1>(  0.04691007705));
//            mPoints.push_back(ChastePoint<1>(  0.23076534495));
//            mPoints.push_back(ChastePoint<1>(            0.5));
//            mPoints.push_back(ChastePoint<1>(  0.76923465505));
//            mPoints.push_back(ChastePoint<1>(  0.95308992295));
//            break;
//        case 8: // 1d, 8 points
//            mWeights.push_back(  0.05061426815);
//            mWeights.push_back(  0.11119051725);
//            mWeights.push_back(  0.15685332295);
//            mWeights.push_back(   0.1813418917);
//            mWeights.push_back(   0.1813418917);
//            mWeights.push_back(  0.15685332295);
//            mWeights.push_back(  0.11119051725);
//            mWeights.push_back(  0.05061426815);
//            mPoints.push_back(ChastePoint<1>(  0.01985507175));
//            mPoints.push_back(ChastePoint<1>(   0.1016667613));
//            mPoints.push_back(ChastePoint<1>(  0.23723379505));
//            mPoints.push_back(ChastePoint<1>(  0.40828267875));
//            mPoints.push_back(ChastePoint<1>(  0.59171732125));
//            mPoints.push_back(ChastePoint<1>(  0.76276620495));
//            mPoints.push_back(ChastePoint<1>(   0.8983332387));
//            mPoints.push_back(ChastePoint<1>(  0.98014492825));
//            break;
         default:
            EXCEPTION("Number of gauss points per dimension not supported.");
    }
}

/**
 * Constructor specialization for 2d.
 *
 * @param numPointsInEachDimension  number of Gauss points in each dimension
 */
template<>
GaussianQuadratureRule<2>::GaussianQuadratureRule(unsigned numPointsInEachDimension)
{
    mNumQuadPoints = numPointsInEachDimension * numPointsInEachDimension;

    mWeights.reserve(mNumQuadPoints);
    mPoints.reserve(mNumQuadPoints);

    switch (numPointsInEachDimension)
    {
        case 1: // 2d, 1 point per dimension
            mWeights.push_back(0.5);
            mPoints.push_back(ChastePoint<2>(0.25,0.5));
            break;

        case 2: // 2d, 2 points per dimension
            mWeights.push_back(0.19716878364870);
            mWeights.push_back(0.19716878364870);
            mWeights.push_back(0.05283121635130);
            mWeights.push_back(0.05283121635130);

            mPoints.push_back(ChastePoint<2>(0.16666666666667,0.21132486540519));
            mPoints.push_back(ChastePoint<2>(0.62200846792815,0.21132486540519));
            mPoints.push_back(ChastePoint<2>(0.04465819873852,0.78867513459481));
            mPoints.push_back(ChastePoint<2>(0.16666666666667,0.78867513459481));
            break;

        case 3: // 2d, 3 points per dimension
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
            EXCEPTION("Number of gauss points per dimension not supported.");
    }
}

/**
 * Constructor specialization for 3d.
 *
 * @param numPointsInEachDimension  number of Gauss points in each dimension
 */
template<>
GaussianQuadratureRule<3>::GaussianQuadratureRule(unsigned numPointsInEachDimension)
{
    mNumQuadPoints = numPointsInEachDimension * numPointsInEachDimension * numPointsInEachDimension;

    mWeights.reserve(mNumQuadPoints);
    mPoints.reserve(mNumQuadPoints);

    switch (numPointsInEachDimension)
    {
        case 1: //3d, 1 point per dimension
            mWeights.push_back(0.12500000000000);
            mPoints.push_back(ChastePoint<3>(0.25000000000000,0.50000000000000,0.12500000000000));
            break;

        case 2: //3d, 2 points per dimension
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

        case 3: //3d, 3 points per dimension
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

        case 4: //3d, 4 points per dimension
            mWeights.push_back(0.00423982561968);
            mWeights.push_back(0.00572288385156);
            mWeights.push_back(0.00281885467361);
            mWeights.push_back(0.00031634320391);
            mWeights.push_back(0.00412036229051);
            mWeights.push_back(0.00556163317318);
            mWeights.push_back(0.00273942929295);
            mWeights.push_back(0.00030742976838);
            mWeights.push_back(0.00099965677330);
            mWeights.push_back(0.00134932898618);
            mWeights.push_back(0.00066462336430);
            mWeights.push_back(0.00007458670588);
            mWeights.push_back(0.00002360309872);
            mWeights.push_back(0.00003185928022);
            mWeights.push_back(0.00001569255698);
            mWeights.push_back(0.00000176108183);
            mWeights.push_back(0.00794866986669);
            mWeights.push_back(0.01072905315027);
            mWeights.push_back(0.00528468555374);
            mWeights.push_back(0.00059306865848);
            mWeights.push_back(0.00772470439029);
            mWeights.push_back(0.01042674628127);
            mWeights.push_back(0.00513578175757);
            mWeights.push_back(0.00057635807584);
            mWeights.push_back(0.00187411992466);
            mWeights.push_back(0.00252967258912);
            mWeights.push_back(0.00124601155388);
            mWeights.push_back(0.00013983242583);
            mWeights.push_back(0.00004425022545);
            mWeights.push_back(0.00005972861231);
            mWeights.push_back(0.00002941983138);
            mWeights.push_back(0.00000330161175);
            mWeights.push_back(0.00794866986669);
            mWeights.push_back(0.01072905315027);
            mWeights.push_back(0.00528468555374);
            mWeights.push_back(0.00059306865848);
            mWeights.push_back(0.00772470439029);
            mWeights.push_back(0.01042674628127);
            mWeights.push_back(0.00513578175757);
            mWeights.push_back(0.00057635807584);
            mWeights.push_back(0.00187411992466);
            mWeights.push_back(0.00252967258912);
            mWeights.push_back(0.00124601155388);
            mWeights.push_back(0.00013983242583);
            mWeights.push_back(0.00004425022545);
            mWeights.push_back(0.00005972861231);
            mWeights.push_back(0.00002941983138);
            mWeights.push_back(0.00000330161175);
            mWeights.push_back(0.00423982561968);
            mWeights.push_back(0.00572288385156);
            mWeights.push_back(0.00281885467361);
            mWeights.push_back(0.00031634320391);
            mWeights.push_back(0.00412036229051);
            mWeights.push_back(0.00556163317318);
            mWeights.push_back(0.00273942929295);
            mWeights.push_back(0.00030742976838);
            mWeights.push_back(0.00099965677330);
            mWeights.push_back(0.00134932898618);
            mWeights.push_back(0.00066462336430);
            mWeights.push_back(0.00007458670588);
            mWeights.push_back(0.00002360309872);
            mWeights.push_back(0.00003185928022);
            mWeights.push_back(0.00001569255698);
            mWeights.push_back(0.00000176108183);

            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.06012499793653));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.04328879995478));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.02132226325621));
            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.00448606527446));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.04328879995478));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.03116707302848));
            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.01535160449661));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.00322987757031));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.02132226325621));
            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.01535160449661));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.00756156217830));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.00159090341870));
            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.00448606527446));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.00322987757031));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.00159090341870));
            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00033471571455));
            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.28577404826889));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.20575161800155));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.10134469352354));
            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.02132226325621));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.20575161800155));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.14813706341321));
            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.07296615908496));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.01535160449661));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.10134469352354));
            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.07296615908496));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.03594009661688));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.00756156217830));
            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.02132226325621));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.01535160449661));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.00756156217830));
            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00159090341870));
            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.58018304432012));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.41772022627335));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.20575161800155));
            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.04328879995478));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.41772022627335));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.30075023588863));
            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.14813706341321));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.03116707302848));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.20575161800155));
            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.14813706341321));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.07296615908496));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.01535160449661));
            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.04328879995478));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.03116707302848));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.01535160449661));
            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00322987757031));
            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.06943184420000,   0.80583209465249));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.06943184420000,   0.58018304432012));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.06943184420000,   0.28577404826889));
            mPoints.push_back(ChastePoint<3>(0.86595709258901,   0.06943184420000,   0.06012499793653));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.33000947820000,   0.58018304432012));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.33000947820000,   0.41772022627335));
            mPoints.push_back(ChastePoint<3>(0.44888729930184,   0.33000947820000,   0.20575161800155));
            mPoints.push_back(ChastePoint<3>(0.62347184427491,   0.33000947820000,   0.04328879995478));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.66999052180000,   0.28577404826889));
            mPoints.push_back(ChastePoint<3>(0.10890625570184,   0.66999052180000,   0.20575161800155));
            mPoints.push_back(ChastePoint<3>(0.22110322249816,   0.66999052180000,   0.10134469352354));
            mPoints.push_back(ChastePoint<3>(0.30709631152509,   0.66999052180000,   0.02132226325621));
            mPoints.push_back(ChastePoint<3>(0.00482078098901,   0.93056815580000,   0.06012499793653));
            mPoints.push_back(ChastePoint<3>(0.02291316667491,   0.93056815580000,   0.04328879995478));
            mPoints.push_back(ChastePoint<3>(0.04651867752509,   0.93056815580000,   0.02132226325621));
            mPoints.push_back(ChastePoint<3>(0.06461106321099,   0.93056815580000,   0.00448606527446));
            break;

        default:
            EXCEPTION("Number of gauss points per dimension not supported.");
    }
}

template<unsigned ELEMENT_DIM>
GaussianQuadratureRule<ELEMENT_DIM>::GaussianQuadratureRule(unsigned numPointsInEachDimension)
{
    EXCEPTION("Gauss points not available for this dimension.");
}

//////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////

template class GaussianQuadratureRule<0>;
template class GaussianQuadratureRule<1>;
template class GaussianQuadratureRule<2>;
template class GaussianQuadratureRule<3>;
template class GaussianQuadratureRule<4>;
