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


#include "HiornsAirwayWallFactory.hpp"

//Airway generation alpha0 data, generation 0 -> generation 16
const double HiornsAirwayWallFactory::mAlpha0[] = {0.882,0.882,0.686,0.546,0.450,0.370,0.310,0.255,0.213,0.184,0.153,0.125,0.100,0.075,0.057,0.045,0.039};

//The below are in Pascals^-1
const double HiornsAirwayWallFactory::mAlpha0Prime[] = {0.011/98.0665,0.011/98.0665,0.051/98.0665,0.080/98.0665,0.100/98.0665,0.125/98.0665,0.142/98.0665,0.159/98.0665,0.174/98.0665,0.184/98.0665,0.194/98.0665,0.206/98.0665,0.218/98.0665,0.226/98.0665,0.233/98.0665,0.239/98.0665,0.243/98.0665};

const double HiornsAirwayWallFactory::mN1[] = {0.5,0.5,0.6,0.6,0.7,0.8,0.9,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

const double HiornsAirwayWallFactory::mN2[] = {10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,9.0,8.0,8.0,8.0,7.0,7.0};

const double HiornsAirwayWallFactory::mMaxGeneration = 16;

//The below are values that can be returned for particular airway trees

const double HiornsAirwayWallFactory::mmu[] = {64002., 62826., 8099., 3472., 1617., 778., 425., 241., 170., 125., 94.3, 70.3, 52.9, 34., 31.3, 24.4, 18.73};

const double HiornsAirwayWallFactory::mphi1[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

const double HiornsAirwayWallFactory::mphi2[] = {0.7854, 0.7363, 0.6872, 0.6381, 0.5890, 0.54, 0.4909, 0.4418, 0.3927, 0.3436, 0.2945, 0.2454, 0.1963, 0.1473, 0.0982, 0.0491, 0.};

const double HiornsAirwayWallFactory::mC1[] = {179380., 176033., 8531., 2232., 930., 352., 157., 70.9, 32.9, 16.8, 8.9, 4.70, 2.56, 0.84, 0.83, 0.51, 0.299};

const double HiornsAirwayWallFactory::mC2[] = {101.9786, 102.3, 9.31, 2.72, 0.893, 0.4415, 0.2264, 0.1289, 0.07906, 0.04733, 0.02941, 0.018297, 0.0115175, 0.009006, 0.004454616, 0.002782845, 0.00165941};


HiornsAirwayWallFactory::HiornsAirwayWallFactory(bool useStrahlerOrder) : mpWalker(nullptr), mUseStrahlerOrder(useStrahlerOrder)
{}

HiornsAirwayWallFactory::~HiornsAirwayWallFactory()
{
    if (mpWalker)
    {
        delete mpWalker;
    }
}

double HiornsAirwayWallFactory::GetmuForGeneration(unsigned generation)
{
    return mmu[generation];
}
double HiornsAirwayWallFactory::Getphi1ForGeneration(unsigned generation)
{
    return mphi1[generation];
}
double HiornsAirwayWallFactory::Getphi2ForGeneration(unsigned generation)
{
    return mphi2[generation];
}
double HiornsAirwayWallFactory::GetC1ForGeneration(unsigned generation)
{
    return mC1[generation];
}
double HiornsAirwayWallFactory::GetC2ForGeneration(unsigned generation)
{
    return mC2[generation];
}

double HiornsAirwayWallFactory::GetAlpha0ForGeneration(unsigned generation)
{
    return mAlpha0[generation];
}

HiornsAirwayWall* HiornsAirwayWallFactory::CreateAirwayWallForElement(Element<1,3>* pElement)
{
    unsigned order = 0;

    if (mUseStrahlerOrder)
    {
        order = mpWalker->GetElementStrahlerOrder(pElement);
    }
    else
    {
        order = mpWalker->GetElementHorsfieldOrder(pElement);
    }


    //We linearly interpolate the generation data on to the corresponding order
    double generation_factor = (mMaxOrder - order)*mMaxGeneration/(mMaxOrder-1);
    unsigned lower_generation = std::floor(generation_factor);
    unsigned upper_generation = std::ceil(generation_factor);

    assert(lower_generation >= 0.0);
    assert(upper_generation <= 16.0);

    double RIn = pElement->GetAttribute()/2.;  // Measured R values are fully inflated.  As a proxy for uninflated, we simply divide by 2 (for now...)
    double ROut = 1.1*RIn;  // For now, assume ROut is 10% higher than RIn

    double mu = mmu[lower_generation] + (mmu[upper_generation] - mmu[lower_generation] )*(generation_factor - lower_generation);
    double phi1 = mphi1[lower_generation] + (mphi1[upper_generation] - mphi1[lower_generation] )*(generation_factor - lower_generation);
    double phi2 = mphi2[lower_generation] + (mphi2[upper_generation] - mphi2[lower_generation] )*(generation_factor - lower_generation);
    double C1 = mC1[lower_generation] + (mC1[upper_generation] - mC1[lower_generation] )*(generation_factor - lower_generation);
    double C2 = mC2[lower_generation] + (mC2[upper_generation] - mC2[lower_generation] )*(generation_factor - lower_generation);

    HiornsAirwayWall* wall = CreateBasicAirwayWall();

    wall->Setmu(mu);
    wall->Setphi1(phi1);
    wall->Setphi2(phi2);
    wall->SetC1(C1);
    wall->SetC2(C2);
    wall->SetRIn(RIn);
    wall->SetROut(ROut);

    return wall;
}

HiornsAirwayWall* HiornsAirwayWallFactory::CreateBasicAirwayWall()
{
    return new HiornsAirwayWall;
}


double HiornsAirwayWallFactory::GetPleuralPressureForAirway(double time, Element<1,3>* pElement)
{
    return 0.0;
}

/**
* @param pMesh  the mesh for which to create acinar units.
*/
void HiornsAirwayWallFactory::SetMesh(AbstractTetrahedralMesh<1,3>* pMesh)
{
    AbstractAirwayWallFactory::SetMesh(pMesh);

    mpWalker = new AirwayTreeWalker(*pMesh, 0u);

    if (mUseStrahlerOrder)
    {
        mMaxOrder = mpWalker->GetMaxElementStrahlerOrder();
    }
    else
    {
        mMaxOrder = mpWalker->GetMaxElementHorsfieldOrder();
    }
}
