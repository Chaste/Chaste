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


#include "LambertAirwayWallFactory.hpp"

//Airway generation alpha0 data, generation 0 -> generation 16
const double LambertAirwayWallFactory::mAlpha0[] = {0.882,0.882,0.686,0.546,0.450,0.370,0.310,0.255,0.213,0.184,0.153,0.125,0.100,0.075,0.057,0.045,0.039};

//The below are in Pascals^-1
const double LambertAirwayWallFactory::mAlpha0Prime[] = {0.011/98.0665,0.011/98.0665,0.051/98.0665,0.080/98.0665,0.100/98.0665,0.125/98.0665,0.142/98.0665,0.159/98.0665,0.174/98.0665,0.184/98.0665,0.194/98.0665,0.206/98.0665,0.218/98.0665,0.226/98.0665,0.233/98.0665,0.239/98.0665,0.243/98.0665};

const double LambertAirwayWallFactory::mN1[] = {0.5,0.5,0.6,0.6,0.7,0.8,0.9,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

const double LambertAirwayWallFactory::mN2[] = {10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,9.0,8.0,8.0,8.0,7.0,7.0};

const double LambertAirwayWallFactory::mMaxGeneration = 16;

LambertAirwayWallFactory::LambertAirwayWallFactory(bool useStrahlerOrder) : mpWalker(nullptr), mUseStrahlerOrder(useStrahlerOrder)
{}

LambertAirwayWallFactory::~LambertAirwayWallFactory()
{
    if (mpWalker)
    {
        delete mpWalker;
    }
}

double LambertAirwayWallFactory::GetAlpha0ForGeneration(unsigned generation)
{
    return mAlpha0[generation];
}

double LambertAirwayWallFactory::GetAlpha0PrimeForGeneration(unsigned generation)
{
    return mAlpha0Prime[generation];
}

double LambertAirwayWallFactory::GetN1ForGeneration(unsigned generation)
{
    return mN1[generation];
}

double LambertAirwayWallFactory::GetN2ForGeneration(unsigned generation)
{
    return mN2[generation];
}

LambertAirwayWall* LambertAirwayWallFactory::CreateAirwayWallForElement(Element<1,3>* pElement)
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

    double RiMax = pElement->GetAttribute();
    double alpha0 = mAlpha0[lower_generation] + (mAlpha0[upper_generation] - mAlpha0[lower_generation] )*(generation_factor - lower_generation);
    double alpha0prime = mAlpha0Prime[lower_generation] + (mAlpha0Prime[upper_generation] - mAlpha0Prime[lower_generation] )*(generation_factor - lower_generation);
    double N1 = mN1[lower_generation] + (mN1[upper_generation] - mN1[lower_generation] )*(generation_factor - lower_generation);
    double N2 = mN2[lower_generation] + (mN2[upper_generation] - mN2[lower_generation] )*(generation_factor - lower_generation);
    double P1 = alpha0*N1/alpha0prime;

    double lumen_max_area = M_PI*RiMax*RiMax;
    double Ri = std::sqrt(alpha0*lumen_max_area/M_PI);

    LambertAirwayWall* wall = new LambertAirwayWall;
    wall->SetN1(N1);
    wall->SetN2(N2);
    wall->SetP1(P1);
    wall->SetRiMax(RiMax);
    wall->SetRi(Ri);

    return wall;
}

double LambertAirwayWallFactory::GetPleuralPressureForAirway(double time, Element<1,3>* pElement)
{
    return 0.0;
}

/**
* @param pMesh  the mesh for which to create acinar units.
*/
void LambertAirwayWallFactory::SetMesh(AbstractTetrahedralMesh<1,3>* pMesh)
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
