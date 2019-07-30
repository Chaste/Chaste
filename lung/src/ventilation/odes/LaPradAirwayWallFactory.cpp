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


#include "LaPradAirwayWallFactory.hpp"

//Airway generation alpha0 data, generation 0 -> generation 16
const double LaPradAirwayWallFactory::mAlpha0[] = {0.882,0.882,0.686,0.546,0.450,0.370,0.310,0.255,0.213,0.184,0.153,0.125,0.100,0.075,0.057,0.045,0.039};

const double LaPradAirwayWallFactory::mMaxGeneration = 16;

//The below are (very) approximate values from Fig. 2b of LaPrad 2013 - these are in kPa, but will be converted to Pa in code
const double LaPradAirwayWallFactory::mk1[] = {4000.,3850.,3700.,3550.,3400.,3250.,3100.,2950.,2800.,2650.,2500.,2350.,2200.,2050.,1900.,1750.,1600.};

const double LaPradAirwayWallFactory::mk2[] = {1000.,960.,920.,880.,840.,800.,760.,720.,680.,640.,600.,560.,520.,480.,440.,400.,360.};

const double LaPradAirwayWallFactory::mk3[] = {20.,19.4,18.8,18.2,17.6,17.,16.4,15.8,15.2,14.6,14.,13.4,12.8,12.2,11.6,11.,10.4};

LaPradAirwayWallFactory::LaPradAirwayWallFactory(bool useStrahlerOrder) : mpWalker(nullptr), mUseStrahlerOrder(useStrahlerOrder)
{}

LaPradAirwayWallFactory::~LaPradAirwayWallFactory()
{
    if (mpWalker)
    {
        delete mpWalker;
    }
}

double LaPradAirwayWallFactory::Getk1ForGeneration(unsigned generation)
{
    return mk1[generation];
}

double LaPradAirwayWallFactory::Getk2ForGeneration(unsigned generation)
{
    return mk2[generation];
}

double LaPradAirwayWallFactory::Getk3ForGeneration(unsigned generation)
{
    return mk3[generation];
}
double LaPradAirwayWallFactory::GetAlpha0ForGeneration(unsigned generation)
{
    return mAlpha0[generation];
}


LaPradAirwayWall* LaPradAirwayWallFactory::CreateAirwayWallForElement(Element<1,3>* pElement)
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

    double k1 = mk1[lower_generation] + (mk1[upper_generation] - mk1[lower_generation] )*(generation_factor - lower_generation);
    double k2 = mk2[lower_generation] + (mk2[upper_generation] - mk2[lower_generation] )*(generation_factor - lower_generation);
    double k3 = mk3[lower_generation] + (mk3[upper_generation] - mk3[lower_generation] )*(generation_factor - lower_generation);


    LaPradAirwayWall* wall = new LaPradAirwayWall;
    wall->Setk1(k1);
    wall->Setk2(k2);
    wall->Setk3(k3);
    wall->SetRIn(RIn);
    wall->SetROut(ROut);

    return wall;
}

double LaPradAirwayWallFactory::GetPleuralPressureForAirway(double time, Element<1,3>* pElement)
{
    return 0.0;
}

/**
* @param pMesh  the mesh for which to create acinar units.
*/
void LaPradAirwayWallFactory::SetMesh(AbstractTetrahedralMesh<1,3>* pMesh)
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
