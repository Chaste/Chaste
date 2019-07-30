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


#include "LambertAirwayWall.hpp"
#include "MathsCustomFunctions.hpp"
#include <cmath>

LambertAirwayWall::LambertAirwayWall()
{}

LambertAirwayWall::~LambertAirwayWall() {}

void LambertAirwayWall::SetTimestep(double dt) {}

void LambertAirwayWall::SolveAndUpdateState(double tStart, double tEnd)
{
    //NB. This isn't a time dependent model, so we just directly calculate the new radius
    double transpulmonary_pressure = mAirwayPressure - mPleuralPressure;
    double P2 = (mP1*mN2*(mRi*mRi - mRiMax*mRiMax))/(mN1*mRi*mRi);

    double ri2 = 0.0;

    if (transpulmonary_pressure <= 0.0)
    {
        ri2 = mRi*mRi*std::pow(1 - transpulmonary_pressure/mP1, -mN1);
        //std::cout << "here " << transpulmonary_pressure << " " << mP1 << " " <<  transpulmonary_pressure/mP1 << " " << mN1 << " " << ri2 << std::endl;;
    }
    else //transpulmonary_pressure > 0
    {
        ri2 = mRiMax*mRiMax - (mRiMax*mRiMax  - mRi*mRi)*std::pow(1 - transpulmonary_pressure/P2, -mN2);
    }

    mDeformedAirwayRadius = std::sqrt(ri2);
}

double LambertAirwayWall::GetLumenRadius()
{
    return mDeformedAirwayRadius;
}

void LambertAirwayWall::SetAirwayPressure(double pressure)
{
    mAirwayPressure = pressure;
}

void LambertAirwayWall::SetPleuralPressure(double pressure)
{
    mPleuralPressure = pressure;
}

void LambertAirwayWall::SetN1(double n1)
{
    mN1 = n1;
}

void LambertAirwayWall::SetN2(double n2)
{
    mN2 = n2;
}

void LambertAirwayWall::SetP1(double p1)
{
    mP1 = p1;
}

void LambertAirwayWall::SetRiMax(double riMax)
{
    this->mRiMax = riMax;
}

void LambertAirwayWall::SetRi(double Ri)
{
    mRi = Ri;
}
