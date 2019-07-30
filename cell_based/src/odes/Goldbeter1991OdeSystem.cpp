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

#include "Goldbeter1991OdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

Goldbeter1991OdeSystem::Goldbeter1991OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Goldbeter1991OdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - C (Cyclin)
     * 1 - M (CDC-2 Kinase)
     * 2 - X (Cyclin Protease)
     */

//    SetDefaultInitialCondition(0, 0.01); // soon overwritten
//    SetDefaultInitialCondition(1, 0.01); // soon overwritten
//    SetDefaultInitialCondition(2, 0.01); // soon overwritten

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Goldbeter1991OdeSystem::~Goldbeter1991OdeSystem()
{
}

void Goldbeter1991OdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    // state values
    double C = rY[0]; // cyclin
    double M = rY[1]; // kinase
    double X = rY[2]; // protease

    // consts
    const double cell = 1;
    const double VM1 = 3;
    const double VM3 = 1;
    const double Kc = 0.5;

    // calculations
    double reaction1 = cell * 0.025;
    double reaction2 = C * cell * 0.01;
    double reaction3 = C * cell * 0.25 * X * pow(C + 0.02, -1);
    double reaction5 = cell * M * 1.5 * pow(0.005 + M, -1);
    double reaction7 = cell * 0.5 * X * pow(0.005 + X, -1);
    double V3 = M * VM3;
    double V1 = C * VM1 * pow(C + Kc, -1);
    double reaction6 = cell * V3 * (1 + -1 * X) * pow(0.005 + -1 * X + 1, -1);
    double reaction4 = cell * (1 + -1 * M) * V1 * pow(0.005 + -1 * M + 1, -1);

    // ODEs
    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // dC/dt
    rDY[1] = (reaction4 - reaction5) / cell; // dM/dt
    rDY[2] = (reaction6 - reaction7) / cell; // dX/dt
}

template<>
void CellwiseOdeSystemInformation<Goldbeter1991OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Cyclin");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01); // will be filled in later

    this->mVariableNames.push_back("CDC-2 Kinase");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01); // will be filled in later

    this->mVariableNames.push_back("Cyclin Protease");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01); // will be filled in later

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991OdeSystem)

