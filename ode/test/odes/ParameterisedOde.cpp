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

#include "ParameterisedOde.hpp"
#include "OdeSystemInformation.hpp"

bool ParameterisedOde::fakeSecondParameter = false;
bool ParameterisedOde::noParameterDefaults = false;

ParameterisedOde::ParameterisedOde() : AbstractOdeSystem(1) // 1 here is the number of variables
{
    mpSystemInfo = OdeSystemInformation<ParameterisedOde>::Instance();
    ResetToInitialConditions();
    if (!noParameterDefaults)
    {
        mParameters.push_back(0);
        if (fakeSecondParameter)
        {
            mParameters.push_back(-1.0);
        }
    }
}

void ParameterisedOde::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    rDY[0] = mParameters[0];
}

std::vector<double> ParameterisedOde::ComputeDerivedQuantities(double time,
                                                               const std::vector<double>& rState)
{
    std::vector<double> dqs;
    dqs.push_back(2*mParameters[0] + rState[0]);
    return dqs;
}


#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ParameterisedOde)


template<>
void OdeSystemInformation<ParameterisedOde>::Initialise()
{
    this->mSystemName = "ParameterisedOde";

    this->mFreeVariableName = "time";
    this->mFreeVariableUnits = "ms";

    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mParameterNames.push_back("a");
    this->mParameterUnits.push_back("dimensionless");

    this->mDerivedQuantityNames.push_back("2a_plus_y");
    this->mDerivedQuantityUnits.push_back("dimensionless");

    this->mAttributes["attr"] = 1.1;

    this->mInitialised = true;
}
