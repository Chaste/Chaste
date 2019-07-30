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

#include "ParameterisedCvode.hpp"
#include "OdeSystemInformation.hpp"
#include "VectorHelperFunctions.hpp"

#ifdef CHASTE_CVODE
bool ParameterisedCvode::fakeSecondParameter = false;
bool ParameterisedCvode::noParameterDefaults = false;

ParameterisedCvode::ParameterisedCvode() : AbstractCvodeSystem(1) // 1 here is the number of variables
{
    this->mpSystemInfo = OdeSystemInformation<ParameterisedCvode>::Instance();
    Init();
    if (!noParameterDefaults)
    {
        CreateVectorIfEmpty(mParameters,1);
        SetVectorComponent(mParameters, 0, 0.0);
        if (fakeSecondParameter)
        {
            SetVectorComponent(mParameters, 1, -1.0);
        }
    }
}

void ParameterisedCvode::EvaluateYDerivatives(double time, const N_Vector rY, N_Vector rDY)
{
    NV_Ith_S(rDY, 0) = GetVectorComponent(mParameters,0);
}

N_Vector ParameterisedCvode::ComputeDerivedQuantities(double time,
                                                      const N_Vector& rState)
{
    N_Vector derived_quantities = NULL;
    CreateVectorIfEmpty(derived_quantities,1);
    SetVectorComponent(derived_quantities, 0,
                       2*GetVectorComponent(mParameters,0) + GetVectorComponent(rState,0));
    return derived_quantities;
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ParameterisedCvode)

template<>
void OdeSystemInformation<ParameterisedCvode>::Initialise()
{
    this->mSystemName = "ParameterisedCvode";

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
#endif // CHASTE_CVODE
