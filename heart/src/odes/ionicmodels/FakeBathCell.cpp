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

#include "FakeBathCell.hpp"
#include "OdeSystemInformation.hpp"

FakeBathCell::FakeBathCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                           boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCell(pSolver, 1, 0, pIntracellularStimulus)
{
    mpSystemInfo = OdeSystemInformation<FakeBathCell>::Instance();
    Init();
}

FakeBathCell::~FakeBathCell()
{
}

// LCOV_EXCL_START  // This method should never be called (it would be a waste of time!)
void FakeBathCell::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

double FakeBathCell::GetIIonic(const std::vector<double>* pStateVariables)
{
    return 0.0;
}

void FakeBathCell::ComputeExceptVoltage(double tStart, double tEnd)
{
}

double FakeBathCell::GetIntracellularCalciumConcentration()
{
    return 0.0;
}

template<>
void OdeSystemInformation<FakeBathCell>::Initialise(void)
{
    // State variables
    this->mVariableNames.push_back("Fake voltage");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FakeBathCell)
