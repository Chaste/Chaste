/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "CellwiseOdeSystemInformation.hpp"
#include "DeltaNotchEdgeOdeSystem.hpp"

DeltaNotchEdgeOdeSystem::DeltaNotchEdgeOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(2)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<DeltaNotchEdgeOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell edge
     * 1 - Delta concentration for this cell edge
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1.0); // soon overwritten
    SetDefaultInitialCondition(1, 1.0); // soon overwritten

    this->mParameters.push_back(0.5);
    //By default zero. If no interior SRN model is specified, interior delta/notch is zero
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

DeltaNotchEdgeOdeSystem::~DeltaNotchEdgeOdeSystem()
{
}

void DeltaNotchEdgeOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    const double notch = rY[0];
    const double delta = rY[1];
    const double neigh_delta = this->mParameters[0]; // Shorthand for "this->mParameter("neighbor delta");"
    const double interior_delta = this->mParameters[1];
    const double interior_notch = this->mParameters[2];
    // The next two lines define the ODE system by Collier et al. (1996), except that we use neighbour Delta directly
    //That is, neighbour activates notch, which can be degraded and trafficked from cytoplasm to this edge
    rDY[0] = neigh_delta*neigh_delta/(0.01 + neigh_delta*neigh_delta) - notch+0.1*interior_notch;  // d[Notch]/dt
    // Delta is inhibited by notch in this this
    rDY[1] = 1.0/(1.0 + 100.0*notch*notch) - delta+0.1*interior_delta;// d[Delta]/dt
}

template<>
void CellwiseOdeSystemInformation<DeltaNotchEdgeOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mParameterNames.push_back("neighbour delta");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("interior delta");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("interior notch");
    this->mParameterUnits.push_back("non-dim");
    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchEdgeOdeSystem)
