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

#include "DeltaNotchOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

DeltaNotchOdeSystem::DeltaNotchOdeSystem(double meanDelta, std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<DeltaNotchOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell
     * 1 - Delta concentration for this cell
     * 2 - average Delta concentration for this cell's immediate neighbours
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1.0); // soon overwritten
    SetDefaultInitialCondition(1, 1.0); // soon overwritten
    SetDefaultInitialCondition(2, 0.5);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

DeltaNotchOdeSystem::~DeltaNotchOdeSystem()
{
}

void DeltaNotchOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double notch = rY[0];
    double delta = rY[1];
    double mean_delta = rY[2];


    // The next two lines define the ODE system by Collier et al. (1996)
    double dx1 = mean_delta*mean_delta/(0.01 + mean_delta*mean_delta) - notch;
    double dx2 = 1/(1 + 100*notch*notch) - delta;

    rDY[0] = dx1;
    rDY[1] = dx2;
    rDY[2] = 0.0; // don't change the mean Delta level over the course of the mechanics time step
}

template<>
void CellwiseOdeSystemInformation<DeltaNotchOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Mean Delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchOdeSystem)
