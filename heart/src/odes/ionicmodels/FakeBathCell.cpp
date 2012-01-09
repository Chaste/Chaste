/*

Copyright (C) University of Oxford, 2005-2012

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

// This method should never be called, but we implement it with something sensible just in case...
#define COVERAGE_IGNORE
void FakeBathCell::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    rDY[0] = 0.0;
}
#undef COVERAGE_IGNORE

double FakeBathCell::GetIIonic(const std::vector<double>* pStateVariables)
{
    return 0.0;
}

void FakeBathCell::ComputeExceptVoltage(double tStart, double tEnd)
{
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
