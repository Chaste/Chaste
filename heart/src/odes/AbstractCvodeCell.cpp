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

#ifdef CHASTE_CVODE

#include <sstream>
#include <iostream>
#include <cmath>

#include "AbstractCvodeCell.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "VectorHelperFunctions.hpp"


AbstractCvodeCell::AbstractCvodeCell(boost::shared_ptr<AbstractIvpOdeSolver> /* unused */,
                                     unsigned numberOfStateVariables,
                                     unsigned voltageIndex,
                                     boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCellInterface(boost::shared_ptr<AbstractIvpOdeSolver>(), voltageIndex, pIntracellularStimulus),
      AbstractCvodeSystem(numberOfStateVariables),
      mMaxDt(DOUBLE_UNSET)
{
}


AbstractCvodeCell::~AbstractCvodeCell()
{
}

double AbstractCvodeCell::GetVoltage()
{
    assert(mStateVariables);
    return AbstractCvodeSystem::GetAnyVariable(mVoltageIndex);
}

void AbstractCvodeCell::SetVoltage(double voltage)
{
    assert(mStateVariables);
    SetAnyVariable(mVoltageIndex, voltage);
    SetFixedVoltage(voltage);
}


void AbstractCvodeCell::SetMaxTimestep(double maxDt)
{
    SetTimestep(maxDt); // Note: SetTimestep() sets the maximum timestep.
}

void AbstractCvodeCell::SetTimestep(double maxDt)
{
    mMaxDt = maxDt;
}

double AbstractCvodeCell::GetTimestep()
{
    return mMaxDt;
}


void AbstractCvodeCell::SolveAndUpdateState(double tStart, double tEnd)
{
    if (mMaxDt == DOUBLE_UNSET)
    {
        SetTimestep(HeartConfig::Instance()->GetPrintingTimeStep());
    }
    Solve(tStart, tEnd, mMaxDt);
}

OdeSolution AbstractCvodeCell::Compute(double tStart, double tEnd, double tSamp)
{
    if (tSamp == 0.0)
    {
        tSamp = HeartConfig::Instance()->GetPrintingTimeStep();
    }
    if (mMaxDt == DOUBLE_UNSET)
    {
        SetTimestep(tSamp);
    }
    return Solve(tStart, tEnd, mMaxDt, tSamp);
}


void AbstractCvodeCell::ComputeExceptVoltage(double tStart, double tEnd)
{
    double saved_voltage = GetVoltage();

    AbstractCardiacCellInterface::SetVoltageDerivativeToZero(true);
    SolveAndUpdateState(tStart, tEnd);
    AbstractCardiacCellInterface::SetVoltageDerivativeToZero(false);

    SetVoltage(saved_voltage); // In case of naughty models

#ifndef NDEBUG
    // Note that tests which rely on this throwing  (e.g. such-and-such a variable is out of range)
    // ought to be annotated with the NDEBUG macro
    VerifyStateVariables();
#endif // NDEBUG
}


void AbstractCvodeCell::SetVoltageDerivativeToZero(bool clamp)
{
    if (clamp != mSetVoltageDerivativeToZero)
    {
        ResetSolver();
    }
    AbstractCardiacCellInterface::SetVoltageDerivativeToZero(clamp);
}

unsigned AbstractCvodeCell::GetNumberOfStateVariables() const
{
    return AbstractCvodeSystem::GetNumberOfStateVariables();
}

unsigned AbstractCvodeCell::GetNumberOfParameters() const
{
    return AbstractCvodeSystem::GetNumberOfParameters();
}

std::vector<double> AbstractCvodeCell::GetStdVecStateVariables()
{
    std::vector<double> state_variables(GetNumberOfStateVariables());
    CopyToStdVector(AbstractCvodeSystem::rGetStateVariables(), state_variables);
    return state_variables;
}

const std::vector<std::string>& AbstractCvodeCell::rGetStateVariableNames() const
{
    return AbstractCvodeSystem::rGetStateVariableNames();
}

void AbstractCvodeCell::SetStateVariables(const std::vector<double>& rVariables)
{
    N_Vector vars = MakeNVector(rVariables);
    AbstractCvodeSystem::SetStateVariables(vars);
    DeleteVector(vars);
}

void AbstractCvodeCell::SetStateVariables(const N_Vector& rVariables)
{
    AbstractCvodeSystem::SetStateVariables(rVariables);
}

void AbstractCvodeCell::SetStateVariable(unsigned index, double newValue)
{
    AbstractCvodeSystem::SetStateVariable(index, newValue);
}

void AbstractCvodeCell::SetStateVariable(const std::string& rName, double newValue)
{
    AbstractCvodeSystem::SetStateVariable(rName, newValue);
}

double AbstractCvodeCell::GetAnyVariable(const std::string& rName, double time)
{
    return AbstractCvodeSystem::GetAnyVariable(rName,time);
}

double AbstractCvodeCell::GetParameter(const std::string& rParameterName)
{
    return AbstractCvodeSystem::GetParameter(rParameterName);
}

double AbstractCvodeCell::GetParameter(unsigned parameterIndex)
{
    return AbstractCvodeSystem::GetParameter(parameterIndex);
}

void AbstractCvodeCell::SetParameter(const std::string& rParameterName, double value)
{
    AbstractCvodeSystem::SetParameter(rParameterName,value);
}

#endif // CHASTE_CVODE
