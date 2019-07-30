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
#include "AbstractCardiacCell.hpp"

#include <cassert>
#include <iostream>

#include "HeartConfig.hpp"
#include "Exception.hpp"

AbstractCardiacCell::AbstractCardiacCell(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                         unsigned numberOfStateVariables,
                                         unsigned voltageIndex,
                                         boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCellInterface(pOdeSolver, voltageIndex, pIntracellularStimulus),
      AbstractOdeSystem(numberOfStateVariables),
      mDt(HeartConfig::Instance()->GetOdeTimeStep())
{
    // The second clause is to allow for FakeBathCell.
    assert(voltageIndex < mNumberOfStateVariables || mNumberOfStateVariables == 0);
}

AbstractCardiacCell::~AbstractCardiacCell()
{
}

void AbstractCardiacCell::Init()
{
    ResetToInitialConditions();
    mParameters.resize(rGetParameterNames().size());
}

void AbstractCardiacCell::SetTimestep(double dt)
{
    mDt = dt;
}

void AbstractCardiacCell::SolveAndUpdateState(double tStart, double tEnd)
{
    mpOdeSolver->SolveAndUpdateStateVariable(this, tStart, tEnd, mDt);
}

OdeSolution AbstractCardiacCell::Compute(double tStart, double tEnd, double tSamp)
{
    if (tSamp < mDt)
    {
        tSamp = mDt;
    }
    return mpOdeSolver->Solve(this, rGetStateVariables(), tStart, tEnd, mDt, tSamp);
}

void AbstractCardiacCell::ComputeExceptVoltage(double tStart, double tEnd)
{
    double saved_voltage = GetVoltage();

    SetVoltageDerivativeToZero(true);
    mpOdeSolver->SolveAndUpdateStateVariable(this, tStart, tEnd, mDt);
    SetVoltageDerivativeToZero(false);

    SetVoltage(saved_voltage); // In case of naughty models

#ifndef NDEBUG
    // Note that tests which rely on this throwing  (e.g. such-and-such a variable is out of range)
    // ought to be annotated with the NDEBUG macro
    VerifyStateVariables();
#endif // NDEBUG
}

void AbstractCardiacCell::SetVoltage(double voltage)
{
    SetAnyVariable(mVoltageIndex, voltage);
    SetFixedVoltage(voltage);
}

double AbstractCardiacCell::GetVoltage()
{
    return AbstractOdeSystem::GetAnyVariable(mVoltageIndex);
}

unsigned AbstractCardiacCell::GetNumberOfStateVariables() const
{
    return AbstractOdeSystem::GetNumberOfStateVariables();
}

unsigned AbstractCardiacCell::GetNumberOfParameters() const
{
    return AbstractOdeSystem::GetNumberOfParameters();
}

std::vector<double> AbstractCardiacCell::GetStdVecStateVariables()
{
    return AbstractOdeSystem::GetStateVariables();
}

const std::vector<std::string>& AbstractCardiacCell::rGetStateVariableNames() const
{
    return AbstractOdeSystem::rGetStateVariableNames();
}

void AbstractCardiacCell::SetStateVariables(const std::vector<double>& rVariables)
{
    AbstractOdeSystem::SetStateVariables(rVariables);
}

void AbstractCardiacCell::SetStateVariable(unsigned index, double newValue)
{
    AbstractOdeSystem::SetStateVariable(index, newValue);
}

void AbstractCardiacCell::SetStateVariable(const std::string& rName, double newValue)
{
    AbstractOdeSystem::SetStateVariable(rName, newValue);
}

double AbstractCardiacCell::GetAnyVariable(const std::string& rName, double time)
{
    return AbstractOdeSystem::GetAnyVariable(rName, time);
}

double AbstractCardiacCell::GetParameter(const std::string& rParameterName)
{
    return AbstractOdeSystem::GetParameter(rParameterName);
}

double AbstractCardiacCell::GetParameter(unsigned parameterIndex)
{
    return AbstractOdeSystem::GetParameter(parameterIndex);
}

void AbstractCardiacCell::SetParameter(const std::string& rParameterName, double value)
{
    AbstractOdeSystem::SetParameter(rParameterName,value);
}

#include "LuoRudy1991.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
void AbstractCardiacCell::CheckForArchiveFix()
{
    if (dynamic_cast<CellLuoRudy1991FromCellML*>(this) || dynamic_cast<CellLuoRudy1991FromCellMLBackwardEuler*>(this))
    {
        // If you do reach this, then you can try un-commenting the below, but you would be better re-generating your archive with a new version of Chaste.
        NEVER_REACHED;
//        // The LR91 model saved in previous Chaste versions had a different ordering of state variables...
//        // Old is h, j, m, CaI, V, d, f, x
//        // New is V, m, h, j, d, f, X, [Ca]
//        assert(GetNumberOfStateVariables() == 8);
//        unsigned var_index_map[8] = {2, 3, 1, 7, 0, 4, 5, 6};
//        std::vector<double> old_state(this->mStateVariables);
//        for (unsigned i=0; i<8; i++)
//        {
//            this->mStateVariables[var_index_map[i]] = old_state[i];
//        }
//        // It also didn't used to have parameters...
//        this->mParameters.resize(this->rGetParameterNames().size());
//        assert(this->mParameters.size() == 3u);
//        this->mParameters[0] = 0.09;
//        this->mParameters[1] = 23.0;
//        this->mParameters[2] = 0.282;
    }
}


/*
 *  METHODS NEEDED BY FAST CARDIAC CELLS
 */
void AbstractCardiacCell::SetState(CellModelState state)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

void AbstractCardiacCell::SetSlowValues(const std::vector<double> &rSlowValues)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

void AbstractCardiacCell::GetSlowValues(std::vector<double>& rSlowValues)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

bool AbstractCardiacCell::IsFastOnly()
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

unsigned AbstractCardiacCell::GetNumSlowValues()
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}

void AbstractCardiacCell::AdjustOutOfRangeSlowValues(std::vector<double>& rSlowValues)
{
    EXCEPTION("Non fast-slow cell model being used in a fast-slow problem.");
}
