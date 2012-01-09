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
    //Note that tests which rely on this throwing  (e.g. such-and-such a variable is out of range)
    //ought to be anotated with the NDEBUG macro
    VerifyStateVariables();
#endif // NDEBUG
}

void AbstractCardiacCell::SetVoltage(double voltage)
{
    SetAnyVariable(mVoltageIndex, voltage);
}

double AbstractCardiacCell::GetVoltage()
{
    return GetAnyVariable(mVoltageIndex);
}

double AbstractCardiacCell::GetIntracellularCalciumConcentration()
{
    EXCEPTION("AbstractCardiacCell::GetIntracellularCalciumConcentration() called. Either model has no [Ca_i] or method has not been implemented yet");
}



#include "LuoRudy1991.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
void AbstractCardiacCell::CheckForArchiveFix()
{
    if (dynamic_cast<CellLuoRudy1991FromCellML*>(this) || dynamic_cast<CellLuoRudy1991FromCellMLBackwardEuler*>(this))
    {
        // The LR91 model saved in previous Chaste versions had a different ordering of state variables...
        // Old is h, j, m, CaI, V, d, f, x
        // New is V, m, h, j, d, f, X, [Ca]
        assert(GetNumberOfStateVariables() == 8);
        unsigned var_index_map[8] = {2, 3, 1, 7, 0, 4, 5, 6};
        std::vector<double> old_state(this->mStateVariables);
        for (unsigned i=0; i<8; i++)
        {
            this->mStateVariables[var_index_map[i]] = old_state[i];
        }
        // It also didn't use to have parameters...
        this->mParameters.resize(this->rGetParameterNames().size());
        assert(this->mParameters.size() == 2u);
        this->mParameters[0] = 23.0;
        this->mParameters[1] = 0.282;
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
