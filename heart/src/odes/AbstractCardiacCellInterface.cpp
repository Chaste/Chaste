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

#include "HeartConfig.hpp" // First for Boost 1.33/PETSc 2.2
#include "AbstractCardiacCellInterface.hpp"
#include "Exception.hpp"

#include "Citations.hpp"
static PetscBool CardiacCellChasteCite = PETSC_FALSE;
const char CardiacCellChasteCitation[] = "@article{cooper2015cce,\n"
" title = {{Cellular cardiac electrophysiology modelling with Chaste and CellML}},\n"
"  author = {Cooper, Jonathan and Spiteri, Raymond and Mirams, Gary R},\n"
"  year = {2015},\n"
"  journal = {Frontiers in Physiology},\n"
"  pages = {511},\n"
"  url = {http://journal.frontiersin.org/Journal/10.3389/fphys.2014.00511},\n"
"  volume = {5},\n"
"  doi = {10.3389/fphys.2014.00511},\n"
"}\n";

AbstractCardiacCellInterface::AbstractCardiacCellInterface(
            boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
            unsigned voltageIndex,
            boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : mVoltageIndex(voltageIndex),
      mpOdeSolver(pOdeSolver),
      mpIntracellularStimulus(pIntracellularStimulus),
      mSetVoltageDerivativeToZero(false),
      mIsUsedInTissue(false),
      mHasDefaultStimulusFromCellML(false),
      mFixedVoltage(DOUBLE_UNSET)
{
    // Record a reference for the calculations performed using this class,
    // can be extracted with the '-citations' flag as an argument to any executable.
    Citations::Register(CardiacCellChasteCitation, &CardiacCellChasteCite);
}


AbstractCardiacCellInterface::~AbstractCardiacCellInterface()
{
}


unsigned AbstractCardiacCellInterface::GetVoltageIndex()
{
    return mVoltageIndex;
}


void AbstractCardiacCellInterface::SetStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus)
{
    SetIntracellularStimulusFunction(pStimulus);
}


double AbstractCardiacCellInterface::GetStimulus(double time)
{
    return GetIntracellularStimulus(time);
}


void AbstractCardiacCellInterface::SetIntracellularStimulusFunction(boost::shared_ptr<AbstractStimulusFunction> pStimulus)
{
    mpIntracellularStimulus = pStimulus;
}


double AbstractCardiacCellInterface::GetIntracellularStimulus(double time)
{
    return mpIntracellularStimulus->GetStimulus(time);
}


double AbstractCardiacCellInterface::GetIntracellularAreaStimulus(double time)
{
    double stim;
    if (mIsUsedInTissue)
    {
        // Convert from uA/cm^3 to uA/cm^2 by dividing by Am
        stim = GetIntracellularStimulus(time) / HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    }
    else
    {
        stim = GetIntracellularStimulus(time);
    }
    return stim;
}

void AbstractCardiacCellInterface::SetUsedInTissueSimulation(bool tissue)
{
    mIsUsedInTissue = tissue;
}

boost::shared_ptr<RegularStimulus> AbstractCardiacCellInterface::UseCellMLDefaultStimulus()
{
    assert(!mHasDefaultStimulusFromCellML);
    EXCEPTION("This class has no default stimulus from CellML metadata.");
    return boost::shared_ptr<RegularStimulus>();
}

bool AbstractCardiacCellInterface::HasCellMLDefaultStimulus()
{
    return mHasDefaultStimulusFromCellML;
}

boost::shared_ptr<AbstractStimulusFunction> AbstractCardiacCellInterface::GetStimulusFunction()
{
    return mpIntracellularStimulus;
}

// Methods needed by boost serialization.
const boost::shared_ptr<AbstractStimulusFunction> AbstractCardiacCellInterface::GetStimulusFunction() const
{
    return mpIntracellularStimulus;
}

const boost::shared_ptr<AbstractIvpOdeSolver> AbstractCardiacCellInterface::GetSolver() const
{
    return mpOdeSolver;
}

void AbstractCardiacCellInterface::SetSolver(boost::shared_ptr<AbstractIvpOdeSolver> pSolver)
{
    mpOdeSolver = pSolver;
}

void AbstractCardiacCellInterface::SetVoltageDerivativeToZero(bool clamp)
{
    mSetVoltageDerivativeToZero = clamp;
    if (clamp)
    {
        mFixedVoltage = GetVoltage();
    }
}

void AbstractCardiacCellInterface::SetFixedVoltage(double voltage)
{
    mFixedVoltage = voltage;
}

double AbstractCardiacCellInterface::GetIntracellularCalciumConcentration()
{
    EXCEPTION("AbstractCardiacCellInterface::GetIntracellularCalciumConcentration() called. "
              "Either model has no [Ca_i] or method has not been implemented yet");
}


