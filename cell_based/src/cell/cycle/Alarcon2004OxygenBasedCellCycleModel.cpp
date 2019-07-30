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

#include "Alarcon2004OxygenBasedCellCycleModel.hpp"

#include "CellLabel.hpp"

Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeBasedPhaseBasedCellCycleModel(SimulationTime::Instance()->GetTime(), pOdeSolver)
{
    if (!mpOdeSolver)
    {
        mpOdeSolver = CellCycleModelOdeSolver<Alarcon2004OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
    }
    SetDt(0.0001);
}

Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel(const Alarcon2004OxygenBasedCellCycleModel& rModel)
    : AbstractOdeBasedPhaseBasedCellCycleModel(rModel)
{
    /*
     * Initialize only those member variables defined in this class.
     * Create the new cell-cycle model's ODE system and use the current
     * values of the state variables in mpOdeSystem as an initial condition.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */

    assert(rModel.GetOdeSystem());
    bool is_labelled = rModel.mpCell->HasCellProperty<CellLabel>();
    SetOdeSystem(new Alarcon2004OxygenBasedCellCycleOdeSystem(rModel.mpCell->GetCellData()->GetItem("oxygen"), is_labelled));
    SetStateVariables(rModel.GetOdeSystem()->rGetStateVariables());
}

AbstractCellCycleModel* Alarcon2004OxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new Alarcon2004OxygenBasedCellCycleModel(*this);
}

void Alarcon2004OxygenBasedCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedPhaseBasedCellCycleModel::ResetForDivision();
    assert(mpOdeSystem != nullptr);

    // This model needs the protein concentrations and phase resetting to G0/G1.
    // Keep the oxygen concentration the same but reset everything else
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i=0; i<5; i++)
    {
        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
}

void Alarcon2004OxygenBasedCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    assert(mpCell != nullptr);

    bool is_labelled = mpCell->HasCellProperty<CellLabel>();
    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(mpCell->GetCellData()->GetItem("oxygen"), is_labelled);
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    AbstractOdeBasedPhaseBasedCellCycleModel::Initialise();
}

void Alarcon2004OxygenBasedCellCycleModel::AdjustOdeParameters(double currentTime)
{
    // Pass this time step's oxygen concentration into the solver as a constant over this time step
    mpOdeSystem->rGetStateVariables()[5] = mpCell->GetCellData()->GetItem("oxygen");

    // Use whether the cell is currently labelled as another input
    bool is_labelled = mpCell->HasCellProperty<CellLabel>();
    static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->SetIsLabelled(is_labelled);
}

void Alarcon2004OxygenBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeBasedPhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Alarcon2004OxygenBasedCellCycleModel)
