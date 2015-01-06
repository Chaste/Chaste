/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "UblasIncludes.hpp"
#include "DeltaNotchCellCycleModel.hpp"
#include "CellCycleModelOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "Exception.hpp"



DeltaNotchCellCycleModel::DeltaNotchCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : CellCycleModelOdeHandler(DOUBLE_UNSET, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

AbstractCellCycleModel* DeltaNotchCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel(this->mpOdeSolver);

    // Create the new cell-cycle model's ODE system
    double mean_neighbouring_delta = GetMeanNeighbouringDelta();
    p_model->SetOdeSystem(new DeltaNotchOdeSystem);
    p_model->GetOdeSystem()->SetParameter("Mean Delta", mean_neighbouring_delta);

    // Use the current values of the state variables in mpOdeSystem as an initial condition for the new cell-cycle model's ODE system
    assert(mpOdeSystem);
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());

    // Set the values of the new cell-cycle model's member variables
    p_model->SetBirthTime(mBirthTime);
    p_model->SetLastTime(mLastTime);
    p_model->SetDimension(mDimension);
    p_model->SetGeneration(mGeneration);
    p_model->SetMaxTransitGenerations(mMaxTransitGenerations);

    return p_model;
}

void DeltaNotchCellCycleModel::UpdateCellCyclePhase()
{
    assert(SimulationTime::Instance()->IsStartTimeSetUp());
    UpdateDeltaNotch();
    SolveOdeToTime(SimulationTime::Instance()->GetTime());
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void DeltaNotchCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    mpOdeSystem = new DeltaNotchOdeSystem;
    if (mInitialConditions == std::vector<double>())
    {
        mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    }
    else
    {
        mpOdeSystem->SetStateVariables(mInitialConditions);
    }

    StochasticDurationGenerationBasedCellCycleModel::Initialise();

    SetLastTime(mBirthTime);
}


void DeltaNotchCellCycleModel::SetInitialConditions(std::vector<double> initialConditions)
{
    assert(initialConditions.size() == 2);
    mInitialConditions = initialConditions;
}

void DeltaNotchCellCycleModel::UpdateDeltaNotch()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    double mean_delta = mpCell->GetCellData()->GetItem("mean delta");
    mpOdeSystem->SetParameter("Mean Delta", mean_delta);
}

double DeltaNotchCellCycleModel::GetNotch()
{
    assert(mpOdeSystem != NULL);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

double DeltaNotchCellCycleModel::GetDelta()
{
    assert(mpOdeSystem != NULL);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

double DeltaNotchCellCycleModel::GetMeanNeighbouringDelta()
{
    assert(mpOdeSystem != NULL);
    double mean_neighbouring_delta = mpOdeSystem->GetParameter("Mean Delta");
    return mean_neighbouring_delta;
}

void DeltaNotchCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output.

    // Call direct parent class
    StochasticDurationGenerationBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchCellCycleModel)
