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

#include "UblasIncludes.hpp"
#include "DeltaNotchCellCycleModel.hpp"
#include "CellwiseData.hpp"
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
    p_model->SetOdeSystem(new DeltaNotchOdeSystem(mean_neighbouring_delta));

    // Use the current values of the state variables in mpOdeSystem as an initial condition for the new cell-cycle model's ODE system
    assert(mpOdeSystem);
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());

    // Set the values of the new cell-cycle model's member variables
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetBirthTime(mBirthTime);
    p_model->SetLastTime(mLastTime);
    p_model->SetDimension(mDimension);
    p_model->SetGeneration(mGeneration);
    p_model->SetMaxTransitGenerations(mMaxTransitGenerations);
    p_model->SetCellProliferativeType(mCellProliferativeType);

    return p_model;
}

void DeltaNotchCellCycleModel::UpdateCellCyclePhase()
{
    assert(SimulationTime::Instance()->IsStartTimeSetUp());
    UpdateDeltaNotch();
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void DeltaNotchCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    mpOdeSystem = new DeltaNotchOdeSystem;
    if(mInitialConditions == std::vector<double>())
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
    assert(initialConditions.size() == 3);
    mInitialConditions = initialConditions;
}

void DeltaNotchCellCycleModel::UpdateDeltaNotch()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    double mean_delta;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            mean_delta = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            mean_delta = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            mean_delta = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        default:
            NEVER_REACHED;
    }

    mpOdeSystem->rGetStateVariables()[2] = mean_delta;

    SolveOdeToTime(SimulationTime::Instance()->GetTime());
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
    double mean_neighbouring_delta = mpOdeSystem->rGetStateVariables()[2];
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
