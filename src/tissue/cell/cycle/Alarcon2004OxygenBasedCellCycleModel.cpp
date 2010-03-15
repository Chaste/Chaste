/*

Copyright (C) University of Oxford, 2005-2010

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
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"


RungeKutta4IvpOdeSolver Alarcon2004OxygenBasedCellCycleModel::msSolver;

Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel(const Alarcon2004OxygenBasedCellCycleModel& rOtherModel)
    : AbstractOdeBasedCellCycleModelWithStoppingEvent(rOtherModel),
      mDimension(rOtherModel.mDimension)
{
    if (rOtherModel.mpOdeSystem != NULL)
    {
        mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(*static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(rOtherModel.mpOdeSystem));
    }
}


Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel()
    : AbstractOdeBasedCellCycleModelWithStoppingEvent()
{
}


Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                                                                           const CryptCellMutationState& rMutationState,
                                                                           const unsigned& rDimension)
{
    mDimension = rDimension;
    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(rParentProteinConcentrations[5], rMutationState);

    // Set the model to be the same as the parent cell.
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
}


void Alarcon2004OxygenBasedCellCycleModel::ResetForDivision()
{
	AbstractOdeBasedCellCycleModelWithStoppingEvent::ResetForDivision();
    assert(mpOdeSystem!=NULL);

    // This model needs the protein concentrations and phase resetting to G0/G1.
    // Keep the oxygen concentration the same but reset everything else
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i=0; i<5; i++)
    {
        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
}

AbstractCellCycleModel* Alarcon2004OxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new Alarcon2004OxygenBasedCellCycleModel(*this);
}

void Alarcon2004OxygenBasedCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(*mpCell,0), mpCell->GetMutationState());
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(*mpCell,0), mpCell->GetMutationState());
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<DIM>::Instance()->GetValue(*mpCell,0), mpCell->GetMutationState());
            break;
        }
        default:
            NEVER_REACHED;
    }

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}


bool Alarcon2004OxygenBasedCellCycleModel::SolveOdeToTime(double currentTime)
{
    double dt = 0.0001; // the time step must be this small because the Runge-Kutter solver has poor stability

    // Pass this time step's oxygen concentration into the solver as a constant over this timestep
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            mpOdeSystem->rGetStateVariables()[5] = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            mpOdeSystem->rGetStateVariables()[5] = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            mpOdeSystem->rGetStateVariables()[5] = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
            break;
        }
        default:
            NEVER_REACHED;
    }

    // Use the cell's current mutation status as another input
    static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

    msSolver.SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, dt);
    return msSolver.StoppingEventOccurred();
}


double Alarcon2004OxygenBasedCellCycleModel::GetOdeStopTime()
{
    assert(msSolver.StoppingEventOccurred());
    return msSolver.GetStoppingTime();
}

unsigned Alarcon2004OxygenBasedCellCycleModel::GetDimension()
{
    return mDimension;
}



// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)
