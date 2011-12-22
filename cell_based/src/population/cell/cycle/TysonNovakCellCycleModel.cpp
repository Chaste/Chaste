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

#include "TysonNovakCellCycleModel.hpp"

TysonNovakCellCycleModel::TysonNovakCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeBasedCellCycleModel(SimulationTime::Instance()->GetTime(), pOdeSolver)
{
    if (!mpOdeSolver)
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<TysonNovakCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
        mpOdeSolver->SetTolerances(1e-6, 1e-8);
#else
        mpOdeSolver = CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver>::Instance();
        mpOdeSolver->SetSizeOfOdeSystem(6);
        mpOdeSolver->Initialise();
        SetDt(0.1/60.0);
#endif //CHASTE_CVODE
    }
}

void TysonNovakCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    mpOdeSystem = new TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    AbstractCellCycleModel::Initialise();
}

void TysonNovakCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedCellCycleModel::ResetForDivision();

    assert(mpOdeSystem != NULL);

    /**
     * This model needs the protein concentrations and phase resetting to G0/G1.
     *
     * In theory, the solution to the Tyson-Novak equations should exhibit stable
     * oscillations, and we only need to halve the mass of the cell each period.
     *
     * However, the backward Euler solver used to solve the equations
     * currently returns a solution that diverges after long times, so
     * we must reset the initial conditions each period.
     *
     * When running with CVODE however we can use the halving the mass of the cell method.
     */
#ifdef CHASTE_CVODE
    mpOdeSystem->rGetStateVariables()[5] = 0.5*mpOdeSystem->rGetStateVariables()[5];
#else
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
#endif //CHASTE_CVODE
}

void TysonNovakCellCycleModel::InitialiseDaughterCell()
{
    if (mCellProliferativeType == STEM)
    {
        mCellProliferativeType = TRANSIT;
    }
}

AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    TysonNovakCellCycleModel* p_model = new TysonNovakCellCycleModel(mpOdeSolver);

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide, mDt, mpOdeSolver)
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: the member variable mDimension remains unset, since this cell-cycle
     * model does not need to know the spatial dimension, so if we were to call
     * SetDimension() on the new cell-cycle model an exception would be triggered;
     * hence we do not set this member variable.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetDivideTime(mDivideTime);
    p_model->SetFinishedRunningOdes(mFinishedRunningOdes);
    p_model->SetG2PhaseStartTime(mG2PhaseStartTime);
    p_model->SetLastTime(mLastTime);

    /*
     * Create the new cell-cycle model's ODE system and use the current values
     * of the state variables in mpOdeSystem as an initial condition.
     */
    assert(mpOdeSystem);
    p_model->SetOdeSystem(new TysonNovak2001OdeSystem);
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());

    return p_model;
}

double TysonNovakCellCycleModel::GetSDuration()
{
    /**
     * Tyson & Novak pretends it is running ODEs in just G1,
     * but they really represent the whole cell cycle, so
     * we set the other phases to zero.
     */
    return 0.0;
}

double TysonNovakCellCycleModel::GetG2Duration()
{
    /**
     * Tyson & Novak pretends it is running ODEs in just G1,
     * but they really represent the whole cell cycle so
     * we set the other phases to zero.
     */
    return 0.0;
}

double TysonNovakCellCycleModel::GetMDuration()
{
    /**
     * Tyson & Novak pretends it is running ODEs in just G1,
     * but they really represent the whole cell cycle so
     * we set the other phases to zero.
     */
    return 0.0;
}

double TysonNovakCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 1.25;
}

double TysonNovakCellCycleModel::GetAverageStemCellCycleTime()
{
    return 1.25;
}


bool TysonNovakCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

void TysonNovakCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovakCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovakCellCycleModel)
