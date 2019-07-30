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

#include "UblasIncludes.hpp"
#include "SingleOdeWntCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

SingleOdeWntCellCycleModel::SingleOdeWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : CellCycleModelOdeHandler(DOUBLE_UNSET, pOdeSolver),
      mBetaCateninDivisionThreshold(100)// MAGIC NUMBER!
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<SingleOdeWntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

SingleOdeWntCellCycleModel::SingleOdeWntCellCycleModel(const SingleOdeWntCellCycleModel& rModel)
    : SimpleWntCellCycleModel(rModel),
      CellCycleModelOdeHandler(rModel),
      mBetaCateninDivisionThreshold(rModel.mBetaCateninDivisionThreshold)
{
    /*
     * Initialize only those member variables defined in this class.
     * Create the new cell-cycle model's ODE system and use the current values
     * of the state variables in mpOdeSystem as an initial condition.
     *
     * The member variables mUseCellProliferativeTypeDependentG1Duration,
     * mWntStemThreshold, mWntTransitThreshold and mWntLabelledThreshold
     * are initialized in the SimpleWntCellCycleModel constructor.
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
    double wnt_level = rModel.GetWntLevel();
    SetOdeSystem(new Mirams2010WntOdeSystem(wnt_level, rModel.mpCell->GetMutationState()));
    SetStateVariables(rModel.GetOdeSystem()->rGetStateVariables());
}

AbstractCellCycleModel* SingleOdeWntCellCycleModel::CreateCellCycleModel()
{
    return new SingleOdeWntCellCycleModel(*this);
}

void SingleOdeWntCellCycleModel::UpdateCellCyclePhase()
{
    assert(SimulationTime::Instance()->IsStartTimeSetUp());
    SolveOdeToTime(SimulationTime::Instance()->GetTime());
    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
    AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase(); /// Don't call the SimpleWntCellCycleModel - it will overwrite the beta catenin levels.
}

void SingleOdeWntCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    assert(mpCell != nullptr);

    double wnt_level = this->GetWntLevel();
    mpOdeSystem = new Mirams2010WntOdeSystem(wnt_level, mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    // MAGIC NUMBER!
    //mBetaCateninDivisionThreshold = 100.0;

    // This call actually sets up the G1 phase to something sensible (random number generated)
    SimpleWntCellCycleModel::Initialise();

    SetLastTime(mBirthTime);

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}

void SingleOdeWntCellCycleModel::AdjustOdeParameters(double currentTime)
{
    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[2] = this->GetWntLevel();

    // Use the cell's current mutation status as another input
    static_cast<Mirams2010WntOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());
}

void SingleOdeWntCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    if (GetBetaCateninConcentration() < GetBetaCateninDivisionThreshold())
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
         * would be creating a new CellPropertyRegistry. In this case the cell proliferative
         * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
         * would be incorrect. We must therefore access the CellProliferativeType via the cell's
         * CellPropertyCollection.
         */
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    else
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }
}

double SingleOdeWntCellCycleModel::GetBetaCateninConcentration()
{
    return mpOdeSystem->rGetStateVariables()[0] + mpOdeSystem->rGetStateVariables()[1];
}

void SingleOdeWntCellCycleModel::SetBetaCateninDivisionThreshold(double betaCateninDivisionThreshold)
{
    mBetaCateninDivisionThreshold = betaCateninDivisionThreshold;
}

double SingleOdeWntCellCycleModel::GetBetaCateninDivisionThreshold()
{
    return mBetaCateninDivisionThreshold;
}

void SingleOdeWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    SimpleWntCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SingleOdeWntCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(SingleOdeWntCellCycleModel)
