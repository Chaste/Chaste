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

#include "TysonNovakCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

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
        SetDt(0.1/90.0);
#endif //CHASTE_CVODE
    }
}

TysonNovakCellCycleModel::TysonNovakCellCycleModel(const TysonNovakCellCycleModel& rModel)
    : AbstractOdeBasedCellCycleModel(rModel)
{
    /*
     * Initialize only those member variables defined in this class.
     * Create the new cell-cycle model's ODE system and use the current
     * values of the state variables in mpOdeSystem as an initial condition.
     *
     * The member variable mDivideTime is initialized in the
     * AbstractOdeBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that the cell proliferative type is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
    assert(rModel.GetOdeSystem());
    SetOdeSystem(new TysonNovak2001OdeSystem);
    SetStateVariables(rModel.GetOdeSystem()->rGetStateVariables());
}

void TysonNovakCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    mpOdeSystem = new TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    AbstractOdeBasedCellCycleModel::Initialise();
}

void TysonNovakCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedCellCycleModel::ResetForDivision();

    assert(mpOdeSystem != nullptr);

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
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
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
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }
}

AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
{
    return new TysonNovakCellCycleModel(*this);
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
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovakCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovakCellCycleModel)
