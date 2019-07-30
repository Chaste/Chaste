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
#include "AbstractVanLeeuwen2009WntSwatCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

AbstractVanLeeuwen2009WntSwatCellCycleModel::AbstractVanLeeuwen2009WntSwatCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
   : AbstractWntOdeBasedCellCycleModel(pOdeSolver)
{
}

AbstractVanLeeuwen2009WntSwatCellCycleModel::AbstractVanLeeuwen2009WntSwatCellCycleModel(const AbstractVanLeeuwen2009WntSwatCellCycleModel& rModel)
   : AbstractWntOdeBasedCellCycleModel(rModel)
{
    /*
     * The member variables mDivideTime and mG2PhaseStartTime are
     * initialized in the AbstractOdeBasedPhaseBasedCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

void AbstractVanLeeuwen2009WntSwatCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);
    double beta_catenin_level =   mpOdeSystem->rGetStateVariables()[16]
                                + mpOdeSystem->rGetStateVariables()[17]
                                + mpOdeSystem->rGetStateVariables()[18]
                                + mpOdeSystem->rGetStateVariables()[19];

    // For mitogenic stimulus of 1/25.0 in Wnt equations
    if (beta_catenin_level < 10.188)
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

void AbstractVanLeeuwen2009WntSwatCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    assert(mpCell != nullptr);

    double wnt_level = GetWntLevel();

    InitialiseOdeSystem(wnt_level, mpCell->GetMutationState());

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}

void AbstractVanLeeuwen2009WntSwatCellCycleModel::AdjustOdeParameters(double currentTime)
{
    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[21] = GetWntLevel();

    // Use the cell's current mutation status as another input
    static_cast<VanLeeuwen2009WntSwatCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());
}

double AbstractVanLeeuwen2009WntSwatCellCycleModel::GetMembraneBoundBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[13] + mpOdeSystem->rGetStateVariables()[14];
}

double AbstractVanLeeuwen2009WntSwatCellCycleModel::GetCytoplasmicBetaCateninLevel()
{
    return  mpOdeSystem->rGetStateVariables()[7] + mpOdeSystem->rGetStateVariables()[8]
          + mpOdeSystem->rGetStateVariables()[9] + mpOdeSystem->rGetStateVariables()[10]
          + mpOdeSystem->rGetStateVariables()[11];
}

double AbstractVanLeeuwen2009WntSwatCellCycleModel::GetNuclearBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[16] +
           mpOdeSystem->rGetStateVariables()[17] +
           mpOdeSystem->rGetStateVariables()[18] +
           mpOdeSystem->rGetStateVariables()[19];
}

void AbstractVanLeeuwen2009WntSwatCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractWntOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
