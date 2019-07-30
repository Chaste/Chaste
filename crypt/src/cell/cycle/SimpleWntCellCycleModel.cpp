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

#include "SimpleWntCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

SimpleWntCellCycleModel::SimpleWntCellCycleModel()
    : mUseCellProliferativeTypeDependentG1Duration(false),
      mWntStemThreshold(0.8),
      mWntTransitThreshold(0.65),
      mWntLabelledThreshold(0.65)
{
}

SimpleWntCellCycleModel::SimpleWntCellCycleModel(const SimpleWntCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModel(rModel),
     mUseCellProliferativeTypeDependentG1Duration(rModel.mUseCellProliferativeTypeDependentG1Duration),
     mWntStemThreshold(rModel.mWntStemThreshold),
     mWntTransitThreshold(rModel.mWntTransitThreshold),
     mWntLabelledThreshold(rModel.mWntLabelledThreshold)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration and the cell's proliferative type are
     * (re)set as soon as InitialiseDaughterCell() is called on the
     * new cell-cycle model.
     */
}

AbstractCellCycleModel* SimpleWntCellCycleModel::CreateCellCycleModel()
{
    return new SimpleWntCellCycleModel(*this);
}

bool SimpleWntCellCycleModel::GetUseCellProliferativeTypeDependentG1Duration() const
{
    return mUseCellProliferativeTypeDependentG1Duration;
}

void SimpleWntCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}

void SimpleWntCellCycleModel::SetG1Duration()
{
    assert(mpCell != nullptr);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        if (mUseCellProliferativeTypeDependentG1Duration)
        {
            mG1Duration = p_gen->NormalRandomDeviate(GetStemCellG1Duration(), 1.0);
        }
        else
        {
            // Normally stem cells should behave just like transit cells in a Wnt simulation
            mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 1.0);
        }
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 1.0);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}

double SimpleWntCellCycleModel::GetWntLevel() const
{
    assert(mpCell != nullptr);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }
    return level;
}

WntConcentrationType SimpleWntCellCycleModel::GetWntType()
{
    WntConcentrationType wnt_type;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case UNSIGNED_UNSET:
        {
            // If you trip this you have tried to use a simulation without setting the dimension.
            NEVER_REACHED;
        }
        default:
            NEVER_REACHED;
    }
    return wnt_type;
}

void SimpleWntCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold
    double wnt_division_threshold = DBL_MAX;

    // Set up under what level of Wnt stimulus a cell will divide
    if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        wnt_division_threshold = mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
    {
        // should be less than healthy values
        wnt_division_threshold = 0.77*mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
    {
        // less than above value
        wnt_division_threshold = 0.155*mWntTransitThreshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;
    }
    else
    {
        NEVER_REACHED;
    }

    if (mpCell->HasCellProperty<CellLabel>())
    {
        wnt_division_threshold = mWntLabelledThreshold;
    }

    double wnt_level = GetWntLevel();
    WntConcentrationType wnt_type = GetWntType();

    // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt_division_threshold
    if (wnt_level >= wnt_division_threshold)
    {
        // For a RADIAL Wnt type, override the cell type to StemCellProliferativeType if the Wnt stimulus exceeds a higher threshold
        if ((wnt_type == RADIAL) && (wnt_level > mWntStemThreshold))
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
            boost::shared_ptr<AbstractCellProperty> p_stem_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_stem_type);
        }
        else
        {
            boost::shared_ptr<AbstractCellProperty> p_transit_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_transit_type);
        }
    }
    else
    {
        // The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();
}

void SimpleWntCellCycleModel::InitialiseDaughterCell()
{
    WntConcentrationType wnt_type = GetWntType();

    if (wnt_type == RADIAL)
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }

    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}

bool SimpleWntCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double SimpleWntCellCycleModel::GetWntStemThreshold() const
{
    return mWntStemThreshold;
}

void SimpleWntCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double SimpleWntCellCycleModel::GetWntTransitThreshold() const
{
    return mWntTransitThreshold;
}

void SimpleWntCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
    //assert(wntTransitThreshold <= 1.0);
    //assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double SimpleWntCellCycleModel::GetWntLabelledThreshold() const
{
    return mWntLabelledThreshold;
}

void SimpleWntCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
//    assert(wntLabelledThreshold <= 1.0);
//    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}

void SimpleWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<UseCellProliferativeTypeDependentG1Duration>" << mUseCellProliferativeTypeDependentG1Duration << "</UseCellProliferativeTypeDependentG1Duration>\n";
    *rParamsFile << "\t\t\t<WntStemThreshold>" << mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile << "\t\t\t<WntTransitThreshold>" << mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile << "\t\t\t<WntLabelledThreshold>" << mWntLabelledThreshold << "</WntLabelledThreshold>\n";

    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleWntCellCycleModel)
