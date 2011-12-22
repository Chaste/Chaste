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

#include "ContactInhibitionCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "CellwiseData.hpp"

ContactInhibitionCellCycleModel::ContactInhibitionCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mQuiescentVolumeFraction(DOUBLE_UNSET),
      mEquilibriumVolume(DOUBLE_UNSET),
      mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
      mCurrentQuiescentDuration(0.0)
{
}

void ContactInhibitionCellCycleModel::UpdateCellCyclePhase()
{
    if ((mQuiescentVolumeFraction == DOUBLE_UNSET) || (mEquilibriumVolume == DOUBLE_UNSET))
    {
        EXCEPTION("The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");
    }

    // Get cell volume
    double cell_volume;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            cell_volume = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            cell_volume = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            cell_volume = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        default:
            NEVER_REACHED;
    }

    // Removes the cell label
    mpCell->RemoveCellProperty<CellLabel>();

    if (mCurrentCellCyclePhase == G_ONE_PHASE)
    {
        // Update G1 duration based on cell volume
        double dt = SimulationTime::Instance()->GetTimeStep();
        double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

        if (cell_volume < quiescent_volume)
        {
            // Update the duration of the current period of hypoxia
            mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;
            mG1Duration += dt;
            mpCell->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellLabel>());
        }
        else
        {
            // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
            mCurrentQuiescentDuration = 0.0;
            mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
        }
    }

    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    if (mpCell->GetCellCycleModel()->GetCellProliferativeType()==DIFFERENTIATED)
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if ( time_since_birth < GetMDuration() )
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
}

AbstractCellCycleModel* ContactInhibitionCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide, mTimeSpentInG1Phase,
     * mCurrentHypoxicDuration, mCurrentHypoxiaOnsetTime) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetQuiescentVolumeFraction(mQuiescentVolumeFraction);
    p_model->SetEquilibriumVolume(mEquilibriumVolume);
    p_model->SetCurrentQuiescentOnsetTime(mCurrentQuiescentOnsetTime);
    p_model->SetCurrentQuiescentDuration(mCurrentQuiescentDuration);

    return p_model;
}

void ContactInhibitionCellCycleModel::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
}

double ContactInhibitionCellCycleModel::GetQuiescentVolumeFraction()
{
    return mQuiescentVolumeFraction;
}

void ContactInhibitionCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double ContactInhibitionCellCycleModel::GetEquilibriumVolume()
{
    return mEquilibriumVolume;
}

void ContactInhibitionCellCycleModel::SetCurrentQuiescentDuration(double currentQuiescentDuration)
{
    mCurrentQuiescentDuration = currentQuiescentDuration;
}

double ContactInhibitionCellCycleModel::GetCurrentQuiescentDuration()
{
    return mCurrentQuiescentDuration;
}

void ContactInhibitionCellCycleModel::SetCurrentQuiescentOnsetTime(double currentQuiescentOnsetTime)
{
    mCurrentQuiescentOnsetTime = currentQuiescentOnsetTime;
}

double ContactInhibitionCellCycleModel::GetCurrentQuiescentOnsetTime()
{
    return mCurrentQuiescentOnsetTime;
}

void ContactInhibitionCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<QuiescentVolumeFraction>" << mQuiescentVolumeFraction << "</QuiescentVolumeFraction>\n";
    *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ContactInhibitionCellCycleModel)
