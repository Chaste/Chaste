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

#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellPropertyRegistry.hpp"
#include "Exception.hpp"

SimpleOxygenBasedCellCycleModel::SimpleOxygenBasedCellCycleModel()
    : mTimeSpentInG1Phase(0.0),
      mCurrentHypoxicDuration(0.0),
      mHypoxicConcentration(0.4),
      mQuiescentConcentration(1.0),
      mCriticalHypoxicDuration(2.0)
{
    mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
}

double SimpleOxygenBasedCellCycleModel::GetCurrentHypoxicDuration()
{
    return mCurrentHypoxicDuration;
}

double SimpleOxygenBasedCellCycleModel::GetCurrentHypoxiaOnsetTime()
{
    return mCurrentHypoxiaOnsetTime;
}

void SimpleOxygenBasedCellCycleModel::UpdateCellCyclePhase()
{
    // mG1Duration is set when the cell-cycle model is given a cell

    bool cell_is_apoptotic = mpCell->HasCellProperty<ApoptoticCellProperty>();

    if (!cell_is_apoptotic)
    {
        UpdateHypoxicDuration();

        // Get cell's oxygen concentration
        double oxygen_concentration;
        switch (mDimension)
        {
            case 1:
            {
                const unsigned DIM = 1;
                oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
                break;
            }
            case 2:
            {
                const unsigned DIM = 2;
                oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
                break;
            }
            case 3:
            {
                const unsigned DIM = 3;
                oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
                break;
            }

            default:
                NEVER_REACHED;
        }

        AbstractSimpleCellCycleModel::UpdateCellCyclePhase();

        if (mCurrentCellCyclePhase == G_ONE_PHASE)
        {
            // Update G1 duration based on oxygen concentration
            double dt = SimulationTime::Instance()->GetTimeStep();

            if (oxygen_concentration < mQuiescentConcentration)
            {
                mG1Duration += (1 - std::max(oxygen_concentration, 0.0)/mQuiescentConcentration)*dt;
                mTimeSpentInG1Phase += dt;
            }
        }
    }
}

AbstractCellCycleModel* SimpleOxygenBasedCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();

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
    p_model->SetHypoxicConcentration(mHypoxicConcentration);
    p_model->SetQuiescentConcentration(mQuiescentConcentration);
    p_model->SetCriticalHypoxicDuration(mCriticalHypoxicDuration);
    p_model->SetCurrentHypoxiaOnsetTime(mCurrentHypoxiaOnsetTime);

    return p_model;
}

void SimpleOxygenBasedCellCycleModel::UpdateHypoxicDuration()
{
    assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());

    // Get cell's oxygen concentration
    double oxygen_concentration;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        default:
            NEVER_REACHED;
    }

    if (oxygen_concentration < mHypoxicConcentration)
    {
        // Update the duration of the current period of hypoxia
        mCurrentHypoxicDuration = (SimulationTime::Instance()->GetTime() - mCurrentHypoxiaOnsetTime);

        // Include a little bit of stochasticity here
        double prob_of_death = 0.9 - 0.5*(oxygen_concentration/mHypoxicConcentration);
        if (mCurrentHypoxicDuration > mCriticalHypoxicDuration && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
        {
            mpCell->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
        }
    }
    else
    {
        // Reset the cell's hypoxic duration and update the time at which the onset of hypoxia occurs
        mCurrentHypoxicDuration = 0.0;
        mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
    }
}

double SimpleOxygenBasedCellCycleModel::GetHypoxicConcentration()
{
    return mHypoxicConcentration;
}

void SimpleOxygenBasedCellCycleModel::SetHypoxicConcentration(double hypoxicConcentration)
{
    assert(hypoxicConcentration<=1.0);
    assert(hypoxicConcentration>=0.0);
    mHypoxicConcentration = hypoxicConcentration;
}

double SimpleOxygenBasedCellCycleModel::GetQuiescentConcentration()
{
    return mQuiescentConcentration;
}

void SimpleOxygenBasedCellCycleModel::SetQuiescentConcentration(double quiescentConcentration)
{
    assert(quiescentConcentration <= 1.0);
    assert(quiescentConcentration >= 0.0);
    mQuiescentConcentration = quiescentConcentration;
}

double SimpleOxygenBasedCellCycleModel::GetCriticalHypoxicDuration()
{
    return mCriticalHypoxicDuration;
}

void SimpleOxygenBasedCellCycleModel::SetCriticalHypoxicDuration(double criticalHypoxicDuration)
{
    assert(criticalHypoxicDuration >= 0.0);
    mCriticalHypoxicDuration = criticalHypoxicDuration;
}

void SimpleOxygenBasedCellCycleModel::SetCurrentHypoxiaOnsetTime(double currentHypoxiaOnsetTime)
{
    assert(currentHypoxiaOnsetTime >= 0.0);
    mCurrentHypoxiaOnsetTime = currentHypoxiaOnsetTime;
}

void SimpleOxygenBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<HypoxicConcentration>" << mHypoxicConcentration << "</HypoxicConcentration>\n";
    *rParamsFile << "\t\t\t<QuiescentConcentration>" << mQuiescentConcentration << "</QuiescentConcentration>\n";
    *rParamsFile << "\t\t\t<CriticalHypoxicDuration>" << mCriticalHypoxicDuration << "</CriticalHypoxicDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)
