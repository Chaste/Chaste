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

#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "ApoptoticCellProperty.hpp"

SimpleOxygenBasedCellCycleModel::SimpleOxygenBasedCellCycleModel()
    : mCurrentHypoxicDuration(0.0),
      mHypoxicConcentration(0.4),
      mQuiescentConcentration(1.0),
      mCriticalHypoxicDuration(2.0)
{
    mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
}

double SimpleOxygenBasedCellCycleModel::GetCurrentHypoxicDuration() const
{
    return mCurrentHypoxicDuration;
}

double SimpleOxygenBasedCellCycleModel::GetCurrentHypoxiaOnsetTime() const
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
        double oxygen_concentration = mpCell->GetCellData()->GetItem("oxygen");

        AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();

        if (mCurrentCellCyclePhase == G_ONE_PHASE)
        {
            // Update G1 duration based on oxygen concentration
            double dt = SimulationTime::Instance()->GetTimeStep();

            if (oxygen_concentration < mQuiescentConcentration)
            {
                mG1Duration += (1 - std::max(oxygen_concentration, 0.0)/mQuiescentConcentration)*dt;
            }
        }
    }
}

SimpleOxygenBasedCellCycleModel::SimpleOxygenBasedCellCycleModel(const SimpleOxygenBasedCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModel(rModel),
     mCurrentHypoxicDuration(rModel.mCurrentHypoxicDuration),
     mCurrentHypoxiaOnsetTime(rModel.mCurrentHypoxiaOnsetTime),
     mHypoxicConcentration(rModel.mHypoxicConcentration),
     mQuiescentConcentration(rModel.mQuiescentConcentration),
     mCriticalHypoxicDuration(rModel.mCriticalHypoxicDuration)
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
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* SimpleOxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new SimpleOxygenBasedCellCycleModel(*this);
}

void SimpleOxygenBasedCellCycleModel::UpdateHypoxicDuration()
{
    assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());

    // Get cell's oxygen concentration
    double oxygen_concentration = mpCell->GetCellData()->GetItem("oxygen");

    if (oxygen_concentration < mHypoxicConcentration)
    {
        // Update the duration of the current period of hypoxia
        mCurrentHypoxicDuration = (SimulationTime::Instance()->GetTime() - mCurrentHypoxiaOnsetTime);

        // Include a little bit of stochasticity here
        double prob_of_death = 0.9 - 0.5*(oxygen_concentration/mHypoxicConcentration);
        if (mCurrentHypoxicDuration > mCriticalHypoxicDuration && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
        {
            /*
             * This method is usually called within a CellBasedSimulation, after the CellPopulation
             * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
             * CellPropertyRegistry::Instance() here when adding the ApoptoticCellProperty, we would
             * be creating a new CellPropertyRegistry. In this case the ApoptoticCellProperty cell
             * count would be incorrect. We must therefore access the ApoptoticCellProperty via the
             * cell's CellPropertyCollection.
             */
            boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
            mpCell->AddCellProperty(p_apoptotic_property);
        }
    }
    else
    {
        // Reset the cell's hypoxic duration and update the time at which the onset of hypoxia occurs
        mCurrentHypoxicDuration = 0.0;
        mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
    }
}

double SimpleOxygenBasedCellCycleModel::GetHypoxicConcentration() const
{
    return mHypoxicConcentration;
}

void SimpleOxygenBasedCellCycleModel::SetHypoxicConcentration(double hypoxicConcentration)
{
    assert(hypoxicConcentration<=1.0);
    assert(hypoxicConcentration>=0.0);
    mHypoxicConcentration = hypoxicConcentration;
}

double SimpleOxygenBasedCellCycleModel::GetQuiescentConcentration() const
{
    return mQuiescentConcentration;
}

void SimpleOxygenBasedCellCycleModel::SetQuiescentConcentration(double quiescentConcentration)
{
    assert(quiescentConcentration <= 1.0);
    assert(quiescentConcentration >= 0.0);
    mQuiescentConcentration = quiescentConcentration;
}

double SimpleOxygenBasedCellCycleModel::GetCriticalHypoxicDuration() const
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
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)
