/*

Copyright (C) University of Oxford, 2005-2009

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
#include "StochasticOxygenBasedCellCycleModel.hpp"

void StochasticOxygenBasedCellCycleModel::SetG2Duration()
{
    TissueConfig* p_params = TissueConfig::Instance();
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double mean = p_params->GetG2Duration();
    double standard_deviation = 1.0;

    mG2Duration = p_gen->NormalRandomDeviate(mean, standard_deviation);

    // Check that the normal random deviate has not returned a small or negative G2 duration
    if (mG2Duration < p_params->GetMinimumGapDuration())
    {
        #define COVERAGE_IGNORE
        mG2Duration = p_params->GetMinimumGapDuration();
        #undef COVERAGE_IGNORE
    }
}


void StochasticOxygenBasedCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
    SetG2Duration();
}


void StochasticOxygenBasedCellCycleModel::Initialise()
{
    AbstractSimpleCellCycleModel::Initialise();
    SetG2Duration();
}


void StochasticOxygenBasedCellCycleModel::ResetForDivision()
{
    AbstractSimpleCellCycleModel::ResetForDivision();
    SetG2Duration();
}


double StochasticOxygenBasedCellCycleModel::GetG2Duration()
{
    return mG2Duration;
}


StochasticOxygenBasedCellCycleModel::StochasticOxygenBasedCellCycleModel(unsigned dimension)
    : mTimeSpentInG1Phase(0.0),
      mCurrentHypoxicDuration(0.0),
      mDimension(dimension)
{
    mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
}


double StochasticOxygenBasedCellCycleModel::GetCurrentHypoxicDuration()
{
    return mCurrentHypoxicDuration;
}


double StochasticOxygenBasedCellCycleModel::GetCurrentHypoxiaOnsetTime()
{
    return mCurrentHypoxiaOnsetTime;
}


void StochasticOxygenBasedCellCycleModel::UpdateCellCyclePhase()
{
    // mG1Duration is set when the cell cycle model is given a cell

    if (mpCell->GetCellProliferativeType()!=APOPTOTIC)
    {
        UpdateHypoxicDuration();

        // Get cell's oxygen concentration
        double oxygen_concentration;
        switch (mDimension)
        {
            case 1:
            {
                const unsigned DIM = 1;
                oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
                break;
            }
            case 2:
            {
                const unsigned DIM = 2;
                oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
                break;
            }
            case 3:
            {
                const unsigned DIM = 3;
                oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
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
            double quiescent_concentration = TissueConfig::Instance()->GetHepaOneCellQuiescentConcentration();

            if (oxygen_concentration < quiescent_concentration)
            {
                mG1Duration += (1 - std::max(oxygen_concentration,0.0)/quiescent_concentration)*dt;
                mTimeSpentInG1Phase += dt;
            }
        }
    }
}

AbstractCellCycleModel* StochasticOxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new StochasticOxygenBasedCellCycleModel(*this);
}

void StochasticOxygenBasedCellCycleModel::UpdateHypoxicDuration()
{
    assert(mpCell->GetCellProliferativeType()!=APOPTOTIC);
    assert(!mpCell->HasApoptosisBegun());

    // Get cell's oxygen concentration
    double oxygen_concentration;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            oxygen_concentration = CellwiseData<DIM>::Instance()->GetValue(*mpCell, 0);
            break;
        }
        default:
            NEVER_REACHED;
    }

    double hypoxic_concentration = TissueConfig::Instance()->GetHepaOneCellHypoxicConcentration();

    if ( oxygen_concentration < hypoxic_concentration)
    {
        // Update the duration of the current period of hypoxia
        mCurrentHypoxicDuration = (SimulationTime::Instance()->GetTime() - mCurrentHypoxiaOnsetTime);

        // Include a little bit of stochasticity here
        double prob_of_death = 0.9 - 0.5*(oxygen_concentration/hypoxic_concentration);
        if (mCurrentHypoxicDuration > TissueConfig::Instance()->GetCriticalHypoxicDuration() && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
        {
            mpCell->SetCellProliferativeType(APOPTOTIC);
        }
    }
    else
    {
        // Reset the cell's hypoxic duration and update the time at which the onset of hypoxia occurs
        mCurrentHypoxicDuration = 0.0;
        mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetTime();
    }
}


unsigned StochasticOxygenBasedCellCycleModel::GetDimension()
{
    return mDimension;
}
