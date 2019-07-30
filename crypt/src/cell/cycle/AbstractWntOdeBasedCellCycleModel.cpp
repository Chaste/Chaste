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
#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "Exception.hpp"

AbstractWntOdeBasedCellCycleModel::AbstractWntOdeBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeBasedPhaseBasedCellCycleModel(SimulationTime::Instance()->GetTime(), pOdeSolver)
{
}

AbstractWntOdeBasedCellCycleModel::~AbstractWntOdeBasedCellCycleModel()
{
}

AbstractWntOdeBasedCellCycleModel::AbstractWntOdeBasedCellCycleModel(const AbstractWntOdeBasedCellCycleModel& rModel)
   : AbstractOdeBasedPhaseBasedCellCycleModel(rModel)
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

double AbstractWntOdeBasedCellCycleModel::GetWntLevel() const
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

void AbstractWntOdeBasedCellCycleModel::ResetForDivision()
{
    AbstractOdeBasedPhaseBasedCellCycleModel::ResetForDivision();

    assert(mpOdeSystem != nullptr);

    // This model needs the protein concentrations and phase resetting to G0/G1.
    // Keep the Wnt pathway in the same state but reset the cell cycle part
    // Cell cycle is proteins 0 to 4 (first 5 ODEs)
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i=0; i<5; i++)
    {
       mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
}

void AbstractWntOdeBasedCellCycleModel::UpdateCellCyclePhase()
{
    AbstractOdeBasedPhaseBasedCellCycleModel::UpdateCellCyclePhase();
    if (SimulationTime::Instance()->GetTime() == mLastTime
        || GetOdeStopTime() == mLastTime)
    {
        // Should only be called if we're still running ODEs
        UpdateCellProliferativeType();
    }
}

void AbstractWntOdeBasedCellCycleModel::UpdateCellProliferativeType()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);
    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}

double AbstractWntOdeBasedCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 16.0;
}

double AbstractWntOdeBasedCellCycleModel::GetAverageStemCellCycleTime()
{
    return 16.0;
}

bool AbstractWntOdeBasedCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

void AbstractWntOdeBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeBasedPhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
