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

#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"

VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo::VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractVanLeeuwen2009WntSwatCellCycleModel(pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.00005);
#endif //CHASTE_CVODE
    }
}

VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo::VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo(const VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo& rModel)
   : AbstractVanLeeuwen2009WntSwatCellCycleModel(rModel)
{
    /*
     * Initialize only those member variables defined in this class.
     * Create the new cell-cycle model's ODE system and use the current
     * values of the state variables in mpOdeSystem as an initial condition.
     *
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
    assert(rModel.GetOdeSystem());
    double wnt_level = rModel.GetWntLevel();
    InitialiseOdeSystem(wnt_level, rModel.mpCell->GetMutationState());
    SetStateVariables(rModel.GetOdeSystem()->rGetStateVariables());
}

void VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo::InitialiseOdeSystem(double wntConcentration, boost::shared_ptr<AbstractCellMutationState> pMutationState)
{
    mpOdeSystem = new VanLeeuwen2009WntSwatCellCycleOdeSystem(2, wntConcentration,  pMutationState);
}

AbstractCellCycleModel* VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo::CreateCellCycleModel()
{
    return new VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo(*this);
}

void VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractVanLeeuwen2009WntSwatCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo)
