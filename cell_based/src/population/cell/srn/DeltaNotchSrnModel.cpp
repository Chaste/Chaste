/*

Copyright (c) 2005-2015, University of Oxford.
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

//#include "UblasIncludes.hpp"
#include "DeltaNotchSrnModel.hpp"
//#include "CellCycleModelOdeSolver.hpp"
#include "AbstractOdeSrnModel.hpp"
//#include "CvodeAdaptor.hpp"
//#include "Exception.hpp"



DeltaNotchSrnModel::DeltaNotchSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(2, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

AbstractSrnModel* DeltaNotchSrnModel::CreateSrnModel()
{
    // Create a new srn model
    DeltaNotchSrnModel* p_model = new DeltaNotchSrnModel(this->mpOdeSolver);
    // Create the new srn model's ODE system
    p_model->SetOdeSystem(new DeltaNotchOdeSystem);
    // Call super to set current values of the state variables in mpOdeSystem as an initial condition for the new srn model's ODE system
    return AbstractOdeSrnModel::CreateSrnModel(p_model);
}


void DeltaNotchSrnModel::SimulateToCurrentTime()
{
	// Custom behaviour
	UpdateDeltaNotch();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}


void DeltaNotchSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new DeltaNotchOdeSystem);
}


void DeltaNotchSrnModel::UpdateDeltaNotch()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    double mean_delta = mpCell->GetCellData()->GetItem("mean delta");
    mpOdeSystem->SetParameter("Mean Delta", mean_delta);
}

double DeltaNotchSrnModel::GetNotch()
{
    assert(mpOdeSystem != NULL);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

double DeltaNotchSrnModel::GetDelta()
{
    assert(mpOdeSystem != NULL);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

double DeltaNotchSrnModel::GetMeanNeighbouringDelta()
{
    assert(mpOdeSystem != NULL);
    double mean_neighbouring_delta = mpOdeSystem->GetParameter("Mean Delta");
    return mean_neighbouring_delta;
}

void DeltaNotchSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output.

    // Call direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchSrnModel)



