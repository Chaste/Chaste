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

#include "DeltaNotchSrnEdgeModel.hpp"

DeltaNotchSrnEdgeModel::DeltaNotchSrnEdgeModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(2, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchSrnEdgeModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchSrnEdgeModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

DeltaNotchSrnEdgeModel::DeltaNotchSrnEdgeModel(const DeltaNotchSrnEdgeModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    assert(rModel.GetOdeSystem());
    AbstractOdeSystem* p_parent_system(rModel.GetOdeSystem());
    SetOdeSystem(new DeltaNotchEdgeOdeSystem(p_parent_system->rGetStateVariables()));
    for (unsigned int i=0; i < p_parent_system->GetNumberOfParameters(); ++i)
        mpOdeSystem->SetParameter(i, p_parent_system->GetParameter(i));
}

AbstractSrnModel* DeltaNotchSrnEdgeModel::CreateSrnModel()
{
    return new DeltaNotchSrnEdgeModel(*this);
}

void DeltaNotchSrnEdgeModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateDeltaNotch();
    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void DeltaNotchSrnEdgeModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new DeltaNotchEdgeOdeSystem);
}

void DeltaNotchSrnEdgeModel::InitialiseDaughterCell()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);
    mpOdeSystem->SetStateVariable("Notch",0.0);
    mpOdeSystem->SetStateVariable("Delta",0.0);

    mpOdeSystem->SetParameter("neighbour delta",0.0);
    mpOdeSystem->SetParameter("interior delta",0.0);
    mpOdeSystem->SetParameter("interior notch",0.0);
}

void DeltaNotchSrnEdgeModel::UpdateDeltaNotch()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);
    double neigh_delta
    = mpCell->GetCellEdgeData()->GetItem("neighbour delta")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour delta", neigh_delta);

    double interior_delta = mpCell->GetCellData()->GetItem("interior delta");
    mpOdeSystem->SetParameter("interior delta", interior_delta);
    double interior_notch = mpCell->GetCellData()->GetItem("interior notch");
    mpOdeSystem->SetParameter("interior notch", interior_notch);
}

double DeltaNotchSrnEdgeModel::GetNotch()
{
    assert(mpOdeSystem != nullptr);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

void DeltaNotchSrnEdgeModel::SetNotch(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[0] = value;
}

double DeltaNotchSrnEdgeModel::GetDelta()
{
    assert(mpOdeSystem != nullptr);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

void DeltaNotchSrnEdgeModel::SetDelta(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[1] = value;
}

double DeltaNotchSrnEdgeModel::GetNeighbouringDelta() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour delta");
}

double DeltaNotchSrnEdgeModel::GetInteriorDelta() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("interior delta");
}

double DeltaNotchSrnEdgeModel::GetInteriorNotch() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("interior notch");
}


void DeltaNotchSrnEdgeModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

void DeltaNotchSrnEdgeModel::AddSrnQuantities(AbstractSrnModel *p_other_srn,
                                              const double scale)
{
    auto other_srn
    = static_cast<DeltaNotchSrnEdgeModel*>(p_other_srn);
    const double other_delta = other_srn->GetDelta();
    const double other_notch = other_srn->GetNotch();
    const double this_delta = GetDelta();
    const double this_notch = GetNotch();
    SetDelta(this_delta+scale*other_delta);
    SetNotch(this_notch+scale*other_notch);
}

void DeltaNotchSrnEdgeModel::AddShrunkEdgeSrn(AbstractSrnModel *p_shrunk_edge_srn)
{
    // Here we assume that one half of srn quantities are endocytosed and the remaining
    // half are split between neighbouring junctions. Hence we add 1/4 of srn variables
    AddSrnQuantities(p_shrunk_edge_srn, 0.25);
}

void DeltaNotchSrnEdgeModel::AddMergedEdgeSrn(AbstractSrnModel* p_merged_edge_srn)
{
    // Add all srn variables to this edge srn
    AddSrnQuantities(p_merged_edge_srn);
}

void DeltaNotchSrnEdgeModel::SplitEdgeSrn(const double relative_position)
{
    ScaleSrnVariables(relative_position);
}


// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchSrnEdgeModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchSrnEdgeModel)
