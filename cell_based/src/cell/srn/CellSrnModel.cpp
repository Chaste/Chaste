/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "CellSrnModel.hpp"

CellSrnModel::CellSrnModel(const CellSrnModel &rModel)
    : AbstractSrnModel(rModel), 
      mpInteriorSrnModel(nullptr)
{
    /*
     * Edge SRN vector should be empty in the newly created cell. They are added in
     * VertexBasedPopulation::UpdateSrnAfterBirthOrDeath().
     * 
     * Copy interior SRN to the new cell. Interior SRN should have custom implementation 
     * of a new interior SRN creation after cell division.
     */
    if (rModel.mpInteriorSrnModel != nullptr)
    {
        this->SetInteriorSrnModel(boost::shared_ptr<AbstractSrnModel>(rModel.GetInteriorSrn()->CreateSrnModel()));
    }
    mIsEdgeBasedModel = rModel.HasEdgeModel();
    for (auto iter : rModel.mEdgeSrnModels)
    {
        mEdgeSrnModels.push_back(boost::shared_ptr<AbstractSrnModel>(iter->CreateSrnModel()));
    }
}

CellSrnModel::CellSrnModel()
{
    mpInteriorSrnModel = nullptr;
}

CellSrnModel::~CellSrnModel()
{
}

void CellSrnModel::Initialise()
{
    for (auto iter : mEdgeSrnModels)
    {
        iter->Initialise();
    }

    if (mpInteriorSrnModel != nullptr)
    {
        mpInteriorSrnModel->Initialise();
    }
}

void CellSrnModel::ResetForDivision()
{
    /*
     * Making sure that we are at the current time.
     * SimulateToCurrentTime() MUST have been called before in Cell::ReadyToDivide() method
     * so that SRN models should already be simulated up to the current time.
     */
    assert(mSimulatedToTime == SimulationTime::Instance()->GetTime());

    /*
     * Note that edge models follow different rules since the number and the state of edge SRNs
     * after cell divisions depends on local topology. That is custom behavior of edge SRN models 
     * is important to implement correctly.
     */
    for (auto iter : mEdgeSrnModels)
    {
        iter->ResetForDivision();
    }

    if (mpInteriorSrnModel != nullptr)
    {
        mpInteriorSrnModel->ResetForDivision();
    }
}

void CellSrnModel::SimulateToCurrentTime()
{
    for (auto iter : mEdgeSrnModels)
    {
        iter->SimulateToCurrentTime();
    }
    if (mpInteriorSrnModel != nullptr)
    {
        mpInteriorSrnModel->SimulateToCurrentTime();
    }

    mSimulatedToTime = mEdgeSrnModels[0]->GetSimulatedToTime();
}

AbstractSrnModel* CellSrnModel::CreateSrnModel()
{
    return new CellSrnModel(*this);
}

void CellSrnModel::AddEdgeSrn(std::vector<AbstractSrnModelPtr> edgeSrns)
{
    mIsEdgeBasedModel = true;
    mEdgeSrnModels.clear();
    for (unsigned i=0; i<edgeSrns.size(); ++i)
    {
        edgeSrns[i]->SetEdgeModelIndicator(true);
    }
    mEdgeSrnModels = edgeSrns;
}

void CellSrnModel::AddEdgeSrnModel(AbstractSrnModelPtr pEdgeSrn)
{
    mIsEdgeBasedModel = true;
    pEdgeSrn->SetEdgeModelIndicator(true);
    pEdgeSrn->SetEdgeLocalIndex(mEdgeSrnModels.size());
    mEdgeSrnModels.push_back(pEdgeSrn);
}

unsigned CellSrnModel::GetNumEdgeSrn() const
{
    return mEdgeSrnModels.size();
}

AbstractSrnModelPtr CellSrnModel::GetEdgeSrn(unsigned index) const
{
    assert(index < mEdgeSrnModels.size());
    return mEdgeSrnModels[index];
}

const std::vector<AbstractSrnModelPtr>& CellSrnModel::GetEdges() const
{
    return mEdgeSrnModels;
}

void CellSrnModel::SetInteriorSrnModel(AbstractSrnModelPtr pInteriorSrn)
{
    mpInteriorSrnModel = pInteriorSrn;
}

AbstractSrnModelPtr CellSrnModel::GetInteriorSrn() const
{
    return mpInteriorSrnModel;
}

void CellSrnModel::SetCell(CellPtr pCell)
{
    AbstractSrnModel::SetCell(pCell);

    // Make a copy of all SRN models inside the system
    for (auto iter : mEdgeSrnModels)
    {
        iter->SetCell(pCell);
    }

    if (mpInteriorSrnModel != nullptr)
    {
        mpInteriorSrnModel->SetCell(pCell);
    }
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CellSrnModel)
