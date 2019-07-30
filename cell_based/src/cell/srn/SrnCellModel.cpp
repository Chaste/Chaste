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

#include "SrnCellModel.hpp"

SrnCellModel::SrnCellModel(const SrnCellModel &rModel)
    : AbstractSrnModel(rModel), mInteriorSrnModel(nullptr)
{
    // Make a copy of all SRN models inside the system
    for (auto srnModel: rModel.mEdgeSrnModels)
    {
        this->AddEdgeSrn(boost::shared_ptr<AbstractSrnModel>(srnModel->CreateSrnModel()));
    }
    if (rModel.mInteriorSrnModel!=nullptr)
        this->SetInteriorSrnModel(boost::shared_ptr<AbstractSrnModel>(rModel.mInteriorSrnModel->CreateSrnModel()));
}


SrnCellModel::SrnCellModel()
{
    mInteriorSrnModel = nullptr;

}

SrnCellModel::~SrnCellModel()
{

}

void SrnCellModel::Initialise()
{
    for (auto edgeModel : mEdgeSrnModels)
    {
        edgeModel->Initialise();
    }

    if (mInteriorSrnModel != nullptr)
    {
        mInteriorSrnModel->Initialise();
    }
}

void SrnCellModel::SimulateToCurrentTime()
{
    for (auto srnModel : mEdgeSrnModels)
    {
        srnModel->SimulateToCurrentTime();
    }
    if (mInteriorSrnModel != nullptr)
    {
        mInteriorSrnModel->SimulateToCurrentTime();
    }
}

AbstractSrnModel* SrnCellModel::CreateSrnModel()
{
    return new SrnCellModel(*this);
}

void SrnCellModel::AddEdgeSrn(std::vector<AbstractSrnModelPtr> edgeSrn)
{
    mIsEdgeBasedModel = true;
    for (unsigned int i=0; i<edgeSrn.size(); ++i)
        edgeSrn[i]->SetEdgeModelIndicator(true);
    mEdgeSrnModels = edgeSrn;
}

void SrnCellModel::AddEdgeSrn(AbstractSrnModelPtr edgeSrn)
{
    mIsEdgeBasedModel = true;
    edgeSrn->SetEdgeModelIndicator(true);
    edgeSrn->SetEdgeLocalIndex(mEdgeSrnModels.size());
    mEdgeSrnModels.push_back(edgeSrn);
}

void SrnCellModel::InsertEdgeSrn(unsigned index, AbstractSrnModelPtr edgeSrn)
{
    mIsEdgeBasedModel = true;
    edgeSrn->SetEdgeModelIndicator(true);
    mEdgeSrnModels.insert(mEdgeSrnModels.begin() + index, edgeSrn);
}

AbstractSrnModelPtr SrnCellModel::RemoveEdgeSrn(unsigned index)
{
    auto edgeSrn = mEdgeSrnModels[index];
    mEdgeSrnModels.erase(mEdgeSrnModels.begin() + index);
    return edgeSrn;
}

unsigned SrnCellModel::GetNumEdgeSrn()
{
    return mEdgeSrnModels.size();
}

AbstractSrnModelPtr SrnCellModel::GetEdgeSrn(unsigned index)
{
    assert(index < mEdgeSrnModels.size());
    return mEdgeSrnModels[index];
}

const std::vector<AbstractSrnModelPtr>& SrnCellModel::GetEdges()
{
    return mEdgeSrnModels;
}

void SrnCellModel::SetInteriorSrnModel(AbstractSrnModelPtr interiorSrn)
{
    mInteriorSrnModel = interiorSrn;
}

AbstractSrnModelPtr SrnCellModel::GetInteriorSrn()
{
    return mInteriorSrnModel;
}

void SrnCellModel::SetCell(CellPtr pCell)
{
    AbstractSrnModel::SetCell(pCell);

    // Make a copy of all SRN models inside the system
    for (auto srnModel : mEdgeSrnModels)
    {
        srnModel->SetCell(pCell);
    }

    if (mInteriorSrnModel != nullptr)
    {
        mInteriorSrnModel->SetCell(pCell);
    }
}



// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SrnCellModel)
