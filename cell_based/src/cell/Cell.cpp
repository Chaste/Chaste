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

#include "Cell.hpp"

#include "ApoptoticCellProperty.hpp"
#include "CellAncestor.hpp"
#include "CellId.hpp"
#include "CellLabel.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "NullSrnModel.hpp"
#include "SmartPointers.hpp"

Cell::Cell(boost::shared_ptr<AbstractCellProperty> pMutationState,
           AbstractCellCycleModel* pCellCycleModel,
           AbstractSrnModel* pSrnModel,
           bool archiving,
           CellPropertyCollection cellPropertyCollection)
    : mCanDivide(false),
      mCellPropertyCollection(cellPropertyCollection),
      mpCellCycleModel(pCellCycleModel),
      mpSrnModel(pSrnModel),
      mDeathTime(DBL_MAX), // This has to be initialised for archiving
      mStartOfApoptosisTime(DBL_MAX),
      mApoptosisTime(0.25), // cell takes 15 min to fully undergo apoptosis
      mUndergoingApoptosis(false),
      mIsDead(false),
      mIsLogged(false)
{
    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        EXCEPTION("Cell is setting up a cell-cycle model but SimulationTime has not been set up");
    }

    if (pCellCycleModel == nullptr)
    {
        EXCEPTION("Cell-cycle model is null");
    }

    mpCellCycleModel->SetCell(CellPtr(this, null_deleter()));

    // Create a null srn model if none given
    if (pSrnModel == nullptr)
    {
        pSrnModel = new NullSrnModel;
        mpSrnModel = pSrnModel;
    }

    mpSrnModel->SetCell(CellPtr(this, null_deleter()));

    if (!mCellPropertyCollection.HasPropertyType<CellId>())
    {
        // Set cell identifier this will be called all the time unless the constructor is called through archiving
        MAKE_PTR(CellId, p_cell_id);
        p_cell_id->AssignCellId();
        mCellPropertyCollection.AddProperty(p_cell_id);
    }

    if (!pMutationState->IsSubType<AbstractCellMutationState>())
    {
        EXCEPTION("Attempting to create cell with a cell mutation state that is not a subtype of AbstractCellMutationState");
    }

    if (!mCellPropertyCollection.HasProperty(pMutationState))
    {
        mCellPropertyCollection.AddProperty(pMutationState);
    }

    if (!mCellPropertyCollection.HasPropertyType<CellData>())
    {
        // Add empty cell data
        MAKE_PTR(CellData, p_cell_data);
        mCellPropertyCollection.AddProperty(p_cell_data);
    }

    /*
     * If a cell proliferative type was not passed in via the input
     * argument cellPropertyCollection (for example as in the case
     * of a daughter cell being created following division) then add
     * add a 'default' cell proliferative type to the cell property
     * collection. This ensures that the method GetCellProliferativeType()
     * always returns a valid proliferative type.
     */
    if (!mCellPropertyCollection.HasPropertyType<AbstractCellProliferativeType>())
    {
        mCellPropertyCollection.AddProperty(CellPropertyRegistry::Instance()->Get<DefaultCellProliferativeType>());
    }

    if (!archiving)
    {
        // Increment cell count for each cell property in mCellPropertyCollection
        for (CellPropertyCollection::Iterator property_iter = mCellPropertyCollection.Begin();
             property_iter != mCellPropertyCollection.End();
             ++property_iter)
        {
            (*property_iter)->IncrementCellCount();
        }
    }
}

Cell::~Cell()
{
    if (!mIsDead)
    {
        Kill();
    }
    delete mpCellCycleModel;
    delete mpSrnModel;
}

void Cell::SetCellProliferativeType(boost::shared_ptr<AbstractCellProperty> pProliferativeType)
{
    if (!pProliferativeType->IsSubType<AbstractCellProliferativeType>())
    {
        EXCEPTION("Attempting to give cell a cell proliferative type that is not a subtype of AbstractCellProliferativeType");
    }

    boost::shared_ptr<AbstractCellProliferativeType> p_old_proliferative_type = GetCellProliferativeType();

    p_old_proliferative_type->DecrementCellCount();
    mCellPropertyCollection.RemoveProperty(p_old_proliferative_type);

    AddCellProperty(pProliferativeType);
}

boost::shared_ptr<AbstractCellProliferativeType> Cell::GetCellProliferativeType()  const
{
    CellPropertyCollection proliferative_type_collection = mCellPropertyCollection.GetPropertiesType<AbstractCellProliferativeType>();

    /*
     * Note: In its current form the code requires each cell to have exactly
     * one proliferative type. This is reflected in the assertion below. If a user
     * wishes to include cells with multiple proliferative types, each possible
     * combination must be created as a separate proliferative type class.
     */
    assert(proliferative_type_collection.GetSize() == 1);

    return boost::static_pointer_cast<AbstractCellProliferativeType>(proliferative_type_collection.GetProperty());
}

void Cell::SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel)
{
    if (mpCellCycleModel != pCellCycleModel)
    {
        delete mpCellCycleModel;
    }
    mpCellCycleModel = pCellCycleModel;
    mpCellCycleModel->SetCell(CellPtr(this, null_deleter()));
}

AbstractCellCycleModel* Cell::GetCellCycleModel() const
{
    return mpCellCycleModel;
}

void Cell::InitialiseCellCycleModel()
{
    mpCellCycleModel->Initialise();
}

void Cell::SetSrnModel(AbstractSrnModel* pSrnModel)
{
    if (mpSrnModel != pSrnModel)
    {
        delete mpSrnModel;
    }
    mpSrnModel = pSrnModel;
    mpSrnModel->SetCell(CellPtr(this, null_deleter()));
}

AbstractSrnModel* Cell::GetSrnModel() const
{
    return mpSrnModel;
}

void Cell::InitialiseSrnModel()
{
    mpSrnModel->Initialise();
}

double Cell::GetAge() const
{
    return mpCellCycleModel->GetAge();
}

double Cell::GetBirthTime() const
{
    return mpCellCycleModel->GetBirthTime();
}

void Cell::SetBirthTime(double birthTime)
{
    mpCellCycleModel->SetBirthTime(birthTime);
}

void Cell::SetMutationState(boost::shared_ptr<AbstractCellProperty> pMutationState)
{
    if (!pMutationState->IsSubType<AbstractCellMutationState>())
    {
        EXCEPTION("Attempting to give cell a cell mutation state that is not a subtype of AbstractCellMutationState");
    }

    boost::shared_ptr<AbstractCellMutationState> p_old_mutation_state = GetMutationState();
    p_old_mutation_state->DecrementCellCount();
    mCellPropertyCollection.RemoveProperty(p_old_mutation_state);

    AddCellProperty(pMutationState);
}

boost::shared_ptr<AbstractCellMutationState> Cell::GetMutationState() const
{
    CellPropertyCollection mutation_state_collection = mCellPropertyCollection.GetPropertiesType<AbstractCellMutationState>();

    /*
     * Note: In its current form the code requires each cell to have exactly
     * one mutation state. This is reflected in the assertion below. If a user
     * wishes to include cells with multiple mutation states, each possible
     * combination must be created as a separate mutation state class.
     */
    assert(mutation_state_collection.GetSize() == 1);

    return boost::static_pointer_cast<AbstractCellMutationState>(mutation_state_collection.GetProperty());
}

boost::shared_ptr<CellData> Cell::GetCellData() const
{
    CellPropertyCollection cell_data_collection = mCellPropertyCollection.GetPropertiesType<CellData>();

    /*
     * Note: In its current form the code requires each cell to have exactly
     * one CellData object. This is reflected in the assertion below.
     */
    assert(cell_data_collection.GetSize() <= 1);

    return boost::static_pointer_cast<CellData>(cell_data_collection.GetProperty());
}

bool Cell::HasCellVecData() const
{
    return mCellPropertyCollection.HasPropertyType<CellVecData>();
}

boost::shared_ptr<CellVecData> Cell::GetCellVecData() const
{
    assert(HasCellVecData());

    CellPropertyCollection cell_data_collection = mCellPropertyCollection.GetPropertiesType<CellVecData>();

    /*
     * Note: In its current form the code requires each cell to have exactly
     * one CellVecData object. This is reflected in the assertion below.
     */
    assert(cell_data_collection.GetSize() <= 1);

    return boost::static_pointer_cast<CellVecData>(cell_data_collection.GetProperty());
}

CellPropertyCollection& Cell::rGetCellPropertyCollection()
{
    return mCellPropertyCollection;
}

const CellPropertyCollection& Cell::rGetCellPropertyCollection() const
{
    return mCellPropertyCollection;
}

void Cell::AddCellProperty(const boost::shared_ptr<AbstractCellProperty>& rProperty)
{
    // Note: if the cell already has the specified property, no action is taken
    if (!mCellPropertyCollection.HasProperty(rProperty))
    {
        mCellPropertyCollection.AddProperty(rProperty);
        rProperty->IncrementCellCount();
    }
}

void Cell::SetLogged()
{
    mIsLogged = true;
}

bool Cell::IsLogged()
{
    return mIsLogged;
}

void Cell::StartApoptosis(bool setDeathTime)
{
    assert(!IsDead());

    if (mUndergoingApoptosis)
    {
        EXCEPTION("StartApoptosis() called when already undergoing apoptosis");
    }
    mUndergoingApoptosis = true;
    mStartOfApoptosisTime = SimulationTime::Instance()->GetTime();
    if (setDeathTime)
    {
        mDeathTime = mStartOfApoptosisTime + mApoptosisTime;
    }
    else
    {
        mDeathTime = DBL_MAX;
    }
    AddCellProperty(mCellPropertyCollection.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
}

bool Cell::HasApoptosisBegun() const
{
    return mUndergoingApoptosis;
}

double Cell::GetStartOfApoptosisTime() const
{
    return mStartOfApoptosisTime;
}

double Cell::GetApoptosisTime() const
{
    return mApoptosisTime;
}

void Cell::SetApoptosisTime(double apoptosisTime)
{
    assert(apoptosisTime > 0.0);
    mApoptosisTime = apoptosisTime;
}

double Cell::GetTimeUntilDeath() const
{
    if (!mUndergoingApoptosis || mDeathTime==DBL_MAX)
    {
        EXCEPTION("Shouldn't be checking time until apoptosis as it isn't set");
    }

    return mDeathTime - SimulationTime::Instance()->GetTime();
}

bool Cell::IsDead()
{
    if (mUndergoingApoptosis && !mIsDead)
    {
        double sloppy_death_time = mDeathTime - DBL_EPSILON * mApoptosisTime;
        if (SimulationTime::Instance()->GetTime() >= sloppy_death_time )
        {
            this->Kill();
        }
    }
    return mIsDead;
}

void Cell::Kill()
{
    // Decrement cell count for each cell property in mCellPropertyCollection
    for (CellPropertyCollection::Iterator property_iter = mCellPropertyCollection.Begin();
         property_iter != mCellPropertyCollection.End();
         ++property_iter)
    {
        (*property_iter)->DecrementCellCount();
    }
    mIsDead = true;
}

void Cell::SetAncestor(boost::shared_ptr<AbstractCellProperty> pCellAncestor)
{
    if (!pCellAncestor->IsSubType<CellAncestor>())
    {
        EXCEPTION("Attempting to give cell a cell ancestor which is not a CellAncestor");
    }

    // You can only set ancestors once
    CellPropertyCollection ancestor_collection = mCellPropertyCollection.GetPropertiesType<CellAncestor>();
    if (ancestor_collection.GetSize() == 0)
    {
        AddCellProperty(pCellAncestor);
    }
    else
    {
        // Overwrite the CellAncestor
        RemoveCellProperty<CellAncestor>();
        AddCellProperty(pCellAncestor);
    }
}

unsigned Cell::GetAncestor() const
{
    CellPropertyCollection ancestor_collection = mCellPropertyCollection.GetPropertiesType<CellAncestor>();

    assert(ancestor_collection.GetSize() <= 1);
    if (ancestor_collection.GetSize() == 0)
    {
        return UNSIGNED_UNSET;
    }

    boost::shared_ptr<CellAncestor> p_ancestor = boost::static_pointer_cast<CellAncestor>(ancestor_collection.GetProperty());

    return p_ancestor->GetAncestor();
}

unsigned Cell::GetCellId() const
{
    CellPropertyCollection cell_id_collection = mCellPropertyCollection.GetPropertiesType<CellId>();

    assert(cell_id_collection.GetSize() == 1);

    boost::shared_ptr<CellId> p_cell_id = boost::static_pointer_cast<CellId>(cell_id_collection.GetProperty());

    return p_cell_id->GetCellId();
}

bool Cell::ReadyToDivide()
{
    assert(!IsDead());
    if (mUndergoingApoptosis || HasCellProperty<ApoptoticCellProperty>())
    {
        return false;
    }

    // NOTE - we run the SRN model here first before the CCM
    mpSrnModel->SimulateToCurrentTime();
    // This in turn runs any simulations within the CCM through ReadyToDivide();
    mCanDivide = mpCellCycleModel->ReadyToDivide();

    return mCanDivide;
}

CellPtr Cell::Divide()
{
    // Check we're allowed to divide
    assert(!IsDead());
    assert(mCanDivide);
    mCanDivide = false;

    // Reset properties of parent cell
    mpCellCycleModel->ResetForDivision();
    mpSrnModel->ResetForDivision();

    // Create copy of cell property collection to modify for daughter cell
    CellPropertyCollection daughter_property_collection = mCellPropertyCollection;

    // Remove the CellId from the daughter cell, as a new one will be assigned in the constructor
    daughter_property_collection.RemoveProperty<CellId>();

    // Copy all cell data (note we create a new object not just copying the pointer)
    assert(daughter_property_collection.HasPropertyType<CellData>());

    // Get the existing copy of the cell data and remove it from the daughter cell
    boost::shared_ptr<CellData> p_cell_data = GetCellData();
    daughter_property_collection.RemoveProperty(p_cell_data);

    // Create a new cell data object using the copy constructor and add this to the daughter cell
    MAKE_PTR_ARGS(CellData, p_daughter_cell_data, (*p_cell_data));
    daughter_property_collection.AddProperty(p_daughter_cell_data);

    // Copy all cell Vec data (note we create a new object not just copying the pointer)
    if (daughter_property_collection.HasPropertyType<CellVecData>())
    {
        // Get the existing copy of the cell data and remove it from the daughter cell
        boost::shared_ptr<CellVecData> p_cell_vec_data = GetCellVecData();
        daughter_property_collection.RemoveProperty(p_cell_vec_data);

        // Create a new cell data object using the copy constructor and add this to the daughter cell
        MAKE_PTR_ARGS(CellVecData, p_daughter_cell_vec_data, (*p_cell_vec_data));
        daughter_property_collection.AddProperty(p_daughter_cell_vec_data);
    }

    // Create daughter cell with modified cell property collection
    CellPtr p_new_cell(new Cell(GetMutationState(), mpCellCycleModel->CreateCellCycleModel(), mpSrnModel->CreateSrnModel(), false, daughter_property_collection));

    // Initialise properties of daughter cell
    p_new_cell->GetCellCycleModel()->InitialiseDaughterCell();
    p_new_cell->GetSrnModel()->InitialiseDaughterCell();

    // Set the daughter cell to inherit the apoptosis time of the parent cell
    p_new_cell->SetApoptosisTime(mApoptosisTime);

    return p_new_cell;
}
