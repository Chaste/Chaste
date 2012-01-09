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

#include <sstream>
#include <cassert>

#include "AbstractParameterisedSystem.hpp"

#include "Exception.hpp"
#include "VectorHelperFunctions.hpp"


AbstractUntemplatedParameterisedSystem::AbstractUntemplatedParameterisedSystem(unsigned numberOfStateVariables)
    : mNumberOfStateVariables(numberOfStateVariables)
{
}


template<typename VECTOR>
AbstractParameterisedSystem<VECTOR>::AbstractParameterisedSystem(unsigned numberOfStateVariables)
    : AbstractUntemplatedParameterisedSystem(numberOfStateVariables)
{
    InitialiseEmptyVector(mParameters);
    InitialiseEmptyVector(mStateVariables);
}


AbstractUntemplatedParameterisedSystem::~AbstractUntemplatedParameterisedSystem()
{
}


boost::shared_ptr<const AbstractOdeSystemInformation> AbstractUntemplatedParameterisedSystem::GetSystemInformation() const
{
    assert(mpSystemInfo);
    return mpSystemInfo;
}


std::string AbstractUntemplatedParameterisedSystem::GetSystemName() const
{
    return GetSystemInformation()->GetSystemName();
}

template<typename VECTOR>
std::string AbstractParameterisedSystem<VECTOR>::DumpState(const std::string& message)
{
    return GetStateMessage(message, mStateVariables);
}

template<typename VECTOR>
std::string AbstractParameterisedSystem<VECTOR>::DumpState(const std::string& message,
                                                           VECTOR Y)
{
    return GetStateMessage(message, Y);
}

template<typename VECTOR>
std::string AbstractParameterisedSystem<VECTOR>::GetStateMessage(const std::string& message, VECTOR Y)
{
    std::stringstream res;
    res << message << "\nState:\n";
    assert(rGetStateVariableNames().size()==GetVectorSize(Y));
    for (unsigned i=0; i<GetVectorSize(Y); i++)
    {
        res << "\t" << rGetStateVariableNames()[i] << ":" << GetVectorComponent(Y, i) << "\n";
    }
    return res.str();
}

//
// State variable methods
//

unsigned AbstractUntemplatedParameterisedSystem::GetNumberOfStateVariables() const
{
    return mNumberOfStateVariables;
}

template<typename VECTOR>
VECTOR& AbstractParameterisedSystem<VECTOR>::rGetStateVariables()
{
    return mStateVariables;
}

template<typename VECTOR>
VECTOR AbstractParameterisedSystem<VECTOR>::GetStateVariables()
{
    return CopyVector(mStateVariables);
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetStateVariables(const VECTOR& rStateVariables)
{
    if ( mNumberOfStateVariables != GetVectorSize(rStateVariables) )
    {
        EXCEPTION("The size of the passed in vector must be that of the number of state variables.");
    }
    CreateVectorIfEmpty(mStateVariables, mNumberOfStateVariables);
    for (unsigned i=0; i<mNumberOfStateVariables; i++)
    {
        SetVectorComponent(mStateVariables, i, GetVectorComponent(rStateVariables, i));
    }
}

template<typename VECTOR>
double AbstractParameterisedSystem<VECTOR>::GetStateVariable(unsigned index) const
{
    if (index >= mNumberOfStateVariables)
    {
        EXCEPTION("The index passed in must be less than the number of state variables.");
    }
    return GetVectorComponent(mStateVariables, index);
}

template<typename VECTOR>
double AbstractParameterisedSystem<VECTOR>::GetStateVariable(const std::string& rName) const
{
    return GetStateVariable(GetStateVariableIndex(rName));
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetStateVariable(unsigned index, double newValue)
{
    if ( mNumberOfStateVariables <= index )
    {
        EXCEPTION("The index passed in must be less than the number of state variables.");
    }
    SetVectorComponent(mStateVariables, index, newValue);
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetStateVariable(const std::string& rName, double newValue)
{
    SetStateVariable(GetStateVariableIndex(rName), newValue);
}


const std::vector<std::string>& AbstractUntemplatedParameterisedSystem::rGetStateVariableNames() const
{
    return GetSystemInformation()->rGetStateVariableNames();
}

const std::vector<std::string>& AbstractUntemplatedParameterisedSystem::rGetStateVariableUnits() const
{
    return GetSystemInformation()->rGetStateVariableUnits();
}

unsigned AbstractUntemplatedParameterisedSystem::GetStateVariableIndex(const std::string& rName) const
{
    return GetSystemInformation()->GetStateVariableIndex(rName);
}

bool AbstractUntemplatedParameterisedSystem::HasStateVariable(const std::string& rName) const
{
    return GetSystemInformation()->HasStateVariable(rName);
}

std::string AbstractUntemplatedParameterisedSystem::GetStateVariableUnits(unsigned index) const
{
    return GetSystemInformation()->GetStateVariableUnits(index);
}

//
// Initial condition methods
//

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetDefaultInitialConditions(const VECTOR& rInitialConditions)
{
    if (GetVectorSize(rInitialConditions) != mNumberOfStateVariables)
    {
        EXCEPTION("The number of initial conditions must be that of the number of state variables.");
    }
    assert(mpSystemInfo);
    std::vector<double> inits;
    CopyToStdVector(rInitialConditions, inits);
    mpSystemInfo->SetDefaultInitialConditions(inits);
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetDefaultInitialCondition(unsigned index, double initialCondition)
{
    if (index >= mNumberOfStateVariables)
    {
        EXCEPTION("Index is greater than the number of state variables.");
    }
    assert(mpSystemInfo);
    mpSystemInfo->SetDefaultInitialCondition(index, initialCondition);
}

template<typename VECTOR>
VECTOR AbstractParameterisedSystem<VECTOR>::GetInitialConditions() const
{
    assert(mpSystemInfo);
    VECTOR v;
    InitialiseEmptyVector(v);
    CreateVectorIfEmpty(v, mNumberOfStateVariables);
    CopyFromStdVector(mpSystemInfo->GetInitialConditions(), v);
    return v;
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::ResetToInitialConditions()
{
    VECTOR inits = GetInitialConditions();
    SetStateVariables(inits);
    DeleteVector(inits);
}

//
// Parameter methods
//

unsigned AbstractUntemplatedParameterisedSystem::GetNumberOfParameters() const
{
    return GetSystemInformation()->rGetParameterNames().size();
}

template<typename VECTOR>
double AbstractParameterisedSystem<VECTOR>::GetParameter(unsigned index) const
{
    if (index >= GetVectorSize(mParameters))
    {
        EXCEPTION("The index passed in must be less than the number of parameters.");
    }
    return GetVectorComponent(mParameters, index);
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetParameter(unsigned index, double value)
{
    if (index >= GetVectorSize(mParameters))
    {
        EXCEPTION("The index passed in must be less than the number of parameters.");
    }
    SetVectorComponent(mParameters, index, value);
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetParameter(const std::string& rName, double value)
{
    SetVectorComponent(mParameters, GetParameterIndex(rName), value);
}

template<typename VECTOR>
double AbstractParameterisedSystem<VECTOR>::GetParameter(const std::string& rName) const
{
    return GetParameter(GetParameterIndex(rName));
}


const std::vector<std::string>& AbstractUntemplatedParameterisedSystem::rGetParameterNames() const
{
    return GetSystemInformation()->rGetParameterNames();
}

const std::vector<std::string>& AbstractUntemplatedParameterisedSystem::rGetParameterUnits() const
{
    return GetSystemInformation()->rGetParameterUnits();
}

unsigned AbstractUntemplatedParameterisedSystem::GetParameterIndex(const std::string& rName) const
{
    return GetSystemInformation()->GetParameterIndex(rName);
}

bool AbstractUntemplatedParameterisedSystem::HasParameter(const std::string& rName) const
{
    return GetSystemInformation()->HasParameter(rName);
}

std::string AbstractUntemplatedParameterisedSystem::GetParameterUnits(unsigned index) const
{
    return GetSystemInformation()->GetParameterUnits(index);
}

//
// "Any variable" methods
//

template<typename VECTOR>
double AbstractParameterisedSystem<VECTOR>::GetAnyVariable(unsigned index, double time,
                                                           VECTOR* pDerivedQuantities)
{
    if (index < mNumberOfStateVariables)
    {
        return GetVectorComponent(mStateVariables, index);
    }
    else if (index - mNumberOfStateVariables < GetVectorSize(mParameters))
    {
        return GetVectorComponent(mParameters, index - mNumberOfStateVariables);
    }
    else
    {
        unsigned offset = mNumberOfStateVariables + GetVectorSize(mParameters);
        if (index - offset < GetNumberOfDerivedQuantities())
        {
            VECTOR dqs;
            if (pDerivedQuantities == NULL)
            {
                dqs = ComputeDerivedQuantitiesFromCurrentState(time);
                pDerivedQuantities = &dqs;
            }
            double value = GetVectorComponent(*pDerivedQuantities, index - offset);
            if (pDerivedQuantities == &dqs)
            {
                DeleteVector(dqs);
            }
            return value;
        }
        else
        {
            EXCEPTION("Invalid index passed to GetAnyVariable.");
        }
    }
}

template<typename VECTOR>
double AbstractParameterisedSystem<VECTOR>::GetAnyVariable(const std::string& rName,
                                                           double time,
                                                           VECTOR* pDerivedQuantities)
{
    return GetAnyVariable(GetAnyVariableIndex(rName), time, pDerivedQuantities);
}


unsigned AbstractUntemplatedParameterisedSystem::GetAnyVariableIndex(const std::string& rName) const
{
    return GetSystemInformation()->GetAnyVariableIndex(rName);
}

bool AbstractUntemplatedParameterisedSystem::HasAnyVariable(const std::string& rName) const
{
    return GetSystemInformation()->HasAnyVariable(rName);
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetAnyVariable(unsigned index, double value)
{
    if (index < mNumberOfStateVariables)
    {
        SetVectorComponent(mStateVariables, index, value);
    }
    else if (index - mNumberOfStateVariables < GetVectorSize(mParameters))
    {
        SetVectorComponent(mParameters, index - mNumberOfStateVariables, value);
    }
    else
    {
        EXCEPTION("Cannot set the value of a derived quantity, or invalid index.");
    }
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::SetAnyVariable(const std::string& rName, double value)
{
    SetAnyVariable(GetAnyVariableIndex(rName), value);
}


std::string AbstractUntemplatedParameterisedSystem::GetAnyVariableUnits(unsigned index) const
{
    return GetSystemInformation()->GetAnyVariableUnits(index);
}

std::string AbstractUntemplatedParameterisedSystem::GetAnyVariableUnits(const std::string& rName) const
{
    return GetAnyVariableUnits(GetAnyVariableIndex(rName));
}

//
// "Derived quantities" methods
//

unsigned AbstractUntemplatedParameterisedSystem::GetNumberOfDerivedQuantities() const
{
    return GetSystemInformation()->rGetDerivedQuantityNames().size();
}

template<typename VECTOR>
VECTOR AbstractParameterisedSystem<VECTOR>::ComputeDerivedQuantities(double time,
                                                                     const VECTOR& rState)
{
    EXCEPTION("This ODE system does not define derived quantities.");
}

template<typename VECTOR>
VECTOR AbstractParameterisedSystem<VECTOR>::ComputeDerivedQuantitiesFromCurrentState(double time)
{
    return this->ComputeDerivedQuantities(time, mStateVariables);
}


const std::vector<std::string>& AbstractUntemplatedParameterisedSystem::rGetDerivedQuantityNames() const
{
    return GetSystemInformation()->rGetDerivedQuantityNames();
}

const std::vector<std::string>& AbstractUntemplatedParameterisedSystem::rGetDerivedQuantityUnits() const
{
    return GetSystemInformation()->rGetDerivedQuantityUnits();
}

unsigned AbstractUntemplatedParameterisedSystem::GetDerivedQuantityIndex(const std::string& rName) const
{
    return GetSystemInformation()->GetDerivedQuantityIndex(rName);
}

bool AbstractUntemplatedParameterisedSystem::HasDerivedQuantity(const std::string& rName) const
{
    return GetSystemInformation()->HasDerivedQuantity(rName);
}

std::string AbstractUntemplatedParameterisedSystem::GetDerivedQuantityUnits(unsigned index) const
{
    return GetSystemInformation()->GetDerivedQuantityUnits(index);
}


unsigned AbstractUntemplatedParameterisedSystem::GetNumberOfAttributes() const
{
    return GetSystemInformation()->GetNumberOfAttributes();
}

bool AbstractUntemplatedParameterisedSystem::HasAttribute(const std::string& rName) const
{
    return GetSystemInformation()->HasAttribute(rName);
}

double AbstractUntemplatedParameterisedSystem::GetAttribute(const std::string& rName) const
{
    return GetSystemInformation()->GetAttribute(rName);
}



////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class AbstractParameterisedSystem<std::vector<double> >;
#ifdef CHASTE_CVODE
template class AbstractParameterisedSystem<N_Vector>;
#endif
