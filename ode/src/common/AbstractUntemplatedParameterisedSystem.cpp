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

#include "AbstractUntemplatedParameterisedSystem.hpp"


#include "Exception.hpp"
#include "VectorHelperFunctions.hpp"


AbstractUntemplatedParameterisedSystem::AbstractUntemplatedParameterisedSystem(unsigned numberOfStateVariables)
    : mNumberOfStateVariables(numberOfStateVariables)
{
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

//
// State variable methods
//

unsigned AbstractUntemplatedParameterisedSystem::GetNumberOfStateVariables() const
{
    return mNumberOfStateVariables;
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
// Parameter methods
//

unsigned AbstractUntemplatedParameterisedSystem::GetNumberOfParameters() const
{
    return GetSystemInformation()->rGetParameterNames().size();
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

unsigned AbstractUntemplatedParameterisedSystem::GetAnyVariableIndex(const std::string& rName) const
{
    return GetSystemInformation()->GetAnyVariableIndex(rName);
}

bool AbstractUntemplatedParameterisedSystem::HasAnyVariable(const std::string& rName) const
{
    return GetSystemInformation()->HasAnyVariable(rName);
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
