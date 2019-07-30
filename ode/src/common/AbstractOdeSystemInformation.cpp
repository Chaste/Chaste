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


#include <cassert>
#include <algorithm>

#include "AbstractOdeSystemInformation.hpp"
#include "Exception.hpp"

AbstractOdeSystemInformation::AbstractOdeSystemInformation()
    : mInitialised(false)
{
}

AbstractOdeSystemInformation::~AbstractOdeSystemInformation()
{
}

std::string AbstractOdeSystemInformation::GetSystemName() const
{
    return mSystemName;
}

std::string AbstractOdeSystemInformation::GetFreeVariableName() const
{
    return mFreeVariableName;
}

std::string AbstractOdeSystemInformation::GetFreeVariableUnits() const
{
    return mFreeVariableUnits;
}

void AbstractOdeSystemInformation::SetDefaultInitialConditions(const std::vector<double>& rInitialConditions)
{
    assert(mInitialised);
    mInitialConditions = rInitialConditions;
}

void AbstractOdeSystemInformation::SetDefaultInitialCondition(unsigned index, double initialCondition)
{
    assert(mInitialised);
    mInitialConditions.at(index) = initialCondition;
}

std::vector<double> AbstractOdeSystemInformation::GetInitialConditions() const
{
    assert(mInitialised);
    return mInitialConditions;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetStateVariableNames() const
{
    assert(mInitialised);
    return mVariableNames;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetStateVariableUnits() const
{
    assert(mInitialised);
    return mVariableUnits;
}

unsigned AbstractOdeSystemInformation::GetStateVariableIndex(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mVariableNames.begin(), mVariableNames.end(), rName);
    if (it == mVariableNames.end())
    {
        EXCEPTION("No state variable named '" + rName + "'.");
    }
    return (unsigned)(it - mVariableNames.begin());
}

bool AbstractOdeSystemInformation::HasStateVariable(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mVariableNames.begin(), mVariableNames.end(), rName);
    return (it != mVariableNames.end());
}

std::string AbstractOdeSystemInformation::GetStateVariableUnits(unsigned index) const
{
    assert(mInitialised);
    if (index >= mVariableUnits.size())
    {
        EXCEPTION("The index passed in must be less than the number of state variables.");
    }
    return mVariableUnits[index];
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetParameterNames() const
{
    assert(mInitialised);
    return mParameterNames;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetParameterUnits() const
{
    assert(mInitialised);
    return mParameterUnits;
}

unsigned AbstractOdeSystemInformation::GetParameterIndex(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mParameterNames.begin(), mParameterNames.end(), rName);
    if (it == mParameterNames.end())
    {
        EXCEPTION("No parameter named '" + rName + "'.");
    }
    return (unsigned)(it - mParameterNames.begin());
}

bool AbstractOdeSystemInformation::HasParameter(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mParameterNames.begin(), mParameterNames.end(), rName);
    return (it != mParameterNames.end());
}

std::string AbstractOdeSystemInformation::GetParameterUnits(unsigned index) const
{
    assert(mInitialised);
    if (index >= mParameterUnits.size())
    {
        EXCEPTION("The index passed in must be less than the number of parameters.");
    }
    return mParameterUnits[index];
}

unsigned AbstractOdeSystemInformation::GetNumberOfParameters() const
{
    assert(mInitialised);
    return mParameterUnits.size();
}

unsigned AbstractOdeSystemInformation::GetAnyVariableIndex(const std::string& rName) const
{
    assert(mInitialised);
    if (HasStateVariable(rName))
    {
        return GetStateVariableIndex(rName);
    }
    else if (HasParameter(rName))
    {
        return mVariableNames.size() + GetParameterIndex(rName);
    }
    else if (HasDerivedQuantity(rName))
    {
        return mVariableNames.size() + mParameterNames.size() + GetDerivedQuantityIndex(rName);
    }
    else
    {
        EXCEPTION("No state variable, parameter, or derived quantity named '" + rName + "'.");
    }
}


bool AbstractOdeSystemInformation::HasAnyVariable(const std::string& rName) const
{
    assert(mInitialised);
    return (HasStateVariable(rName) || HasParameter(rName) || HasDerivedQuantity(rName));
}

std::string AbstractOdeSystemInformation::GetAnyVariableUnits(unsigned index) const
{
    assert(mInitialised);
    if (index < mVariableUnits.size())
    {
        return mVariableUnits[index];
    }
    else
    {
        unsigned offset = mVariableUnits.size();
        if (index - offset < mParameterUnits.size())
        {
            return mParameterUnits[index - offset];
        }
        else
        {
            offset += mParameterUnits.size();
            if (index - offset < mDerivedQuantityUnits.size())
            {
                return mDerivedQuantityUnits[index - offset];
            }
            else
            {
                EXCEPTION("Invalid index passed to GetAnyVariableUnits.");
            }
        }
    }
}


const std::vector<std::string>& AbstractOdeSystemInformation::rGetDerivedQuantityNames() const
{
    assert(mInitialised);
    return mDerivedQuantityNames;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetDerivedQuantityUnits() const
{
    assert(mInitialised);
    return mDerivedQuantityUnits;
}

unsigned AbstractOdeSystemInformation::GetDerivedQuantityIndex(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mDerivedQuantityNames.begin(), mDerivedQuantityNames.end(), rName);
    if (it == mDerivedQuantityNames.end())
    {
        EXCEPTION("No derived quantity named '" + rName + "'.");
    }
    return (unsigned)(it - mDerivedQuantityNames.begin());
}

bool AbstractOdeSystemInformation::HasDerivedQuantity(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mDerivedQuantityNames.begin(), mDerivedQuantityNames.end(), rName);
    return (it != mDerivedQuantityNames.end());
}

std::string AbstractOdeSystemInformation::GetDerivedQuantityUnits(unsigned index) const
{
    assert(mInitialised);
    if (index >= mDerivedQuantityUnits.size())
    {
        EXCEPTION("The index passed in must be less than the number of derived quantities.");
    }
    return mDerivedQuantityUnits[index];
}

unsigned AbstractOdeSystemInformation::GetNumberOfDerivedQuantities() const
{
    assert(mInitialised);
    return mDerivedQuantityUnits.size();
}

unsigned AbstractOdeSystemInformation::GetNumberOfAttributes() const
{
    assert(mInitialised);
    return mAttributes.size();
}

bool AbstractOdeSystemInformation::HasAttribute(const std::string& rName) const
{
    assert(mInitialised);
    return (mAttributes.find(rName) != mAttributes.end());
}

double AbstractOdeSystemInformation::GetAttribute(const std::string& rName) const
{
    assert(mInitialised);
    std::map<std::string, double>::const_iterator it = mAttributes.find(rName);
    if (it == mAttributes.end())
    {
        EXCEPTION("No attribute '" + rName + "' found.");
    }
    return it->second;
}
