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

#include "AbstractLookupTableCollection.hpp"

#include <cmath>

AbstractLookupTableCollection::AbstractLookupTableCollection()
    : mDt(0.0)
{
}

std::vector<std::string> AbstractLookupTableCollection::GetKeyingVariableNames() const
{
    return mKeyingVariableNames;
}

unsigned AbstractLookupTableCollection::GetNumberOfTables(const std::string& rKeyingVariableName) const
{
    return mNumberOfTables[GetTableIndex(rKeyingVariableName)];
}

void AbstractLookupTableCollection::GetTableProperties(const std::string& rKeyingVariableName, double& rMin, double& rStep, double& rMax) const
{
    unsigned i = GetTableIndex(rKeyingVariableName);
    rMin = mTableMins[i];
    rStep = mTableSteps[i];
    rMax = mTableMaxs[i];
}

void AbstractLookupTableCollection::SetTableProperties(const std::string& rKeyingVariableName, double min, double step, double max)
{
    // Check inputs
    unsigned num_steps = (unsigned) ((max-min)/step+0.5);
    ///\todo remove magic number? (#1884)
    if (fabs(max - min - num_steps*step) > 1e-10)
    {
        EXCEPTION("Table step size does not divide range between table limits.");
    }
    // Set state
    unsigned i = GetTableIndex(rKeyingVariableName);
    if ((min != mTableMins[i]) || (step != mTableSteps[i]) || (max != mTableMaxs[i]))
    {
        mNeedsRegeneration[i] = true;
    }
    mTableMins[i] = min;
    mTableSteps[i] = step;
    mTableStepInverses[i] = 1/step;
    mTableMaxs[i] = max;
}

void AbstractLookupTableCollection::SetTimestep(double dt)
{
    if (mDt != dt)
    {
        mNeedsRegeneration.assign(mNeedsRegeneration.size(), true);
    }
    mDt = dt;
}

unsigned AbstractLookupTableCollection::GetTableIndex(const std::string& rKeyingVariableName) const
{
    unsigned i=0;
    for (; i<mKeyingVariableNames.size(); i++)
    {
        if (mKeyingVariableNames[i] == rKeyingVariableName)
        {
            break;
        }
    }
    if (i == mKeyingVariableNames.size())
    {
        EXCEPTION("Lookup table keying variable '" + rKeyingVariableName + "' does not exist.");
    }
    return i;
}

AbstractLookupTableCollection::~AbstractLookupTableCollection()
{
}

const char* AbstractLookupTableCollection::EventHandler::EventName[] =  {"GenTables"};
