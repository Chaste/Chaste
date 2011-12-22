/*

Copyright (C) University of Oxford, 2005-2011

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
