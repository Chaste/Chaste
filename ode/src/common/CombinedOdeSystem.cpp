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


#include "CombinedOdeSystem.hpp"
#include "CombinedOdeSystemInformation.hpp"



CombinedOdeSystem::CombinedOdeSystem(std::vector<AbstractOdeSystem*> odeSystems)
    : AbstractOdeSystem(0) // will be set properly below
{
    mOdeSystems = odeSystems;
    for (unsigned i=0; i<odeSystems.size(); i++)
    {
        mNumberOfStateVariables += odeSystems[i]->GetNumberOfStateVariables();
    }
    mpSystemInfo = CombinedOdeSystemInformation::Instance(odeSystems);
    ResetToInitialConditions();

    // Set up working memory
    unsigned num_systems = odeSystems.size();
    mWorkingStateVars.resize(num_systems);
    mWorkingDerivs.resize(num_systems);
    unsigned offset = 0;
    for (unsigned i=0; i<num_systems; i++)
    {
        unsigned num_vars = odeSystems[i]->GetNumberOfStateVariables();
        mWorkingStateVars[i].resize(num_vars);
        mWorkingDerivs[i].resize(num_vars);
        mOffsets.push_back(offset);
        offset += num_vars;
    }
}


void CombinedOdeSystem::Configure(
        const std::map<unsigned, unsigned>& rVariableParameterMap,
        AbstractOdeSystem* pVariableOdeSystem,
        AbstractOdeSystem* pParameterOdeSystem)
{
    struct VariableParameterMap new_map;
    new_map.theMap = rVariableParameterMap;
    unsigned var_system_index = 0;
    while (var_system_index < mOdeSystems.size() && mOdeSystems[var_system_index] != pVariableOdeSystem)
    {
        ++var_system_index;
    }
    new_map.pVariableOdeSystemIndex = var_system_index;
    new_map.pParameterOdeSystem = pParameterOdeSystem;
    mVariableParameterMaps.push_back(new_map);
}


void CombinedOdeSystem::EvaluateYDerivatives(
        double time,
        const std::vector<double>& rY,
        std::vector<double>& rDY)
{
    // Copy rY to subsystems
    for (unsigned i=0; i<mOdeSystems.size(); i++)
    {
        unsigned offset = mOffsets[i];
        for (unsigned j=0; j<mOdeSystems[i]->GetNumberOfStateVariables(); j++)
        {
            mWorkingStateVars[i][j] = rY[offset + j];
        }
    }

    // Set parameter values
    for (unsigned i=0; i<mVariableParameterMaps.size(); i++)
    {
        std::map<unsigned, unsigned>& r_var_param_map = mVariableParameterMaps[i].theMap;
        // Iterate through map
        for (std::map<unsigned, unsigned>::iterator iter = r_var_param_map.begin();
             iter != r_var_param_map.end();
             ++iter)
        {
            double value = mWorkingStateVars[mVariableParameterMaps[i].pVariableOdeSystemIndex][iter->first];
            mVariableParameterMaps[i].pParameterOdeSystem->SetParameter(iter->second, value);
        }
    }

    // Call EvaluateYDerivatives on subsystems
    for (unsigned i=0; i<mOdeSystems.size(); i++)
    {
        mOdeSystems[i]->EvaluateYDerivatives(time, mWorkingStateVars[i], mWorkingDerivs[i]);
    }

    // Copy derivatives to rDY
    for (unsigned i=0; i<mOdeSystems.size(); i++)
    {
        unsigned offset = mOffsets[i];
        for (unsigned j=0; j<mOdeSystems[i]->GetNumberOfStateVariables(); j++)
        {
            rDY[offset + j] = mWorkingDerivs[i][j];
        }
    }
}
