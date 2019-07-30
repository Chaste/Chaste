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
