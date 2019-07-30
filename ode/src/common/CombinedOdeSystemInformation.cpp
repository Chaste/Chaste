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


#include "CombinedOdeSystemInformation.hpp"

#include <cassert>

boost::shared_ptr<CombinedOdeSystemInformation> CombinedOdeSystemInformation::Instance(const std::vector<AbstractOdeSystem*>& rSubsystems)
{
    // Get the information for the subsystems
    std::vector<boost::shared_ptr<const AbstractOdeSystemInformation> > info_vec;
    info_vec.reserve(rSubsystems.size());
    for (unsigned i=0; i<rSubsystems.size(); i++)
    {
        info_vec.push_back(rSubsystems[i]->GetSystemInformation());
    }

    boost::shared_ptr<CombinedOdeSystemInformation> p_inst;

    // Search to see if we have an information object for this sequence of
    // subsystems already.
    for (unsigned i=0; i<msInstances.size(); i++)
    {
        if (info_vec.size() == msInstances[i].subsystemInformation.size())
        {
            bool equal = true;
            for (unsigned j=0; j<info_vec.size(); j++)
            {
                if (msInstances[i].subsystemInformation[j] != info_vec[j])
                {
                    equal = false;
                    break;
                }
            }
            if (equal)
            {
                p_inst = msInstances[i].pInfoInstance;
                break;
            }
        }
    }

    // Create a new object if needed
    if (!p_inst)
    {
        p_inst.reset(new CombinedOdeSystemInformation(info_vec));
        struct InstancePointers inst;
        inst.subsystemInformation = info_vec;
        inst.pInfoInstance = p_inst;
        msInstances.push_back(inst);
    }

    return p_inst;
}

CombinedOdeSystemInformation::CombinedOdeSystemInformation(const std::vector<boost::shared_ptr<const AbstractOdeSystemInformation> >& rSubsystemInfo)
{
    // Figure out our size
    unsigned total_system_size = 0u;
    for (unsigned i=0; i<rSubsystemInfo.size(); i++)
    {
        total_system_size += rSubsystemInfo[i]->rGetStateVariableNames().size();
    }
    mVariableNames.reserve(total_system_size);
    mVariableUnits.reserve(total_system_size);
    mInitialConditions.reserve(total_system_size);

    // Set up our info from the subsystems
    for (unsigned i=0; i<rSubsystemInfo.size(); i++)
    {
        std::vector<double> inits = rSubsystemInfo[i]->GetInitialConditions();
        const std::vector<std::string>& names = rSubsystemInfo[i]->rGetStateVariableNames();
        const std::vector<std::string>& units = rSubsystemInfo[i]->rGetStateVariableUnits();
        unsigned system_size = names.size();
        assert(inits.size() == system_size);
        assert(units.size() == system_size);

        for (unsigned j=0; j<system_size; j++)
        {
            mVariableNames.push_back(names[j]);
            mVariableUnits.push_back(units[j]);
            mInitialConditions.push_back(inits[j]);
        }
    }

    mInitialised = true;
}

// LCOV_EXCL_START
void CombinedOdeSystemInformation::Initialise()
{
    // does nothing; work done in constructor
    // but we need the method because it is pure in our base class
}

/** Definition of the instance static member. */
std::vector<struct CombinedOdeSystemInformation::InstancePointers> CombinedOdeSystemInformation::msInstances;
// LCOV_EXCL_STOP
