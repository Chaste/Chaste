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

#define COVERAGE_IGNORE
void CombinedOdeSystemInformation::Initialise()
{
    // does nothing; work done in constructor
    // but we need the method because it is pure in our base class
}

/** Definition of the instance static member. */
std::vector<struct CombinedOdeSystemInformation::InstancePointers> CombinedOdeSystemInformation::msInstances;
#undef COVERAGE_IGNORE
