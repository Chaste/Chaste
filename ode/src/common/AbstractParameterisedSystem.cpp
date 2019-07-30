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

#include <sstream>
#include <cassert>

#include "AbstractParameterisedSystem.hpp"

#include "Exception.hpp"
#include "VectorHelperFunctions.hpp"


template<typename VECTOR>
AbstractParameterisedSystem<VECTOR>::AbstractParameterisedSystem(unsigned numberOfStateVariables)
    : AbstractUntemplatedParameterisedSystem(numberOfStateVariables)
{
    InitialiseEmptyVector(mParameters);
    InitialiseEmptyVector(mStateVariables);
}

template<typename VECTOR>
std::string AbstractParameterisedSystem<VECTOR>::DumpState(const std::string& rMessage)
{
    return GetStateMessage(rMessage, mStateVariables);
}

template<typename VECTOR>
std::string AbstractParameterisedSystem<VECTOR>::DumpState(const std::string& rMessage,
                                                           VECTOR Y)
{
    return GetStateMessage(rMessage, Y);
}

template<typename VECTOR>
std::string AbstractParameterisedSystem<VECTOR>::DumpState(const std::string& rMessage,
                                                           VECTOR Y,
                                                           double time)
{
    std::stringstream extra_message;
    extra_message << std::endl << "At independent variable (usually time) = " << time;
    std::string new_message = rMessage + extra_message.str();
    return GetStateMessage(new_message, Y);
}

template<typename VECTOR>
std::string AbstractParameterisedSystem<VECTOR>::GetStateMessage(const std::string& rMessage, VECTOR Y)
{
    std::stringstream res;
    res << rMessage << std::endl << "State:" << std::endl;
    assert(rGetStateVariableNames().size()==GetVectorSize(Y));
    const std::vector<std::string>& r_units = rGetStateVariableUnits();
    for (unsigned i=0; i<GetVectorSize(Y); i++)
    {
        res << "\t" << rGetStateVariableNames()[i] << ":" << GetVectorComponent(Y, i);
        if (!r_units.empty())
        {
            res << " " << r_units[i];
        }
        res << std::endl;
    }
    return res.str();
}

template<typename VECTOR>
void AbstractParameterisedSystem<VECTOR>::CheckParametersOnLoad(const std::vector<double>& rParameters, const std::vector<std::string>& rParameterNames)
{
    if (GetVectorSize(mParameters) != rGetParameterNames().size())
    {
        // Subclass constructor didn't give default values, so we need the archive to provide them all
        if (rParameterNames.size() != rGetParameterNames().size())
        {
            EXCEPTION("Number of ODE parameters in archive does not match number in class.");
        }
        CreateVectorIfEmpty(mParameters,rGetParameterNames().size());
    }

    // Check whether the archive specifies parameters that don't appear in this class,
    // and create a map from archive index to local index
    std::vector<unsigned> index_map(rParameterNames.size());
    for (unsigned i=0; i<rParameterNames.size(); ++i)
    {
        index_map[i] = find(rGetParameterNames().begin(), rGetParameterNames().end(), rParameterNames[i])
                       - rGetParameterNames().begin();
        if (index_map[i] == rGetParameterNames().size())
        {
            EXCEPTION("Archive specifies a parameter '" + rParameterNames[i] + "' which does not appear in this class.");
        }
    }

    for (unsigned i=0; i<rParameterNames.size(); ++i)
    {
        SetVectorComponent(mParameters,index_map[i],rParameters[i]);
    }

    // Paranoia check
    assert(GetVectorSize(mParameters) == rGetParameterNames().size());
}

//
// State variable methods
//

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
    if (mNumberOfStateVariables != GetVectorSize(rStateVariables))
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
    if (mNumberOfStateVariables <= index)
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
            if (pDerivedQuantities == nullptr)
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

//
// "Derived quantities" methods
//

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


//////////////// Explicit instantiation//////////////

template class AbstractParameterisedSystem<std::vector<double> >;
#ifdef CHASTE_CVODE
template class AbstractParameterisedSystem<N_Vector>;
#endif
