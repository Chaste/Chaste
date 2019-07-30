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

#include "AbstractCardiacCellWithModifiers.hpp"
#include "AbstractCvodeCellWithDataClamp.hpp"
#include "DummyModifier.hpp"

template<class CARDIAC_CELL>
void AbstractCardiacCellWithModifiers<CARDIAC_CELL>::AddModifier(std::string modifierName, boost::shared_ptr<AbstractModifier>& pModifier)
{
    mModifiersMap[modifierName] = &pModifier;
    pModifier = boost::shared_ptr<AbstractModifier>(new DummyModifier()); // This modifier always returns what is passed in.
}

template<class CARDIAC_CELL>
AbstractCardiacCellWithModifiers<CARDIAC_CELL>::AbstractCardiacCellWithModifiers(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                 unsigned numberOfStateVariables,
                                 unsigned voltageIndex,
                                 boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : CARDIAC_CELL(pOdeSolver, numberOfStateVariables, voltageIndex, pIntracellularStimulus)
{
    mModifiersMap.clear();
}

template<class CARDIAC_CELL>
boost::shared_ptr<AbstractModifier> AbstractCardiacCellWithModifiers<CARDIAC_CELL>::GetModifier(const std::string& rModifierName)
{
    if (mModifiersMap.find(rModifierName) == mModifiersMap.end())
    {
        EXCEPTION("There is no modifier called " + rModifierName + " in this model.");
    }
    return *(mModifiersMap[rModifierName]);
}

template<class CARDIAC_CELL>
bool AbstractCardiacCellWithModifiers<CARDIAC_CELL>::HasModifier(const std::string& rModifierName) const
{
    return !(mModifiersMap.find(rModifierName) == mModifiersMap.end());
}

template<class CARDIAC_CELL>
void AbstractCardiacCellWithModifiers<CARDIAC_CELL>::SetModifier(const std::string& rModifierName, boost::shared_ptr<AbstractModifier>& pNewModifier)
{
    if (mModifiersMap.find(rModifierName) == mModifiersMap.end())
    {
        EXCEPTION("There is no modifier called " + rModifierName + " in this model.");
    }
    *(mModifiersMap[rModifierName]) = pNewModifier;
}

// Explicit Instantiation
template class AbstractCardiacCellWithModifiers<AbstractCardiacCell>;
#ifdef CHASTE_CVODE
template class AbstractCardiacCellWithModifiers<AbstractCvodeCell>;
template class AbstractCardiacCellWithModifiers<AbstractCvodeCellWithDataClamp>;
#endif


