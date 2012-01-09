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

#include "AbstractCardiacCellWithModifiers.hpp"
#include "Modifiers.hpp"

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
#endif


