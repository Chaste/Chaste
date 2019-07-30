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

#ifndef ABSTRACTCARDIACCELLWITHMODIFIERS_HPP_
#define ABSTRACTCARDIACCELLWITHMODIFIERS_HPP_

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>

#include "AbstractCvodeCell.hpp"
#include "AbstractCvodeCellWithDataClamp.hpp"
#include "AbstractCardiacCell.hpp"

#include "AbstractModifier.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

/**
 * A base class for cardiac cells that have been altered to include calls to subclasses
 * of AbstractModifier when computing their derivatives.
 */
template<class CARDIAC_CELL>
class AbstractCardiacCellWithModifiers : public CARDIAC_CELL
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<CARDIAC_CELL>(*this);

        // This is a bit unusual - as it contains pointers to member variables in the subclasses.
        // So we deal with it specially on reloading below.
        // The map is always set up in the same way, so should be in the same order, but we'll
        // archive the names as well for testing on load to be on the safe side.
        std::map<std::string, boost::shared_ptr<AbstractModifier>* >::const_iterator iter;
        for (iter = mModifiersMap.begin(); iter!= mModifiersMap.end(); ++iter)
        {
            const boost::shared_ptr<AbstractModifier>* p_to_smart_pointer = (*iter).second;
            const boost::shared_ptr<AbstractModifier> p_modifier = *p_to_smart_pointer;
            archive & (*iter).first; // Name of the modifier
            archive & p_modifier;    // Modifier
        }
    }
    /**
     * Unarchive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<CARDIAC_CELL>(*this);

        // We have made new smart pointers to dummy modifiers via the constructor, so instead of overwriting
        // these pointers with new pointers (leaving a memory leak), go through and make the existing pointers
        // point to the correct modifiers again!
        std::map<std::string, boost::shared_ptr<AbstractModifier>* >::iterator iter;
        for (iter = mModifiersMap.begin(); iter!= mModifiersMap.end(); ++iter)
        {
            boost::shared_ptr<AbstractModifier>* p_to_constructed_smart_pointer = (*iter).second;
            boost::shared_ptr<AbstractModifier> p_loaded;
            std::string modifier_name;
            archive & modifier_name;  // Name of the modifier
            archive & p_loaded;       // Modifier

            // Paranoia check that this is the modifier we think it is.
            if ((*iter).first != modifier_name)
            {
                NEVER_REACHED;
                // You're in trouble.
                // If this breaks, perhaps change this loop to do the right number and
                // use modifier_name as mModifiersMap key instead.
            }

            // Set the constructed smart pointer to be the same as the loaded one.
            *(p_to_constructed_smart_pointer) = p_loaded;
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /** A map between a string description and the location of the relevant modifier in concrete classes. */
    std::map<std::string, boost::shared_ptr<AbstractModifier>* > mModifiersMap;

protected:
    /**
     * Add a new modifier - should only be called by the subclass constructors.
     * Each modifier pointer is set to a #DummyModifier by this method.
     * @param modifierName  The name which will act as a 'key' for this modifier.
     * @param pModifier  The pointer to the modifier in the concrete class.
     */
    void AddModifier(std::string modifierName, boost::shared_ptr<AbstractModifier>& pModifier);

public:
    /**
     * Create a new cardiac cell.
     *
     * This calls the main CARDIAC_CELL constructor, but also supplies modifiers for
     * working with metadata-enabled CellML files.
     *
     * @param pOdeSolver  the ODE solver to use when simulating this cell
     * @param numberOfStateVariables  the size of the ODE system modelling this cell
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    AbstractCardiacCellWithModifiers(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                     unsigned numberOfStateVariables,
                                     unsigned voltageIndex,
                                     boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Get access to a modifier
     *
     * @param rModifierName  The oxmeta name of the modifier to fetch.
     * @return a pointer to the specified modifier
     */
    boost::shared_ptr<AbstractModifier> GetModifier(const std::string& rModifierName);

    /**
     * Check for the presence of a modifier. To avoid the necessity of using
     * try...catch around every call to GetModifier().
     *
     * @param rModifierName  The oxmeta name of a modifier
     * @return whether or not the cell has one of these modifiers.
     */
    bool HasModifier(const std::string& rModifierName) const;

    /**
     * Set a new modifier
     *
     * @param rModifierName  The oxmeta name of the modifier to replace.
     * @param pNewModifier  The new modifier object to use.
     */
    void SetModifier(const std::string& rModifierName, boost::shared_ptr<AbstractModifier>& pNewModifier);
};

// Special case of archiving an abstract class that's templated over Chaste classes.
// See the comments in the top of global/src/checkpointing/ClassIsAbstract.hpp
namespace boost {
namespace serialization {

    template<class C>
    struct is_abstract<AbstractCardiacCellWithModifiers<C> >
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN

    template<class C>
    struct is_abstract<const AbstractCardiacCellWithModifiers<C> >
        TEMPLATED_CLASS_IS_ABSTRACT_DEFN
}}


#endif // ABSTRACTCARDIACCELLWITHMODIFIERS_HPP_
