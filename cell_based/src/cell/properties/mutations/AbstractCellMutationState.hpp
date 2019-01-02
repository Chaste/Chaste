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

#ifndef ABSTRACTCELLMUTATIONSTATE_HPP_
#define ABSTRACTCELLMUTATIONSTATE_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Base class for cell mutation states.
 *
 * Each Cell has a (shared pointer to a) mutation state instance, which will
 * be an instance of a subclass of this class.  When setting up a CellBasedSimulation,
 * the user must specify a list of AbstractCellMutationState instances, which represent
 * the possible mutations that can occur in the simulation (including WildTypeCellMutationState).
 * This provides a registry of available mutation states, and cells will point to
 * one of these objects.
 *
 * The mutation state objects keep track of the number of cells in a given mutation
 * state, as well as what colour should be used by the visualizer to display cells
 * in each state.
 */
class AbstractCellMutationState : public AbstractCellProperty
{
private:

    /**
     * Colour for use by visualizer.
     */
    unsigned mColour;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

    /**
     * Default constructor needs to be defined for archiving, but never actually used,
     * since subclasses call the normal constructor.
     */
    AbstractCellMutationState();

public:

    /**
     * Constructor.
     *
     * @param colour  what colour cells with this mutation state should be in the visualizer
     */
    AbstractCellMutationState(unsigned colour);

    /**
     * Virtual destructor, to make this class polymorphic.
     */
    virtual ~AbstractCellMutationState();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

#endif /* ABSTRACTCELLMUTATIONSTATE_HPP_ */
