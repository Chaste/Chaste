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
     * Get #mColour.
     */
    unsigned GetColour() const;
};

#endif /* ABSTRACTCELLMUTATIONSTATE_HPP_ */
