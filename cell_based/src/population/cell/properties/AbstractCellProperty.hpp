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

#ifndef ABSTRACTCELLPROPERTY_HPP_
#define ABSTRACTCELLPROPERTY_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include "Identifiable.hpp"

/**
 * Base class for cell properties.
 *
 * Each cell has a collection of cell properties, which may express such concepts as
 * mutation states, inherited labels, or per-cell data.
 */
class AbstractCellProperty : public Identifiable
{
private:

    /**
     * The number of cells with this cell property.
     */
    unsigned mCellCount;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * We don't have any member variables yet, but let's make things easy for the future.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCellCount;
    }

public:
    /**
     * Default constructor.
     */
    AbstractCellProperty();

    /**
     * Virtual destructor, to make this class polymorphic.
     */
    virtual ~AbstractCellProperty();

    /**
     * Test whether this property is a particular class.  This tests for exact
     * run-time class identity, and doesn't match subclasses.
     *
     * The current implementation requires that all cell properties have a default
     * constructor.  For this method to be efficient, the default constructor must
     * also be cheap.
     *
     * It should be called like:
     *   bool healthy = p_property->IsType<HealthyMutationState>();
     */
    template<class CLASS>
    bool IsType() const
    {
        CLASS ref_obj;
        return IsSame(&ref_obj);
    }

    /**
     * Test whether this property is an instance of a particular class, or any subclass
     * of that class.
     *
     * It should be called like:
     *   bool is_mutation = p_property->IsSubType<AbstractCellMutationState>();
     */
    template<class BASECLASS>
    bool IsSubType() const
    {
        // We put a const_cast in here so users don't have to worry about whether the
        // property object they are testing is const or not.
        BASECLASS* p_subclass = dynamic_cast<BASECLASS*>(const_cast<AbstractCellProperty*>(this));
        return (p_subclass != NULL);
    }

    /**
     * Determine whether this property is the same class as another.
     * @param pOther  the property to compare against.
     */
    bool IsSame(const AbstractCellProperty* pOther) const;

    /**
     * Determine whether this property is the same as another.
     * @param pOther  the property to compare against.
     */
    bool IsSame(boost::shared_ptr<const AbstractCellProperty> pOther) const;

    /**
     * Increment #mCellCount.
     */
    void IncrementCellCount();

    /**
     * Decrement #mCellCount.
     */
    void DecrementCellCount();

    /**
     * Get #mCellCount
     */
    unsigned GetCellCount() const;

};

#endif /* ABSTRACTCELLPROPERTY_HPP_ */
