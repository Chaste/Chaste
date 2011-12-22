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

#ifndef CELLPROPERTYCOLLECTION_HPP_
#define CELLPROPERTYCOLLECTION_HPP_

#include <set>
#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/set.hpp>

#include "AbstractCellProperty.hpp"
#include "CellPropertyRegistry.hpp"
#include "Exception.hpp"

/**
 * Cell property collection class.
 *
 * Contains methods for accessing and interrogating a set of cell properties.
 */
class CellPropertyCollection
{
private:
    /** The type of container used to store properties */
    typedef std::set<boost::shared_ptr<AbstractCellProperty> > CollectionType;

    /** Type of a const iterator over the container */
    typedef CollectionType::const_iterator ConstIteratorType;

    /** Type of an iterator over the container */
    typedef CollectionType::iterator IteratorType;

    /** The properties stored in this collection. */
    CollectionType mProperties;

    /** Cell property registry. */
    CellPropertyRegistry* mpCellPropertyRegistry;

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
        archive & mProperties;
        // archive & mpCellPropertyRegistry; Not required as archived by the CellPopulation.
    }

public:
    /**
     * Create an empty collection of cell properties.
     */
    CellPropertyCollection();

    /**
     * Add a new property to this collection.
     *
     * @param rProp  the property to add
     */
    void AddProperty(const boost::shared_ptr<AbstractCellProperty>& rProp);

    /**
     * @return a pointer to the CellPropertyRegistry (as assigned to this cell in the AbstractCellPopulation constructor).
     */
    CellPropertyRegistry* GetCellPropertyRegistry();

    /**
     * Set the CellPropertyRegistry for this cell
     *
     * @param pRegistry  The cell property registry (assigned in the AbstractCellPopulation constructor).
     */
    void SetCellPropertyRegistry(CellPropertyRegistry* pRegistry);

    /**
     * Test whether this collection contains the given property @b object.
     *
     * @param rProp  the property to compare against
     */
    bool HasProperty(const boost::shared_ptr<AbstractCellProperty>& rProp) const;

    /**
     * Test whether the collection contains a property that has the exact type CLASS.
     *
     * Should be used like
     *   bool healthy = collection.HasProperty<WildTypeCellMutationState>();
     */
    template<typename CLASS>
    bool HasProperty() const
    {
        for (ConstIteratorType it = mProperties.begin(); it != mProperties.end(); ++it)
        {
            if ((*it)->IsType<CLASS>())
            {
                return true;
            }
        }
        return false;
    }

    /**
     * Test whether the collection contains a property that inherits from BASECLASS.
     *
     * Should be used like
     *   collection.HasPropertyType<AbstractCellMutationState>();
     */
    template<typename BASECLASS>
    bool HasPropertyType() const
    {
        for (ConstIteratorType it = mProperties.begin(); it != mProperties.end(); ++it)
        {
            if ((*it)->IsSubType<BASECLASS>())
            {
                return true;
            }
        }
        return false;
    }

    /**
     * Remove a property of the given type.
     */
    template<typename CLASS>
    void RemoveProperty()
    {
        for (IteratorType it = mProperties.begin(); it != mProperties.end(); ++it)
        {
            if ((*it)->IsType<CLASS>())
            {
                mProperties.erase(it);
                return;
            }
        }
        EXCEPTION("Collection does not contain the given property type.");
    }

    /**
     * Remove the given property from this collection.
     *
     * @param rProp  the property to remove
     */
    void RemoveProperty(const boost::shared_ptr<AbstractCellProperty>& rProp);

    /**
     * Get the size of this container.
     */
    unsigned GetSize() const;

    /**
     * An iterator type over this collection.
     * Don't rely on the particular implementation of the iterator.
     */
    typedef CollectionType::iterator Iterator;

    /**
     * Get an Iterator to the start of this collection.
     */
    Iterator Begin();

    /**
     * Get an Iterator to one past the end of this collection.
     */
    Iterator End();

    /**
     * If this collection contains a single property, then return it.
     * Otherwise, throws an exception.
     */
    boost::shared_ptr<AbstractCellProperty> GetProperty() const;

    /**
     * Get a sub-collection containing all our properties that are instances
     * of the given class.
     */
    template<typename CLASS>
    CellPropertyCollection GetProperties() const
    {
        CellPropertyCollection result;
        for (ConstIteratorType it = mProperties.begin(); it != mProperties.end(); ++it)
        {
            if ((*it)->IsType<CLASS>())
            {
                result.AddProperty(*it);
            }
        }
        return result;
    }

    /**
     * Get a sub-collection containing all our properties that are instances
     * of the given class or any of its subclasses.
     */
    template<typename BASECLASS>
    CellPropertyCollection GetPropertiesType() const
    {
        CellPropertyCollection result;
        for (ConstIteratorType it = mProperties.begin(); it != mProperties.end(); ++it)
        {
            if ((*it)->IsSubType<BASECLASS>())
            {
                result.AddProperty(*it);
            }
        }
        return result;
    }
};

#endif /* CELLPROPERTYCOLLECTION_HPP_ */
