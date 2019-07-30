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
     * Save/load our member variables.
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
     * @return whether this collection contains the given property @b object.
     *
     * @param rProp  the property to compare against
     */
    bool HasProperty(const boost::shared_ptr<AbstractCellProperty>& rProp) const;

    /**
     * @return whether the collection contains a property that has the exact type CLASS.
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
     * @return whether the collection contains a property that inherits from BASECLASS.
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
     * Remove a single property of the given type.
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
     * @return the size of this container.
     */
    unsigned GetSize() const;

    /**
     * An iterator type over this collection.
     * Don't rely on the particular implementation of the iterator.
     */
    typedef CollectionType::iterator Iterator;

    /**
     * @return an Iterator to the start of this collection.
     */
    Iterator Begin();

    /**
     * @return an Iterator to one past the end of this collection.
     */
    Iterator End();

    /**
     * @return If this collection contains a single property, then return it.
     * Otherwise, throws an exception.
     */
    boost::shared_ptr<AbstractCellProperty> GetProperty() const;

    /**
     * @return a sub-collection containing all our properties that are instances
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
     * @return a sub-collection containing all our properties that are instances
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
