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

#ifndef CELLPROPERTYREGISTRY_HPP_
#define CELLPROPERTYREGISTRY_HPP_

#include <boost/shared_ptr.hpp>
#include <vector>

#include "AbstractCellProperty.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A singleton registry of available cell properties.
 */
class CellPropertyRegistry
{
public:
    /**
     * The main interface to this class: get a particular cell property object.
     * Use like:
     *    boost::shared_ptr<AbstractCellProperty> p_property(
                CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
     */
    template<class SUBCLASS>
    boost::shared_ptr<AbstractCellProperty> Get();

    /**
     * Get the single instance of the registry.
     */
    static CellPropertyRegistry* Instance();

    /**
     * Get a list of the cell properties registered.
     */
    const std::vector<boost::shared_ptr<AbstractCellProperty> >& rGetAllCellProperties();

    /**
     * Clear all registered cell proeprties.
     */
    void Clear();

    /**
     * Take ownership of the current registry. Calling Instance after this will
     * create a new registry. The caller takes responsibility for freeing the
     * returned registry when finished with it.
     *
     * This method is intended for use by CellBasedSimulation, so that we can have
     * multiple concurrent simulations, each with their own registry.
     */
    CellPropertyRegistry* TakeOwnership();

    /**
     * Specify the ordering in which cell properties should be returned by rGetAllCellProperties().
     * The provided ordering must include all cell properties in the registry. Once an ordering
     * has been specified, attempts to get a cell property which is not in the ordering will
     * throw an exception.
     *
     * @param rOrdering  vector of cell properties in the desired order
     */
    void SpecifyOrdering(const std::vector<boost::shared_ptr<AbstractCellProperty> >& rOrdering);

    /**
     * @return whether an ordering has been specified.
     */
    bool HasOrderingBeenSpecified();

private:

    /**
     * Default constructor.
     */
    CellPropertyRegistry();

    /**
     * Copy constructor.
     */
    CellPropertyRegistry(const CellPropertyRegistry&);

    /**
     * Overloaded assignment operator.
     */
    CellPropertyRegistry& operator= (const CellPropertyRegistry&);

    /**
     * A pointer to the singleton instance of this class.
     */
    static CellPropertyRegistry* mpInstance;

    /**
     * The cell properties in the registry.
     */
    std::vector<boost::shared_ptr<AbstractCellProperty> > mCellProperties;

    /** Whether an ordering has been set up */
    bool mOrderingHasBeenSpecified;

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
        archive & mCellProperties;
        archive & mOrderingHasBeenSpecified;
    }
};

template<class SUBCLASS>
boost::shared_ptr<AbstractCellProperty> CellPropertyRegistry::Get()
{
    boost::shared_ptr<AbstractCellProperty> p_property;
    for (unsigned i=0; i<mCellProperties.size(); i++)
    {
        if (mCellProperties[i]->IsType<SUBCLASS>())
        {
            p_property = mCellProperties[i];
            break;
        }
    }
    if (!p_property)
    {
        // Create a new cell property
        p_property.reset(new SUBCLASS);
        mCellProperties.push_back(p_property);
    }
    return p_property;
}

#endif /* CELLPROPERTYREGISTRY_HPP_ */
