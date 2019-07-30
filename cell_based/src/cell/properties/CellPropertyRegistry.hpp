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
     * The main interface to this class.
     * @return a particular cell property object.
     * Use like:
     *    boost::shared_ptr<AbstractCellProperty> p_property(
                CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
     */
    template<class SUBCLASS>
    boost::shared_ptr<AbstractCellProperty> Get();

    /**
     * @return the single instance of the registry.
     */
    static CellPropertyRegistry* Instance();

    /**
     * @return a list of the cell properties registered.
     */
    const std::vector<boost::shared_ptr<AbstractCellProperty> >& rGetAllCellProperties();

    /**
     * Clear all registered cell properties.
     */
    void Clear();

    /**
     * Take ownership of the current registry. Calling Instance after this will
     * create a new registry. The caller takes responsibility for freeing the
     * returned registry when finished with it.
     *
     * This method is intended for use by CellBasedSimulation, so that we can have
     * multiple concurrent simulations, each with their own registry.
     * @return registry
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
     * @return reference by convention
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
