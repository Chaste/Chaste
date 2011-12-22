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

#ifndef ABSTRACTDYNAMICALLYLOADABLEENTITY_HPP_
#define ABSTRACTDYNAMICALLYLOADABLEENTITY_HPP_

#include "DynamicCellModelLoader.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>

/**
 * A mixin class for things that get loaded dynamically to maintain an association between instance objects and the shared library
 * they have been loaded from
 */
class AbstractDynamicallyLoadableEntity
{
private:

    /** The loader for our shared object file */
    DynamicCellModelLoader* mpLoader;

    friend class boost::serialization::access;
    /**
     * Save the path to the loadable module.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        const std::string so_path = GetLoader()->GetLoadableModulePath();
        archive & so_path;
    }
    /**
     * Load the path to the loadable module, and set our loader from the registry.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        std::string so_path;
        archive & so_path;
        SetLoader(DynamicModelLoaderRegistry::Instance()->GetLoader(so_path));
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    /** Virtual destructor to ensure we're polymorphic */
    virtual ~AbstractDynamicallyLoadableEntity();

    /**
     * @return a shared pointer to the loader
     */
    const DynamicCellModelLoader* GetLoader() const;

    /**
     * Should only be called by a dynamic cell model loader
     *
     * @param pLoader a shared pointer to the loader
     */
    void SetLoader(DynamicCellModelLoader* pLoader);
};

#endif /*ABSTRACTDYNAMICALLYLOADABLEENTITY_HPP_*/


