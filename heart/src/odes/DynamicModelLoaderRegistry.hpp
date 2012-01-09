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

#ifndef DYNAMICMODELLOADERREGISTRY_HPP_
#define DYNAMICMODELLOADERREGISTRY_HPP_

#include <memory>
#include <string>
#include <map>
#include "DynamicCellModelLoader.hpp"
#include "FileFinder.hpp"

/**
 * When loading cell models dynamically, the loader object needs to be alive for as long
 * as cells created by it are alive.  Unfortunately, the HeartConfigRelatedCellFactory
 * is typically destroyed as soon as the cells have been created, prior to the simulation.
 * Hence, this class provides a static registry to keep track of the model loaders used.
 * It also ensures we don't load a given .so more than once.
 */
class DynamicModelLoaderRegistry
{
public:
    /**
     * Get the single instance of the registry.
     */
    static DynamicModelLoaderRegistry* Instance();

    /**
     * Get the loader for the given .so file.
     * @param rPath  absolute path to the .so
     */
    DynamicCellModelLoader* GetLoader(const std::string& rPath);

    /**
     * Get the loader for the given .so file.
     * @param rFileFinder  finder for the .so file
     */
    DynamicCellModelLoader* GetLoader(const FileFinder& rFileFinder);

    /**
     * Destructor closes all loaded .so files.
     */
    ~DynamicModelLoaderRegistry();

private:
    /**
     * Loaders for shared-library cell models.
     * Map is from absolute path of the library, to loader object.
     */
    std::map<std::string, DynamicCellModelLoader*> mLoaders;

    /** The single instance of this class. */
    static std::auto_ptr<DynamicModelLoaderRegistry> mpInstance;

    /**
     * Private constructor; all access should be via Instance().
     */
    DynamicModelLoaderRegistry();

    /**
     * Copy constructor.
     */
    DynamicModelLoaderRegistry(const DynamicModelLoaderRegistry&);

    /**
     * Overloaded assignment operator.
     */
    DynamicModelLoaderRegistry& operator= (const DynamicModelLoaderRegistry&);

};

#endif /*DYNAMICMODELLOADERREGISTRY_HPP_*/
