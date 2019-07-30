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

#ifndef DYNAMICMODELLOADERREGISTRY_HPP_
#define DYNAMICMODELLOADERREGISTRY_HPP_

#include <memory>
#include <string>
#include <map>
#include <set>
#include <boost/utility.hpp>
#include <boost/weak_ptr.hpp>
#include "DynamicCellModelLoader.hpp"
#include "FileFinder.hpp"

typedef boost::weak_ptr<DynamicCellModelLoader> DynamicCellModelLoaderWeakPtr;

/**
 * This class provides a static registry to keep track of the cell model loaders used,
 * hence ensuring that we don't have a given .so loaded more than once at any given point.
 */
class DynamicModelLoaderRegistry : private boost::noncopyable
{
public:
    /**
     * @return the single instance of the registry.
     */
    static DynamicModelLoaderRegistry* Instance();

    /**
     * @return the loader for the given .so file.
     * @param rPath  absolute path to the .so
     */
    DynamicCellModelLoaderPtr GetLoader(const std::string& rPath);

    /**
     * @return the loader for the given .so file.
     * @param rFileFinder  finder for the .so file
     */
    DynamicCellModelLoaderPtr GetLoader(const FileFinder& rFileFinder);

    /**
     * Schedule the given loader for deletion prior to loading any new .so.
     *
     * @param pLoader  the loader to schedule for deletion
     */
    void ScheduleForDeletion(DynamicCellModelLoaderPtr pLoader);

private:
    /**
     * Loaders for shared-library cell models.
     * Weak pointers are used so that the registry doesn't keep loaders alive when all
     * cells created from them have been destroyed.
     * Map is from absolute path of the library, to loader object.
     */
    std::map<std::string, DynamicCellModelLoaderWeakPtr> mLoaders;

    /** Loaders to be deleted before creating any new ones. */
    std::set<DynamicCellModelLoaderPtr> mDeletableLoaders;

    /** The single instance of this class. */
    static std::shared_ptr<DynamicModelLoaderRegistry> mpInstance;

    /**
     * Private constructor; all access should be via Instance().
     */
    DynamicModelLoaderRegistry();
};

#endif /*DYNAMICMODELLOADERREGISTRY_HPP_*/
