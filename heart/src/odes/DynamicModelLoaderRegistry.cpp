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

#include "DynamicModelLoaderRegistry.hpp"


DynamicModelLoaderRegistry* DynamicModelLoaderRegistry::Instance()
{
    if (!mpInstance.get())
    {
        mpInstance.reset(new DynamicModelLoaderRegistry());
    }
    return mpInstance.get();
}

DynamicCellModelLoader* DynamicModelLoaderRegistry::GetLoader(const std::string& rPath)
{
    // Have we opened this library already?
    std::map<std::string, DynamicCellModelLoader*>::iterator it = mLoaders.find(rPath);
    if (it == mLoaders.end())
    {
        // No
        mLoaders[rPath] = new DynamicCellModelLoader(rPath);
    }
    return mLoaders[rPath];
}

DynamicCellModelLoader* DynamicModelLoaderRegistry::GetLoader(const FileFinder& rFileFinder)
{
    return GetLoader(rFileFinder.GetAbsolutePath());
}

DynamicModelLoaderRegistry::~DynamicModelLoaderRegistry()
{
    for (std::map<std::string, DynamicCellModelLoader*>::iterator it=mLoaders.begin();
         it != mLoaders.end();
         ++it)
    {
        delete (it->second);
    }
    mLoaders.clear();
}

DynamicModelLoaderRegistry::DynamicModelLoaderRegistry()
{
}

/** The loader registry instance */
std::auto_ptr<DynamicModelLoaderRegistry> DynamicModelLoaderRegistry::mpInstance;
