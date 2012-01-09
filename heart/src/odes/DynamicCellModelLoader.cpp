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

#include "DynamicCellModelLoader.hpp"
#include "AbstractDynamicallyLoadableEntity.hpp"

#include <dlfcn.h>  // For dealing with .so files

#include "Exception.hpp"

DynamicCellModelLoader::DynamicCellModelLoader(const std::string& rLoadableModulePath)
    : mLoadableModulePath(rLoadableModulePath)
{
    // Open the .so file
    mpDynamicModule = dlopen(rLoadableModulePath.c_str(), RTLD_NOW);
    if (mpDynamicModule == NULL)
    {
        EXCEPTION("Unable to load .so file '" + rLoadableModulePath + "': " + dlerror());
    }

    // Find the cell creation function
    dlerror(); // Reset errors
    void* p_creation_function = dlsym(mpDynamicModule, "MakeCardiacCell");
    // Check it exists
    const char* p_error = dlerror();
    if (p_error)
    {
        EXCEPTION("Failed to load cell creation function from .so file '" + rLoadableModulePath
                  + "': " + p_error);
    }
    // Cast to the right type
    mpCreationFunction = reinterpret_cast<CellCreationFunctionType*>(p_creation_function);
}

DynamicCellModelLoader::~DynamicCellModelLoader()
{
    if (mpDynamicModule)
    {
        dlclose(mpDynamicModule);
    }
}

AbstractCardiacCellInterface* DynamicCellModelLoader::CreateCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                                                        boost::shared_ptr<AbstractStimulusFunction> pStimulus)
{
    AbstractCardiacCellInterface* p_cell = (*mpCreationFunction)(pSolver, pStimulus);

    AbstractDynamicallyLoadableEntity* p_entity = dynamic_cast<AbstractDynamicallyLoadableEntity*>(p_cell);

    if (p_entity != NULL)
    {
        p_entity->SetLoader(this);
    }
    return p_cell;
}

const std::string DynamicCellModelLoader::GetLoadableModulePath() const
{
    return mLoadableModulePath;
}
