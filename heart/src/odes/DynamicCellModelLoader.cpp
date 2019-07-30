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

#include "DynamicCellModelLoader.hpp"
#include "AbstractDynamicallyLoadableEntity.hpp"

#include <dlfcn.h>  // For dealing with .so files
#include <cassert>

#include "Exception.hpp"

DynamicCellModelLoaderPtr DynamicCellModelLoader::Create(const std::string& rLoadableModulePath)
{
    DynamicCellModelLoaderPtr p_loader(new DynamicCellModelLoader(rLoadableModulePath));
    return p_loader;
}


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
    assert(p_entity != NULL);
    p_entity->SetLoader(shared_from_this());

    return p_cell;
}

const std::string DynamicCellModelLoader::GetLoadableModulePath() const
{
    return mLoadableModulePath;
}
