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

#ifndef DYNAMICCELLMODELLOADER_HPP_
#define DYNAMICCELLMODELLOADER_HPP_

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "AbstractCardiacCellInterface.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"

// Forward reference
class DynamicCellModelLoader;

/** The main type for dealing with a loader instance. */
typedef boost::shared_ptr<DynamicCellModelLoader> DynamicCellModelLoaderPtr;

/**
 * This class takes care of loading cell models at run-time from .so files.
 *
 * Instantiate it using the factory function Create with the path of a suitable .so file,
 * then call CreateCell to create individual cells.
 */
class DynamicCellModelLoader : public boost::enable_shared_from_this<DynamicCellModelLoader>
{
public:
    /**
     * Create a cell model loader by opening a loadable module (.so file) containing
     * a cell model.
     *
     * @note The loader object must remain alive for as long as you want to use cells
     * created with it, or you'll get a segfault.  Hence the constructor is private,
     * and this factory function must be used to create a loader.
     *
     * @param rLoadableModulePath  path to .so file
     * @return loader object
     */
    static DynamicCellModelLoaderPtr Create(const std::string& rLoadableModulePath);

    /**
     * Destructor.  Closes the .so file.
     */
    ~DynamicCellModelLoader();

    /**
     * @return a newly created cardiac cell from this dynamic module.
     *
     * The caller takes responsibility for deleting the cell when it's finished with.
     *
     * @param pSolver  ODE solver used to simulate the cell
     * @param pStimulus  intracellular stimulus
     */
    AbstractCardiacCellInterface* CreateCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                                             boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /**
     * @return the absolute path to the .so file we have loaded
     */
    const std::string GetLoadableModulePath() const;

private:
    /**
     * Private constructor to ensure we're always stored in a shared pointer.
     *
     * @param rLoadableModulePath  path to .so file
     */
    DynamicCellModelLoader(const std::string& rLoadableModulePath);

    /** Handle for the loaded .so file */
    void* mpDynamicModule;

    /** Type of the cell creation function in the .so files
     *
     * @param pSolver  ODE solver used to simulate the cell
     * @param pStimulus  intracellular stimulus
     */
    typedef AbstractCardiacCellInterface* CellCreationFunctionType(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                                                                   boost::shared_ptr<AbstractStimulusFunction> pStimulus);

    /** Our cell creation function */
    CellCreationFunctionType* mpCreationFunction;

    /** Absolute path to the .so file we have loaded. */
    std::string mLoadableModulePath;
};

#endif /* DYNAMICCELLMODELLOADER_HPP_ */
