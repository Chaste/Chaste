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

#ifndef DYNAMICCELLMODELLOADER_HPP_
#define DYNAMICCELLMODELLOADER_HPP_

#include <string>
#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCellInterface.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractStimulusFunction.hpp"

/**
 * This class takes care of loading cell models at run-time from .so files.
 *
 * Instantiate it with the path of a suitable .so file, then call CreateCell
 * to create individual cells.
 */
class DynamicCellModelLoader
{
public:
    /**
     * Create a cell model loader by opening a loadable module (.so file) containing
     * a cell model.
     *
     * @note This loader object must remain alive for as long as you want to use cells
     * created with it, or you'll get a segfault.
     *
     * @param rLoadableModulePath  path to .so file
     */
    DynamicCellModelLoader(const std::string& rLoadableModulePath);

    /**
     * Destructor.  Closes the .so file.
     */
    ~DynamicCellModelLoader();

    /**
     * Create a new cardiac cell from this dynamic module.
     *
     * The caller takes responsibility for deleting the cell when it's finished with.
     *
     * @note This loader object must remain alive for as long as you want to use the cell,
     * or you'll get a segfault.
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
