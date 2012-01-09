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

#ifndef DYNAMICLOADINGHELPERFUNCTIONS_HPP_
#define DYNAMICLOADINGHELPERFUNCTIONS_HPP_

#include <string>
#include <vector>

#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "Exception.hpp"
#include "DynamicCellModelLoader.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "SimpleStimulus.hpp"
#include "AbstractDynamicallyLoadableEntity.hpp"

/**
 * Create a PyCml options file for the given model.
 *
 * @param rHandler  where to create the file
 * @param rModelName  base name of the model file (which will be "rModelName.cellml")
 * @param rArgs  extra command-line arguments for the model conversion
 * @param rExtraXml  any extra XML to go in the config file (e.g. LT settings)
 */
void CreateOptionsFile(const OutputFileHandler& rHandler,
                       const std::string& rModelName,
                       const std::vector<std::string>& rArgs,
                       const std::string& rExtraXml="")
{
    if (PetscTools::AmMaster())
    {
        out_stream p_optfile = rHandler.OpenOutputFile(rModelName + "-conf.xml");
        (*p_optfile) << "<?xml version='1.0'?>" << std::endl
                     << "<pycml_config>" << std::endl
                     << "<command_line_args>" << std::endl;
        for (unsigned i=0; i<rArgs.size(); i++)
        {
            (*p_optfile) << "<arg>" << rArgs[i] << "</arg>" << std::endl;
        }
        (*p_optfile) << "</command_line_args>" << std::endl
                     << rExtraXml
                     << "</pycml_config>" << std::endl;
        p_optfile->close();
    }
    PetscTools::Barrier("CreateOptionsFile");
}

/**
 * Copy a file.
 *
 * @param rDestDir  the folder to copy to
 * @param rSourceFile  the file to copy
 */
void CopyFile(const OutputFileHandler& rDestDir,
              const FileFinder& rSourceFile)
{
    rDestDir.CopyFileTo(rSourceFile);
}

/**
 * Create a new cell object from a cell model loader.
 *
 * @param rLoader  the loader
 * @param magnitude  the magnitude of the stimulus to apply, in uA/cm^2 (duration of 2ms and when of 50ms are fixed)
 */
AbstractCardiacCellInterface* CreateCellWithStandardStimulus(DynamicCellModelLoader& rLoader,
                                                             double magnitude=-25.5) // uA/cm^2
{
    // Set stimulus
    double duration  = 2.0; // ms
    double when = 50.0; // ms
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new SimpleStimulus(magnitude, duration, when));
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

    // Load the cell model dynamically
    AbstractCardiacCellInterface* p_cell = rLoader.CreateCell(p_solver, p_stimulus);

    // Simple sanity checks
    AbstractDynamicallyLoadableEntity* p_entity = dynamic_cast<AbstractDynamicallyLoadableEntity*>(p_cell);
    TS_ASSERT(p_entity != NULL);
    if (p_entity != NULL)
    {
        TS_ASSERT_EQUALS(&rLoader, p_entity->GetLoader());
    }

    return p_cell;
}

#endif // DYNAMICLOADINGHELPERFUNCTIONS_HPP_
