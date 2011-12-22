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

#ifndef CELLMLTOSHAREDLIBRARYCONVERTER_HPP_
#define CELLMLTOSHAREDLIBRARYCONVERTER_HPP_

#include <string>
#include <vector>

#include "DynamicCellModelLoader.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"

/**
 * This class encapsulates all the complexity needed to generate a loadable module from
 * a CellML file.
 */
class CellMLToSharedLibraryConverter
{
public:
    /**
     * Create a converter.
     *
     * @param preserveGeneratedSources  whether to save copies of generated C++
     *    source files in the directory containing the .cellml file.
     * @param component  the name of the Chaste component (or project) in which
     *    to build the loadable module (if required).  Allows projects to have
     *    specialised base classes for dynamically loaded cell models.
     */
    CellMLToSharedLibraryConverter(bool preserveGeneratedSources=false,
                                   std::string component="heart");

    /**
     * Get a loadable module from the given file, and return a loader for it.
     * The file can be a .so, in which case there isn't much to do, just create
     * the loader.  The interesting case comes when it is a .cellml file.  If
     * the file has any other extension, an exception is thrown.
     *
     * @param rFilePath  the model to load
     * @param isCollective  whether this method is being called collectively.
     *   If it is not, then we require the .so to already exist, rather than
     *   trying to avoid race conditions.
     *
     * @note If you do not pass isCollective=false, must be called collectively.
     */
    DynamicCellModelLoader* Convert(const FileFinder& rFilePath,
                                    bool isCollective=true);

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
                           const std::string& rExtraXml="");

private:
    /**
     * Helper method performing the actual conversion of a .cellml file to a .so.
     *
     * @note Must be called collectively.
     *
     * @param rCellmlFullPath  full path to the .cellml file
     * @param rCellmlFolder  folder containing the CellML file, with trailing slash
     * @param rModelLeafName  leaf name of the CellML file, minus extension (but including the .)
     */
    void ConvertCellmlToSo(const std::string& rCellmlFullPath,
                           const std::string& rCellmlFolder,
                           const std::string& rModelLeafName);

    /** Whether to save copies of generated C++ source files. */
    bool mPreserveGeneratedSources;

    /** Which component to build the loadable module in. */
    std::string mComponentName;
};

#endif /*CELLMLTOSHAREDLIBRARYCONVERTER_HPP_*/
