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
    /** Test gets access to the msSoSuffix variable */
    friend class TestDynamicallyLoadedCellModels;
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
     * @return a loadable module from the given file, and return a loader for it.
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
    DynamicCellModelLoaderPtr Convert(const FileFinder& rFilePath,
                                      bool isCollective=true);

    /**
     * Create a PyCml options file for the given model.
     *
     * @param rHandler  where to create the file
     * @param rModelName  base name of the model file (which will be "rModelName.cellml")
     * @param rArgs  extra command-line arguments for the model conversion
     * @param rExtraXml  any extra XML to go in the config file (e.g. LT settings)
     */
    static void CreateOptionsFile(const OutputFileHandler& rHandler,
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
     */
    void ConvertCellmlToSo(const std::string& rCellmlFullPath,
                           const std::string& rCellmlFolder);

    /** Whether to save copies of generated C++ source files. */
    bool mPreserveGeneratedSources;

    /** Which component to build the loadable module in. */
    std::string mComponentName;

    /** The .so suffix is nearly always "so" (as you might expect).  On Mac OSX this is redefined to "dylib" */
    static const std::string msSoSuffix;
};

#endif /*CELLMLTOSHAREDLIBRARYCONVERTER_HPP_*/
