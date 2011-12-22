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

#ifndef ABSTRACTHDF5CONVERTER_HPP_
#define ABSTRACTHDF5CONVERTER_HPP_

#include <string>
#include "AbstractTetrahedralMesh.hpp"
#include "OutputFileHandler.hpp"
#include "Hdf5DataReader.hpp"

/**
 * The derived children of this class convert output from Hdf5 format to
 * a range of other formats for postprocessing.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractHdf5Converter
{
protected:

    /** Pointer to reader of the file to be converted. */
    Hdf5DataReader* mpReader;

    /** Number of variables to output. Read from the reader. */
    unsigned mNumVariables;

    /** Base name for the files: [basename].vtu, [basename].dat etc.*/
    std::string mFileBaseName;

    /** Pointer to a mesh. */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /** Initialised as the directory in which to store the results. */
    OutputFileHandler* mpOutputFileHandler;

    /**
     * Get the subdirectory in which the converted output is stored,
     * relative to the input directory.
     */
    std::string mRelativeSubdirectory;

public:

    /**
     * Constructor, which does the conversion and writes the .info file.
     *
     * @note This method is collective, and must be called by al processes.
     *
     * @param inputDirectory The input directory, relative to CHASTE_TEST_OUTPUT, where the .h5 file has been written
     * @param fileBaseName The base name of the data file.
     * @param pMesh Pointer to the mesh.
     * @param subdirectoryName name for the output directory to be created (relative to inputDirectory)
     */
    AbstractHdf5Converter(std::string inputDirectory,
                          std::string fileBaseName,
                          AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                          std::string subdirectoryName);

    /**
     * Destructor.
     */
    ~AbstractHdf5Converter();

    /**
     * Get the relative path of the subdirectory in which the converted output is stored.
     */
    std::string GetSubdirectory();
};

#endif /*ABSTRACTHDF5CONVERTER_HPP_*/
