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

#ifndef ABSTRACTMESHWRITER_HPP_
#define ABSTRACTMESHWRITER_HPP_

#include <string>
#include <vector>
#include <iomanip>

#include "OutputFileHandler.hpp"
#include "AbstractMeshReader.hpp"

/**
 * An abstract mesh writer class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractMeshWriter
{
protected: // Give access of these members to subclasses

    OutputFileHandler* mpOutputFileHandler; /**< Output file handler */
    std::string mBaseName; /**< Base name for the input files */

    AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>* mpMeshReader; /**< Writer by default writes from a reader (for conversion).  If this pointer is non-null, data can be copied straight across*/

    unsigned mNumNodes; /**< Total number of nodes in mesh/mesh-reader*/
    unsigned mNumElements; /**< Total number of elements in mesh/mesh-reader*/
    unsigned mNumBoundaryElements; /**< Total number of boundary elements in mesh/mesh-reader*/
    unsigned mNumCableElements; /**< Total number of cable elements in mesh/mesh-reader*/

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    AbstractMeshWriter(const std::string& rDirectory,
                       const std::string& rBaseName,
                       const bool clearOutputDir=true);

    /**
     * Destructor.
     */
    virtual ~AbstractMeshWriter();

    /**
     * Return the full path to the directory where meshes will be written.
     */
    std::string GetOutputDirectory();

    /**
     * Get the number of nodes in the mesh.
     */
    virtual unsigned GetNumNodes();

    /**
     * Get the number of elements in the mesh.
     */
    unsigned GetNumElements();

    /**
     * Get the number of boundary elements in the mesh.
     */
    unsigned GetNumBoundaryFaces();

    /**
     * Get the number of cable elements in the mesh.
     */
    unsigned GetNumCableElements();

    /**
     * @return the coordinates of the next node to be written to file
     */
    virtual std::vector<double> GetNextNode();

    /**
     * @return the data (indices/attributes) of the next element to be written to file
     */
    virtual ElementData GetNextElement();

    /**
     * @return the data (indices/attributes) of the next face to be written to file
     */
    virtual ElementData GetNextBoundaryElement();

    /**
     * @return the data (indices/attributes) of the next cable element to be written to file
     */
    virtual ElementData GetNextCableElement();

    /**
     * Write mesh data to files.
     * This method must be overridden in concrete classes.
     */
    virtual void WriteFiles()=0;

    /**
     * Read in a mesh and write it to file.
     *
     * @param rMeshReader the mesh reader
     */
    void WriteFilesUsingMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader);
};

#endif /*ABSTRACTMESHWRITER_HPP_*/
