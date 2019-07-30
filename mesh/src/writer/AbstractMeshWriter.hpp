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
     * @return the full path to the directory where meshes will be written.
     */
    std::string GetOutputDirectory();

    /**
     * @return the number of nodes in the mesh.
     */
    virtual unsigned GetNumNodes();

    /**
     * @return the number of elements in the mesh.
     */
    unsigned GetNumElements();

    /**
     * @return the number of boundary elements in the mesh.
     */
    unsigned GetNumBoundaryFaces();

    /**
     * @return the number of cable elements in the mesh.
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
