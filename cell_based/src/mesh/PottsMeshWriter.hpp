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

#ifndef POTTSMESHWRITER_HPP_
#define POTTSMESHWRITER_HPP_

// Forward declaration prevents circular include chain
template<unsigned SPACE_DIM>
class PottsMesh;

#include "PottsMesh.hpp"
#include "AbstractMeshWriter.hpp"
#include "NodeMap.hpp"

template<unsigned SPACE_DIM>
struct MeshPottsWriterIterators;

/**
 * A mesh writer class for Potts meshes.
 */
template<unsigned SPACE_DIM>
class PottsMeshWriter : public AbstractMeshWriter<SPACE_DIM, SPACE_DIM>
{
private:

    /**
     * If writing from a mesh object, the mesh to write to disk.
     * Otherwise NULL.
     */
    PottsMesh<SPACE_DIM>* mpMesh;

    /** Iterators over the mesh */
    MeshPottsWriterIterators<SPACE_DIM>* mpIters;

    /** Track deleted nodes so they don't get written */
    NodeMap* mpNodeMap;

    /** What was the last index written to #mpNodeMap ? */
    unsigned mNodeMapCurrentIndex;

public:

    /**
     * Constructor.
     *
     * @param rDirectory reference to the output directory, relative to where Chaste output is stored
     * @param rBaseName reference to the base name for results files
     * @param clearOutputDir whether to clear the output directory prior to writing files (defaults to true)
     */
    PottsMeshWriter(const std::string& rDirectory,
                     const std::string& rBaseName,
                     const bool clearOutputDir=true);

    /**
     * Destructor.
     */
    ~PottsMeshWriter();

    ///\todo Mesh should be const (#1076)
    /**
     * Write files using a mesh.
     *
     * @param rMesh reference to the Potts mesh
     */
    void WriteFilesUsingMesh(PottsMesh<SPACE_DIM>& rMesh);

    /**
     * @return the coordinates of the next node to be written to file
     */
    std::vector<double> GetNextNode();

    /**
     * @return the data (indices/attributes) of the next element to be written to file
     */
    ElementData GetNextElement();

    /**
     * Write mesh data to files.
     * This method must be overridden in concrete classes.
     */
    void WriteFiles();
};

#endif /*POTTSMESHWRITER_HPP_*/
