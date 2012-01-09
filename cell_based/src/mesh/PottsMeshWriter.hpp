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
 * A mesh writer class for potts-based meshes.
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
     * @param rMesh reference to the potts-based mesh
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
