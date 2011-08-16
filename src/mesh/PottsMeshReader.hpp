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

#ifndef POTTSMESHREADER_HPP_
#define POTTSMESHREADER_HPP_

#include <string>
#include <fstream>
#include <cassert>
#include <vector>

#include "Exception.hpp"
#include "AbstractMeshReader.hpp"

/**
 * Helper structure that stores the nodes and any attribute value
 * associated with a PottsElement.
 */
struct PottsElementData
{
    std::vector<unsigned> NodeIndices; /**< Vector of Node indices owned by the element. */
    std::vector<bool> Orientations; /**< Vector of face orientations (only used in 3D). */
    unsigned AttributeValue; /**< Attribute value associated with the element. */
    unsigned ContainingElement; /**< Only applies to boundary elements: which element contains this boundary element. Only set if reader called with correct params */
};

/**
 * A mesh reader class for potts-based meshes. So far implemented in 2D only.
 */
template<unsigned SPACE_DIM>
class PottsMeshReader : public AbstractMeshReader<SPACE_DIM, SPACE_DIM>
{
private:

    /** The base name for mesh files. */
    std::string mFilesBaseName;

    /** The nodes file for the mesh. */
    std::ifstream mNodesFile;

    /** The elements file for the mesh. */
    std::ifstream mElementsFile;

    /** True if input data is numbered from zero, false otherwise. */
    bool mIndexFromZero;

    /** Number of nodes in the mesh. */
    unsigned mNumNodes;

    /** Number of elements in the mesh. */
    unsigned mNumElements;

    /** Number of nodes read in by the reader. */
    unsigned mNodesRead;

    /** Number of elements read in by the reader. */
    unsigned mElementsRead;

    /** Is the number of attributes stored at each node. */
    unsigned mNumNodeAttributes;

    /** Is the number of attributes stored for each element. */
    unsigned mNumElementAttributes;

    /**
     * Open node and element files.
     */
    void OpenFiles();

    /**
     * Open node file.
     */
    void OpenNodeFile();

    /**
     * Open element file.
     */
    void OpenElementsFile();

    /**
     * Read the file headers to determine node and element numbers and attributes.
     */
    void ReadHeaders();

    /**
     * Close node and element files.
     */
    void CloseFiles();

    /**
     * Get the next line from a given file stream.
     *
     * @param fileStream the file stream
     * @param rawLine the raw line (may contain comments)
     */
    void GetNextLineFromStream(std::ifstream& fileStream, std::string& rawLine);

public:

    /**
     * Constructor.
     *
     * @param pathBaseName the base name for results files
     */
    PottsMeshReader(std::string pathBaseName);

    /**
     * Destructor.
     */
    ~PottsMeshReader()
    {}

    /**
     * @return the number of elements in the mesh.
     */
    unsigned GetNumElements() const;

    /**
     * @return the number of nodes in the mesh.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of attributes in the mesh
     */
    unsigned GetNumElementAttributes() const;

    /**
     * Returns the number of faces in the mesh (synonym of GetNumEdges()).
     */
    unsigned GetNumFaces() const;

    /**
     * Returns a vector of the nodes of each face in turn.
     */
    ElementData GetNextFaceData();

    /**
     * Reset pointers to beginning.
     */
    void Reset();

    /**
     * @return the coordinates of each node in turn.
     */
    std::vector<double> GetNextNode();

    /**
     * @return the nodes of each element (and any attribute infomation, if there is any) in turn.
     */
    ElementData GetNextElementData();

    /**
     * @return the nodes of each element (and any attribute infomation, if there is any) in turn, then its faces.
     *         This method should only be called in 3D.
     */
    PottsElementData GetNextElementDataWithFaces();
};

#endif /*POTTSMESHREADER_HPP_*/
