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

#ifndef ABSTRACTCACHEDMESHREADER_HPP_
#define ABSTRACTCACHEDMESHREADER_HPP_

#include <vector>
#include <string>

#include "AbstractMeshReader.hpp"

/**
 * Abstract mesh reader class, for readers which read and cache the entire
 * mesh in internal storage, for the mesh to use for constructing itself.
 * Concrete readers which will read large, memory-intensive, meshes should
 * inherit from AbstractMeshReader, not this class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractCachedMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
protected:

    unsigned mNumNodeAttributes; /**< Is the number of attributes stored at each node */
    unsigned mMaxNodeBdyMarker; /**< Is the maximum node boundary marker */
    unsigned mNumElementNodes; /**< Is the number of nodes per element*/
    unsigned mNumElementAttributes; /**< Is the number of attributes stored for each element */
    unsigned mMaxFaceBdyMarker; /**< Is the maximum face (or edge) boundary marker */

    std::vector<std::string> mNodeRawData;  /**< Contents of node input file with comments removed */
    std::vector<std::string> mElementRawData;  /**< Contents of element input file with comments removed */
    std::vector<std::string> mFaceRawData;  /**< Contents of face (or edge) input file with comments removed */

    std::vector< std::vector<double> > mNodeData; /**< Is an array of node coordinates ((i,j)th entry is the jth coordinate of node i)*/
    std::vector< std::vector<unsigned> > mElementData; /**< Is an array of the nodes in each element ((i,j)th entry is the jth node of element i) */
    std::vector< std::vector<unsigned> > mFaceData; /**< Is an array of the nodes in each face ((i,j)th entry is the jth node of face i) */

    std::vector< std::vector<double> >::iterator mpNodeIterator; /**< Is an iterator for the node data */
    std::vector< std::vector<unsigned> >::iterator mpElementIterator; /**< Is an iterator for the element data */
    std::vector< std::vector<unsigned> >::iterator mpFaceIterator; /**< Is an iterator for the face data */

    bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */

    /**
     * Reads an input file rFileName, removes comments (indicated by a #) and blank
     * lines and returns a vector of strings. Each string corresponds to one line
     * of the input file.
     *
     * @param rFileName  the name of the file to read from, relative to the output directory
     */
    std::vector<std::string> GetRawDataFromFile(const std::string& rFileName);

public:

    AbstractCachedMeshReader(); /**< Constructor */

    virtual ~AbstractCachedMeshReader()
    {} /**< Destructor. */

    unsigned GetNumElements() const; /**< Returns the number of elements in the mesh. */
    unsigned GetNumNodes() const;    /**< Returns the number of nodes in the mesh. */
    unsigned GetNumFaces() const;    /**< Returns the number of faces in the mesh (synonym of GetNumEdges()) */

    /**
     *  Returns the maximum node index. Used in testing to check that output nodes
     *  are always indexed from zero even if they are input indexed from one.
     */
    unsigned GetMaxNodeIndex();

    /**
     *  Returns the minimum node index. Used in testing to check that output nodes
     *  are always indexed from zero even if they are input indexed from one.
     */
    unsigned GetMinNodeIndex();

    /**
     *  Returns a vector of the coordinates of each node in turn, starting with
     *  node 0 the first time it is called followed by nodes 1, 2, ... , mNumNodes-1.
     */
    std::vector<double> GetNextNode();

    void Reset(); /**< Resets pointers to beginning*/

    /**
     *  Returns a vector of the nodes of each element in turn, starting with
     *  element 0 the first time it is called followed by elements 1, 2, ... ,
     *  mNumElements-1.
     */
    ElementData GetNextElementData();

    /**
     *  Returns a vector of the nodes of each face in turn, starting with face 0 the
     *  first time it is called followed by faces 1, 2, ... , mNumFaces-1.
     *
     *  Is a synonum of GetNextEdge(). The two functions can be used interchangeably,
     *  i.e. they use the same iterator.
     */
    ElementData GetNextFaceData();

};

#endif /*ABSTRACTCACHEDMESHREADER_HPP_*/
