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
     * @return a vector of strings
     * @param rFileName  the name of the file to read from, relative to the output directory
     */
    std::vector<std::string> GetRawDataFromFile(const std::string& rFileName);

public:

    AbstractCachedMeshReader(); /**< Constructor */

    virtual ~AbstractCachedMeshReader()
    {} /**< Destructor. */

    unsigned GetNumElements() const; /**< @return the number of elements in the mesh. */
    unsigned GetNumNodes() const;    /**< @return the number of nodes in the mesh. */
    unsigned GetNumFaces() const;    /**< @return the number of faces in the mesh (synonym of GetNumEdges()) */

    /**
     *  @return the maximum node index. Used in testing to check that output nodes
     *  are always indexed from zero even if they are input indexed from one.
     */
    unsigned GetMaxNodeIndex();

    /**
     *  @return the minimum node index. Used in testing to check that output nodes
     *  are always indexed from zero even if they are input indexed from one.
     */
    unsigned GetMinNodeIndex();

    /**
     *  @return a vector of the coordinates of each node in turn, starting with
     *  node 0 the first time it is called followed by nodes 1, 2, ... , mNumNodes-1.
     */
    std::vector<double> GetNextNode();

    void Reset(); /**< Resets pointers to beginning*/

    /**
     *  @return a vector of the nodes of each element in turn, starting with
     *  element 0 the first time it is called followed by elements 1, 2, ... ,
     *  mNumElements-1.
     */
    ElementData GetNextElementData();

    /**
     *  @return a vector of the nodes of each face in turn, starting with face 0 the
     *  first time it is called followed by faces 1, 2, ... , mNumFaces-1.
     *
     *  Is a synonum of GetNextEdge(). The two functions can be used interchangeably,
     *  i.e. they use the same iterator.
     */
    ElementData GetNextFaceData();
};

#endif /*ABSTRACTCACHEDMESHREADER_HPP_*/
