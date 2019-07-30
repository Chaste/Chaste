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

#ifndef _FEMLABMESHREADER_H_
#define _FEMLABMESHREADER_H_

#include "AbstractCachedMeshReader.hpp"

/**
 * Concrete version of the AbstractCachedMeshReader class.
 * A FemlabMeshReader takes the file names of a set of Femlab mesh files.
 * Once constructed the public methods of the AbstractCachedMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class FemlabMeshReader : public AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:

    /**
     * TokenizeStringsToDoubles is specific to reading node data which came from
     * a Femlab or Matlab PDE toolbox file.
     *
     * Each string is expected to be a series of doubles.
     * @return a vector where each item is a vector of double which represents
     * position.  Indices are implicit in the vector.
     *
     * @param rRawData the node data to be read
     */
    std::vector<std::vector<double> > TokenizeStringsToDoubles(const std::vector<std::string>& rRawData);

    /**
     * TokenizeStringsToInts is for reading element, face or edge data which came from
     * a Femlab or Matlab PDE toolbox file.
     *  Each string is expected to be a series of unsigned which represent:
     *  The first several lines denote the indices of nodes
     *  The rest contains extra information which are ignored currently.
     *  ( In 2-D: 2 indices for an edge, 3 for a triangle)
     *  ( In 3-D: 3 indices for a face, 4 for a tetrahedron)
     * @return a vector where each item is a vector of ints which represents
     * indices of nodes.
     *
     * @param rRawData  the element, face or edge data to be read
     * @param dimensionOfObject  the number of lines of data to be read
     */
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(const std::vector<std::string>& rRawData,
                                                              unsigned dimensionOfObject);

public:

    /**
     * The constructor takes the path to and names of a set of Femlab mesh files
     * (ie. the nodes, elements and faces files (in that order) and allows the data to
     * be queried.
     * Typical use:
     *    AbstractMeshReader* pMeshReader = new FemlabMeshReader("pdes/tests/meshdata/",
     *                                                           "femlab_lshape_nodes.dat",
     *                                                           "femlab_lshape_elements.dat",
     *                                                           "femlab_lshape_edges.dat",);
     *
     * @param rPathBaseName  the base name of the files from which to read the mesh data
     * @param rNodeFileName  the name of the nodes file
     * @param rElementFileName  the name of the elements file
     * @param rEdgeFileName  the name of the edges file
     */
    FemlabMeshReader(const std::string& rPathBaseName,
                     const std::string& rNodeFileName,
                     const std::string& rElementFileName,
                     const std::string& rEdgeFileName);

    /**
     * Destructor
     */
    virtual ~FemlabMeshReader();
};

#endif //_FEMLABMESHREADER_H_
