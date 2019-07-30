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

#ifndef _MEMFEMMESHREADER_HPP_
#define _MEMFEMMESHREADER_HPP_

#include "AbstractCachedMeshReader.hpp"

/**
 * Concrete version of the AbstractCachedMeshReader class.
 * A MemfemMeshReader takes the base name of a set of Memfem
 * mesh files (ie. the path and name of the files without the suffices).
 * Once constructed the public methods of the AbstractCachedMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MemfemMeshReader : public AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:

    /**
     * TokenizeStringsToDoubles is specific to reading node data which came from
     * a Memfem file.
     * Each string is expected to be 3 doubles (representing x,y,z)
     * @return a vector where each item is a vector of doubles which represents
     * position.  Indices are implicit in the vector.
     *
     * @param rRawData  the node data to be read
     */
    std::vector<std::vector<double> > TokenizeStringsToDoubles(const std::vector<std::string>& rRawData);

    /**
     * TokenizeStringsToInts is for reading element or boundary face data which came from
     * a Memfem file.
     *  Each string is expected to be:
     *  3 or 4 node indices
     *  ( 3 indices for a face, 4 for a tetrahedron)
     *  a region marker? (if it's an element)
     *  NB: Region markers are currently ignored.
     * @return a vector where each item is a vector of ints which represents
     * indices of nodes.
     *
     * @param rRawData  the element or boundary face data to be read
     * @param dimensionOfObject  the number of lines of data to be read
     * @param readHeader  whether to read the header
     */
    std::vector<std::vector<unsigned> > TokenizeStringsToInts(const std::vector<std::string>& rRawData,
                                                              unsigned dimensionOfObject,
                                                              bool readHeader);

public:

    /**
     * The constructor takes the base name of a set of Memfem
     * mesh files (ie. the path and name of the files without the suffices)
     * and allows the data to be queried.
     *
     * Typical use:
     *    AbstractMeshReader* pMeshReader = new MemfemMeshReader("pdes/tests/meshdata/Memfem_slab");
     *
     * @param rPathBaseName  the base name of the files from which to read the mesh data
     */
    MemfemMeshReader(const std::string& rPathBaseName);

    /**
     * Destructor.
     */
    virtual ~MemfemMeshReader();
};

#endif //_MEMFEMMESHREADER_HPP_
