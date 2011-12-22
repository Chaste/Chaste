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
     * Return value is a vector where each item is a vector of doubles which represents
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
     * Return value is a vector where each item is a vector of ints which represents
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
