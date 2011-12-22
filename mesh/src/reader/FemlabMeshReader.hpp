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
     * Return value is a vector where each item is a vector of double which represents
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
     * Return value is a vector where each item is a vector of ints which represents
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
