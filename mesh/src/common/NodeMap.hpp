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


#ifndef NODEMAP_HPP_
#define NODEMAP_HPP_

#include <vector>

/**
 * Nodemap class used when remeshing. The map associates the indices of nodes
 * in the old mesh with indices of nodes in the new mesh.
 */
class NodeMap
{
private:

    /** The map is stored as an ordered vector of node indices. */
    std::vector<unsigned> mMap;

public:

    /**
     * Constructor.
     *
     * @param size  the size of the NodeMap
     */
    NodeMap(unsigned size);

    /**
     * Resize the NodeMap.
     *
     * @param size  the new size of the NodeMap
     */
    void Resize(unsigned size);

    /**
     * Reset the NodeMap to the identity map.
     */
    void ResetToIdentity();

    /**
     * Associate a given old index with a new index.
     *
     * @param oldIndex  the old index of a node
     * @param newIndex  the new index of a node
     */
    void SetNewIndex(unsigned oldIndex, unsigned newIndex);

    /**
     * Mark a given old index as 'deleted' by associating it
     * with the new index UINT_MAX.
     *
     * @param index  the old index of a node
     */
    void SetDeleted(unsigned index);

    /**
     * Get whether a given old index is marked as 'deleted'.
     *
     * @param index  the old index of a node
     */
    bool IsDeleted(unsigned index);

    /**
     * Get the new index associated with a given old index.
     *
     * @param oldIndex  the old index of a node
     */
    unsigned GetNewIndex(unsigned oldIndex) const;

    /**
     * Get whether the NodeMap is the identity map.
     */
    bool IsIdentityMap();

    /**
     * Get the size of the NodeMap.
     */
    unsigned Size();

};

#endif /*NODEMAP_HPP_*/
