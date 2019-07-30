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

    /** The map is stored as an ordered vector of node indices.
     *  Note that a vector is more efficient than std::map in this case where
     *  we know the exact size of the map and can benefit from O(1) lookup*/
    std::vector<unsigned> mMap;
    /** Redundant data to provide a shortcut when the map is empty or the identity.
     *  This is initialised to true, since an empty map is equivalent to the identity. */
    bool mIsIdentity;

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
     * @return whether a given old index is marked as 'deleted'.
     *
     * @param index  the old index of a node
     */
    bool IsDeleted(unsigned index);

    /**
     * @return the new index associated with a given old index.
     *
     * @param oldIndex  the old index of a node
     */
    unsigned GetNewIndex(unsigned oldIndex) const;

    /**
     * @return whether the NodeMap is the identity map.
     */
    bool IsIdentityMap();

    /**
     * @return the size of the NodeMap.
     */
    unsigned GetSize();
};

#endif /*NODEMAP_HPP_*/
