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


#ifndef VERTEXELEMENTMAP_HPP_
#define VERTEXELEMENTMAP_HPP_

#include <vector>
#include "Exception.hpp"

/**
 * VertexElementMap class used when remeshing. The map associates the indices
 * of VertexElements in the old mesh with indices of VertexElements in the new mesh.
 */
class VertexElementMap
{
private:

    /** The map is stored as an ordered vector of VertexElement indices. */
    std::vector<unsigned> mMap;

public:

    /**
     * Constructor.
     *
     * @param size  the size of the VertexElementMap
     */
    VertexElementMap(unsigned size);

    /**
     * Resize the VertexElementMap.
     *
     * @param size  the new size of the VertexElementMap
     */
    void Resize(unsigned size);

    /**
     * Reset the VertexElementMap to the identity map.
     */
    void ResetToIdentity();

    /**
     * Associate a given old index with a new index.
     *
     * @param oldIndex  the old index of a VertexElement
     * @param newIndex  the new index of a VertexElement
     */
    void SetNewIndex(unsigned oldIndex, unsigned newIndex);

    /**
     * Mark a given old index as 'deleted' by associating it
     * with the new index UINT_MAX.
     *
     * @param index  the old index of a VertexElement
     */
    void SetDeleted(unsigned index);

    /**
     * @return whether a given old index is marked as 'deleted'.
     *
     * @param index  the old index of a VertexElement
     */
    bool IsDeleted(unsigned index);

    /**
     * @return the new index associated with a given old index.
     *
     * @param oldIndex  the old index of a VertexElement
     */
    unsigned GetNewIndex(unsigned oldIndex) const;

    /**
     * @return whether the VertexElementMap is the identity map.
     */
    bool IsIdentityMap();

    /**
     * @return the size of the VertexElementMap.
     */
    unsigned Size();
};

#endif /*VERTEXELEMENTMAP_HPP_*/
