/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef EDGEHELPER_HPP_
#define EDGEHELPER_HPP_

#include <vector>
#include <map>
#include "Node.hpp"
#include "Edge.hpp"
#include "EdgeOperation.hpp"
#include "EdgeRemapInfo.hpp"

/**
 * Class for facilitating the creation and management of unique edges in a 
 * vertex mesh.
 */
template <unsigned SPACE_DIM>
class EdgeHelper
{
private:

    /**
     * Vector owning the individual edge objects.
     */
    std::vector<std::unique_ptr<Edge<SPACE_DIM>>> mEdges;

    /**
     * Explicit map between the two node global indices and the edge they represent.
     */
    std::map<std::pair<unsigned, unsigned>, Edge<SPACE_DIM>*> mEdgesMap;

    /**
     * Rebuilds node-node to edge map, which is required after removing deleted edges.
     */
    void UpdateEdgesMapKey();

public:

    /**
     * Default constructor.
     */
    EdgeHelper();

    /**
     * Get edge from the node pairs. Construct the edge if it has not been created
     * @param pNodeA pointer to first Node
     * @param pNodeB pointer to second Node
     * 
     * @return the (new) Edge
     */
    Edge<SPACE_DIM>* GetEdgeFromNodes(Node<SPACE_DIM>* pNodeA,
                                      Node<SPACE_DIM>* pNodeB);

    /**
     * Get the edge from the node pairs and add it to element with index elementIndex
     * @param elementIndex the index of an element to which the edge belongs
     * @param pNodeA pointer to first Node
     * @param pNodeB pointer to second Node
     * 
     * @return the Edge
     */
    Edge<SPACE_DIM>* GetEdgeFromNodes(unsigned elementIndex,
                                      Node<SPACE_DIM>* pNodeA,
                                      Node<SPACE_DIM>* pNodeB);

    /**
     * @param index a global Edge index
     * 
     * @return pointer the Edge with this global index
     */
    Edge<SPACE_DIM>* GetEdge(unsigned index) const;

    /**
     * Access operator.
     * 
     * @param index Index of mEdges
     * 
     * @return mEdges[index]
     */
    Edge<SPACE_DIM>* operator[](unsigned index) const;

    /**
     * Remove deleted edges.
     */
    void RemoveDeletedEdges();

    /**
     * @return total number of edges
     */
    unsigned GetNumEdges() const;
};

#endif /* CHASTE_EDGEHELPER_HPP_ */