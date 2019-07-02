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

#ifndef EDGE_HPP_
#define EDGE_HPP_

#include <set>
#include <vector>

#include "Node.hpp"

typedef std::pair<unsigned ,unsigned> UIndexPair;

/**
 * An edge in a finite element mesh
 */
template<unsigned SPACE_DIM>
class Edge
{

private:

    /** Index of this edge within the mesh **/
    unsigned mIndex;

    bool mIsDeleted;

    /** Nodes that form this edge **/
    std::vector<Node<SPACE_DIM>*> mNodes;

    /** Elements that this edge belongs to **/
    std::set<unsigned> mElementIndices;

public:

    ///\todo #2987 Document method
    Edge(unsigned index);

    ///\todo #2987 Document method
    Edge(unsigned index, Node<SPACE_DIM>* pNode0, Node<SPACE_DIM>* pNode1);

    /**
     * Destructor.
     */
    ~Edge();

    ///\todo #2987 Document method
    void MarkDeleted();

    ///\todo #2987 Document method
    bool IsDeleted();

    ///\todo #2987 Document method
    void SetIndex(unsigned index);

    ///\todo #2987 Document method
    unsigned GetIndex();

    ///\todo #2987 Document method
    UIndexPair GetMapIndex();

    ///\todo #2987 Document method
    void RemoveNodes();

    ///\todo #2987 Document method
    void SetNodes(Node<SPACE_DIM>* pNode0, Node<SPACE_DIM>* pNode1);

    ///\todo #2987 Document method
    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode);

    ///\todo #2987 Document method
    Node<SPACE_DIM>* GetNode(unsigned index);

    ///\todo #2987 Document method
    unsigned GetNumNodes();

    ///\todo #2987 Document method
    bool ContainsNode(Node<SPACE_DIM>* pNode);

    ///\todo #2987 Document method
    c_vector<double, SPACE_DIM> rGetCentreLocation();

    ///\todo #2987 Document method
    double rGetLength();

    ///\todo #2987 Document method
    std::set<unsigned> GetOtherElements(unsigned elementIndex);

    ///\todo #2987 Document method
    void AddElement(unsigned elementIndex);

    ///\todo #2987 Document method
    void RemoveElement(unsigned elementIndex);

    ///\todo #2987 Document method
    std::set<unsigned> GetNeighbouringElementIndices();

    ///\todo #2987 Document method
    unsigned GetNumElements();

    ///\todo #2987 Document method
    bool IsEdgeValid();
};

#endif //EDGE_HPP_
