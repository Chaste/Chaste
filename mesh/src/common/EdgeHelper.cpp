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

#include "EdgeHelper.hpp"

template<unsigned int SPACE_DIM>
EdgeHelper<SPACE_DIM>::EdgeHelper()
{
}

template<unsigned int SPACE_DIM>
EdgeHelper<SPACE_DIM>::~EdgeHelper()
{
    for (unsigned int i=0; i<mEdges.size(); ++i)
    {
        delete mEdges[i];
    }
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::Clear()
{
    // Iterate over edges and free the memory
    for (auto edge: mEdges)
    {
        delete edge;
    }
    mEdges.clear();
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::GetEdgeFromNodes(Node<SPACE_DIM> *node0, Node<SPACE_DIM> *node1)
{
    // Swap node so we always have the lower index as node 0
    if (node0->GetIndex() > node1->GetIndex())
    {
        auto swapNode = node0;
        node0 = node1;
        node1 = swapNode;
    }

    auto edgeMapIndices = UIndexPair(node0->GetIndex(), node1->GetIndex());

    // Check that an edge hasn't been created already
    Edge<SPACE_DIM>* edge = nullptr;

    auto edgeItt = mEdgesMap.find(edgeMapIndices);
    if (edgeItt == mEdgesMap.end() || edgeItt->second->IsDeleted())
    {
        edge = new Edge<SPACE_DIM>(mEdges.size(), node0, node1);
        mEdgesMap[edgeMapIndices] = edge;
        mEdges.push_back(edge);
    }
    else
    {
        edge = edgeItt->second;
        edge->SetNodes(node0, node1);
    }

    return edge;
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *
EdgeHelper<SPACE_DIM>::GetEdgeFromNodes(unsigned elementIndex, Node<SPACE_DIM> *node0, Node<SPACE_DIM> *node1)
{
    auto edge = GetEdgeFromNodes(node0, node1);
    edge->AddElement(elementIndex);
    return edge;
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::GetEdge(unsigned index)
{
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::GetEdge(unsigned index) const
{
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::operator[](unsigned index)
{
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::operator[](unsigned index) const
{
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::RemoveDeletedEdges()
{
    // Remove any edges that have been marked for deletion and store all other nodes in a temporary structure
    std::vector<Edge<SPACE_DIM>*> live_edges;
    for (unsigned i=0; i<this->mEdges.size(); i++)
    {
        // An edge can be deleted if it is not contained in any elements
        if (this->mEdges[i]->GetNumElements() == 0)
        {
            delete this->mEdges[i];
        }
        else
        {
            live_edges.push_back(this->mEdges[i]);
        }
    }

    // Repopulate the edges vector
    this->mEdges = live_edges;

    for (unsigned i=0; i<this->mEdges.size(); i++)
    {
        this->mEdges[i]->SetIndex(i);
    }
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::UpdateEdgesMapKey()
{
    mEdgesMap.clear();

    for (auto edge : mEdges)
    {
        mEdgesMap[edge->GetMapIndex()] = edge;
    }
}

template<unsigned int SPACE_DIM>
unsigned EdgeHelper<SPACE_DIM>::GetNumEdges() const {
    return mEdges.size();
}

template class EdgeHelper<1>;
template class EdgeHelper<2>;
template class EdgeHelper<3>;
