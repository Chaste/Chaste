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

#include "EdgeHelper.hpp"

template<unsigned SPACE_DIM>
EdgeHelper<SPACE_DIM>::EdgeHelper()
{
}

template<unsigned SPACE_DIM>
Edge<SPACE_DIM>* EdgeHelper<SPACE_DIM>::GetEdgeFromNodes(Node<SPACE_DIM>* pNodeA,
                                                         Node<SPACE_DIM>* pNodeB)
{
    auto edge_map_indices = Edge<SPACE_DIM>::GenerateMapIndex(pNodeA->GetIndex(), pNodeB->GetIndex());

    // Check that an edge hasn't been created already
    Edge<SPACE_DIM>* p_edge = nullptr;

    auto edge_iter = mEdgesMap.find(edge_map_indices);
    if (edge_iter == mEdgesMap.end() || edge_iter->second->IsDeleted())
    {
        mEdges.push_back(std::make_unique<Edge<SPACE_DIM>> (mEdges.size(), pNodeA, pNodeB));
        p_edge = mEdges.back().get();
        mEdgesMap[edge_map_indices] = p_edge;
    }
    else
    {
        p_edge = edge_iter->second;
        p_edge->SetNodes(pNodeA, pNodeB);
    }

    return p_edge;
}

template<unsigned SPACE_DIM>
Edge<SPACE_DIM> *
EdgeHelper<SPACE_DIM>::GetEdgeFromNodes(unsigned elementIndex,
                                        Node<SPACE_DIM>* pNodeA,
                                        Node<SPACE_DIM>* pNodeB)
{
    auto p_edge = GetEdgeFromNodes(pNodeA, pNodeB);
    p_edge->AddElement(elementIndex);
    return p_edge;
}

template<unsigned SPACE_DIM>
Edge<SPACE_DIM>* EdgeHelper<SPACE_DIM>::GetEdge(unsigned index) const
{
    assert(index < mEdges.size());
    return mEdges[index].get();
}

template<unsigned SPACE_DIM>
Edge<SPACE_DIM>* EdgeHelper<SPACE_DIM>::operator[](unsigned index) const
{
    assert(index < mEdges.size());
    return mEdges[index].get();
}

template<unsigned SPACE_DIM>
void EdgeHelper<SPACE_DIM>::RemoveDeletedEdges()
{
    /*
     * Remove any edges that have been marked for deletion and store all other 
     * nodes in a temporary structure.
     */
    std::vector<std::unique_ptr<Edge<SPACE_DIM> > > live_edges;

    for (auto& p_edge : mEdges)
    {
      // An edge can be deleted if it is not contained in any elements
      if (p_edge->GetNumElements() != 0)
      {
          live_edges.push_back(std::move(p_edge));
      }
      /*
       * Note: there's no need to manually delete the edge, because the 
       * unique_ptr automatically takes care of that when it goes out of scope.
       */
    }

    /*
     * Repopulate the edges vector, using std::move to efficiently transfer 
     * ownership of the pointers.
     */
    mEdges = std::move(live_edges);

    // Reset the edge indices
    for (unsigned i = 0; i < this->mEdges.size(); ++i)
    {
        mEdges[i]->SetIndex(i);
    }

    UpdateEdgesMapKey();
}

template<unsigned SPACE_DIM>
void EdgeHelper<SPACE_DIM>::UpdateEdgesMapKey()
{
    mEdgesMap.clear();

    for (auto& p_edge : mEdges) 
    {
        mEdgesMap[p_edge->GetMapIndex()] = p_edge.get();
    }
}

template<unsigned SPACE_DIM>
unsigned EdgeHelper<SPACE_DIM>::GetNumEdges() const
{
    return mEdges.size();
}

// Explicit instantiation
template class EdgeHelper<1>;
template class EdgeHelper<2>;
template class EdgeHelper<3>;