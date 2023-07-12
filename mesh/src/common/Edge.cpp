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

#include "Edge.hpp"

template <unsigned SPACE_DIM>
Edge<SPACE_DIM>::Edge(unsigned index)
    : mIndex(index),
      mIsDeleted(false)
{
}

template <unsigned SPACE_DIM>
Edge<SPACE_DIM>::Edge(unsigned index,
                      Node<SPACE_DIM>* pNodeA,
                      Node<SPACE_DIM>* pNodeB)
    : mIndex(index),
      mIsDeleted(false)
{
    this->SetNodes(pNodeA, pNodeB);
}

template <unsigned SPACE_DIM>
std::pair<unsigned, unsigned> Edge<SPACE_DIM>::GenerateMapIndex(unsigned indexA,
                                                                unsigned indexB)
{
    return std::make_pair(std::min(indexA, indexB), std::max(indexA, indexB));
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::MarkAsDeleted()
{
    mIsDeleted = true;
}

template<unsigned SPACE_DIM>
bool Edge<SPACE_DIM>::IsDeleted()
{
    return mIsDeleted;
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::SetIndex(unsigned index)
{
    mIndex = index;
}

template<unsigned SPACE_DIM>
unsigned Edge<SPACE_DIM>::GetIndex() const
{
    return mIndex;
}

template<unsigned SPACE_DIM>
std::pair<unsigned ,unsigned> Edge<SPACE_DIM>::GetMapIndex()
{
    assert(mNodes.size() == 2);
    const unsigned index0 = mNodes[0]->GetIndex();
    const unsigned index1 = mNodes[1]->GetIndex();
    return Edge<SPACE_DIM>::GenerateMapIndex(index0, index1);
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::RemoveNodes()
{
    mNodes.clear();
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::SetNodes(Node<SPACE_DIM>* pNodeA,
                               Node<SPACE_DIM>* pNodeB)
{
    // Clear the nodes first
    this->RemoveNodes();

    // Add nodes
    this->mNodes.push_back(pNodeA);
    this->mNodes.push_back(pNodeB);
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::ReplaceNode(Node<SPACE_DIM>* pOldNode,
                                  Node<SPACE_DIM>* pNewNode)
{
    for (unsigned i = 0; i < 2; i++)
    {
        if (this->mNodes[i] == pOldNode)
        {
            this->mNodes[i] = pNewNode;
            return;
        }
    }
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>* Edge<SPACE_DIM>::GetNode(unsigned index) const
{
    return mNodes[index];
}

template<unsigned SPACE_DIM>
unsigned Edge<SPACE_DIM>::GetNumNodes()
{
    return mNodes.size();
}

template<unsigned SPACE_DIM>
bool Edge<SPACE_DIM>::ContainsNode(Node<SPACE_DIM>* pNode) const
{
    for (auto p_node : mNodes)
    {
        if (p_node->GetIndex() == pNode->GetIndex())
        {
            return true;
        }
    }
    return false;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> Edge<SPACE_DIM>::rGetCentreLocation()
{
    assert(mNodes.size() == 2);
    return 0.5*(mNodes[0]->rGetLocation() + mNodes[1]->rGetLocation());
}

template<unsigned SPACE_DIM>
double Edge<SPACE_DIM>::rGetLength()
{
    assert(mNodes.size() == 2);
    return norm_2(mNodes[0]->rGetLocation() - mNodes[1]->rGetLocation());
}

template<unsigned SPACE_DIM>
std::set<unsigned> Edge<SPACE_DIM>::GetOtherElements(unsigned elementIndex)
{
    std::set<unsigned> otherElements;
    for (unsigned elem: mElementIndices)
    {
        if (elem != elementIndex)
        {
            otherElements.insert(elem);
        }
    }
    return otherElements;
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::AddElement(unsigned elementIndex)
{
    mElementIndices.insert(elementIndex);
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::RemoveElement(unsigned elementIndex)
{
    mElementIndices.erase(elementIndex);
}

template<unsigned SPACE_DIM>
std::set<unsigned> Edge<SPACE_DIM>::GetNeighbouringElementIndices()
{
    std::set<unsigned> neighbouring_element_indices;

    auto elem_indices0 = mNodes[0]->rGetContainingElementIndices();
    auto elem_indices1 = mNodes[1]->rGetContainingElementIndices();

    std::set_intersection(elem_indices0.begin(), elem_indices0.end(),
                          elem_indices1.begin(), elem_indices1.end(),
                          std::inserter(neighbouring_element_indices, neighbouring_element_indices.begin()));

    return neighbouring_element_indices;
}

template<unsigned SPACE_DIM>
unsigned Edge<SPACE_DIM>::GetNumElements()
{
    return mElementIndices.size();
}

template<unsigned SPACE_DIM>
bool Edge<SPACE_DIM>::IsBoundaryEdge() const
{
    return mElementIndices.size()<=1;
}

template <unsigned SPACE_DIM>
bool Edge<SPACE_DIM>::operator==(const Edge<SPACE_DIM>& rEdge) const
{
    return this->ContainsNode(rEdge.GetNode(0)) && this->ContainsNode(rEdge.GetNode(1));
}

// Explicit instantiation
template class Edge<1>;
template class Edge<2>;
template class Edge<3>;