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

#include "Edge.hpp"

template<unsigned SPACE_DIM>
Edge<SPACE_DIM>::Edge(unsigned index)
    : mIndex(index),
      mIsDeleted(false)
{
    this->mIndex = index;
}

template<unsigned SPACE_DIM>
Edge<SPACE_DIM>::Edge(unsigned index, Node<SPACE_DIM>* pNode0, Node<SPACE_DIM>* pNode1)
    : mIndex(index),
      mIsDeleted(false)
{
    this->SetNodes(pNode0, pNode1);
}

template<unsigned SPACE_DIM>
Edge<SPACE_DIM>::~Edge()
{
    mNodes.clear();
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::MarkDeleted()
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
unsigned Edge<SPACE_DIM>::GetIndex()
{
    return mIndex;
}

template<unsigned SPACE_DIM>
UIndexPair Edge<SPACE_DIM>::GetMapIndex()
{
    assert(mNodes.size() == 2);
    auto index0 = mNodes[0]->GetIndex();
    auto index1 = mNodes[1]->GetIndex();
    if (index0 > index1)
    {
        auto indexSwap = index0;
        index0 = index1;
        index1 = indexSwap;
    }

    return UIndexPair(index0, index1);
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::RemoveNodes()
{
    mNodes.clear();
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::SetNodes(Node<SPACE_DIM>* pNode0, Node<SPACE_DIM>* pNode1)
{
    // Clear the nodes first
    this->RemoveNodes();

    // Add nodes
    this->mNodes.push_back(pNode0);
    this->mNodes.push_back(pNode1);
}

template<unsigned SPACE_DIM>
void Edge<SPACE_DIM>::ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
{
    for (unsigned i = 0; i < mNodes.size(); i++)
    {
        if (this->mNodes[i] == pOldNode)
        {
            this->mNodes[i] = pNewNode;
            break;
        }
    }
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>* Edge<SPACE_DIM>::GetNode(unsigned index)
{
    return mNodes[index];
}

template<unsigned SPACE_DIM>
unsigned Edge<SPACE_DIM>::GetNumNodes()
{
    return mNodes.size();
}

template<unsigned SPACE_DIM>
bool Edge<SPACE_DIM>::ContainsNode(Node<SPACE_DIM>* pNode)
{
    for (auto node : mNodes)
    {
        if (node->GetIndex() == pNode->GetIndex())
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
    std::set<unsigned> currentElem;
    currentElem.insert(elementIndex);
    std::set_difference(currentElem.begin(), currentElem.end(),
            mElementIndices.begin(), mElementIndices.end(),
            std::inserter(otherElements, otherElements.begin()));

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
bool Edge<SPACE_DIM>::IsEdgeValid()
{
    // MUST have 2 existing nodes to form an edge
    if (mNodes.size() != 2)
    {
        printf("[Error] Edge %i - has less than two nodes\n", this->mIndex);
        return false;
    }

    // Nodes should not be nullptr
    for (auto node : mNodes)
    {
        if (node == nullptr)
        {
            printf("[Error] Edge %i - has a nullptr node\n", this->mIndex);
            return false;
        }
    }

    // Can't have associated elements if we're less than 2D
    if (SPACE_DIM <= 1 && mElementIndices.size() > 0)
    {
        printf("[Error] Edge %i - Can't have an associated element if less than 2D\n", this->mIndex);
        return false;
    }

    // An edge can only have a maximum of two elements in 2D
    if (SPACE_DIM == 2 && mElementIndices.size() > 2)
    {
        printf("[Error] Edge %i - an edge can only have a maximum of two neighbouring elements in 2D\n", this->mIndex);
        return false;
    }

    auto neighbour_indices = GetNeighbouringElementIndices();
    if (neighbour_indices != mElementIndices)
    {
        printf("[Error] Edge %i - the neighbouring elements in mElementIndices does not match with node\n", this->mIndex);
        return false;
    }

    return true;
}

// Explicit instantiation
template class Edge<1>;
template class Edge<2>;
template class Edge<3>;
