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
#include "MutableElement.hpp"
#include "RandomNumberGenerator.hpp"
#include <cassert>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableElement<ELEMENT_DIM, SPACE_DIM>::MutableElement(unsigned index)
    : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableElement<ELEMENT_DIM, SPACE_DIM>::MutableElement(unsigned index,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
    : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
    if (SPACE_DIM == ELEMENT_DIM)
    {
        RegisterWithNodes();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableElement<ELEMENT_DIM, SPACE_DIM>::~MutableElement()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableElement<ELEMENT_DIM, SPACE_DIM>::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddElement(this->mIndex);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableElement<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;

    // Update nodes in the element so they know they are not contained by it
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveElement(this->mIndex);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableElement<ELEMENT_DIM, SPACE_DIM>::ResetIndex(unsigned index)
{
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
       this->mNodes[i]->RemoveElement(this->mIndex);
    }
    this->mIndex = index;
    RegisterWithNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableElement<ELEMENT_DIM, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Update the node at this location
    this->mNodes[rIndex] = pNode;

    // Add element to this node
    this->mNodes[rIndex]->AddElement(this->mIndex);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableElement<ELEMENT_DIM, SPACE_DIM>::DeleteNode(const unsigned& rIndex)
{
    assert(rIndex < this->mNodes.size());

    // Remove element from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Remove the node at rIndex (removes node from element)
    this->mNodes.erase(this->mNodes.begin() + rIndex);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableElement<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex)
{
    /**
     * When constructing a VertexMesh as the Voronoi dual to a Delaunay mesh,
     * each MutableElement is initially constructed without nodes. We therefore
     * require the two cases below.
     */
    if (this->mNodes.empty())
    {
        // Populate mNodes with pNode
        this->mNodes.push_back(pNode);

        // Add element to this node
        this->mNodes[0]->AddElement(this->mIndex);
    }
    else
    {
        assert(rIndex < this->mNodes.size());

        // Add pNode to rIndex+1 element of mNodes pushing the others up
        this->mNodes.insert(this->mNodes.begin() + rIndex+1,  pNode);

        // Add element to this node
        this->mNodes[rIndex+1]->AddElement(this->mIndex);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableElement<ELEMENT_DIM, SPACE_DIM>::GetNodeLocalIndex(unsigned globalIndex) const
{
    unsigned local_index = UINT_MAX;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->GetNodeGlobalIndex(i) == globalIndex)
        {
            local_index = i;
        }
    }
    return local_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableElement<ELEMENT_DIM, SPACE_DIM>::IsElementOnBoundary() const
{
    bool is_element_on_boundary = false;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->GetNode(i)->IsBoundaryNode())
        {
            is_element_on_boundary = true;
            break;
        }
    }
    return is_element_on_boundary;
}

//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
MutableElement<1, SPACE_DIM>::MutableElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : AbstractElement<1, SPACE_DIM>(index, rNodes)
{
    // Sanity checking
    assert(this->mNodes.size() == 2);
    assert(SPACE_DIM > 0);
}

template<unsigned SPACE_DIM>
MutableElement<1, SPACE_DIM>::~MutableElement()
{
}

template<unsigned SPACE_DIM>
void MutableElement<1, SPACE_DIM>::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddElement(this->mIndex);
    }
}

template<unsigned SPACE_DIM>
void MutableElement<1, SPACE_DIM>::MarkAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;

    // Update nodes in the element so they know they are not contained by it
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveElement(this->mIndex);
    }
}

template <unsigned SPACE_DIM>
void MutableElement<1, SPACE_DIM>::ResetIndex(unsigned index)
{
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
       this->mNodes[i]->RemoveElement(this->mIndex);
    }
    this->mIndex = index;
    RegisterWithNodes();
}

template<unsigned SPACE_DIM>
void MutableElement<1, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Update the node at this location
    this->mNodes[rIndex] = pNode;

    // Add element to this node
    this->mNodes[rIndex]->AddElement(this->mIndex);
}

template<unsigned SPACE_DIM>
void MutableElement<1, SPACE_DIM>::DeleteNode(const unsigned& rIndex)
{
    assert(rIndex < this->mNodes.size());

    // Remove element from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Remove the node at rIndex (removes node from element)
    this->mNodes.erase(this->mNodes.begin() + rIndex);
}

template<unsigned SPACE_DIM>
void MutableElement<1, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex)
{
    assert(rIndex < this->mNodes.size());

    // Add pNode to rIndex+1 element of mNodes pushing the others up
    this->mNodes.insert(this->mNodes.begin() + rIndex+1,  pNode);

    // Add element to this node
    this->mNodes[rIndex+1]->AddElement(this->mIndex);
}

template<unsigned SPACE_DIM>
unsigned MutableElement<1, SPACE_DIM>::GetNodeLocalIndex(unsigned globalIndex) const
{
    unsigned local_index = UINT_MAX;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->GetNodeGlobalIndex(i) == globalIndex)
        {
            local_index = i;
        }
    }
    return local_index;
}

template<unsigned SPACE_DIM>
bool MutableElement<1, SPACE_DIM>::IsElementOnBoundary() const
{
    bool is_element_on_boundary = false;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->GetNode(i)->IsBoundaryNode())
        {
            is_element_on_boundary = true;
            break;
        }
    }
    return is_element_on_boundary;
}

// Explicit instantiation
template class MutableElement<1,1>;
template class MutableElement<1,2>;
template class MutableElement<1,3>;
template class MutableElement<2,2>;
template class MutableElement<2,3>;
template class MutableElement<3,3>;
