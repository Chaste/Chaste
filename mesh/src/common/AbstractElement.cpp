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

#include "AbstractElement.hpp"

#include "MathsCustomFunctions.hpp"
#include "Exception.hpp"

#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractElement<ELEMENT_DIM, SPACE_DIM>::AbstractElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : mNodes(rNodes),
      mEdgeHelper(nullptr),
      mIndex(index),
      mIsDeleted(false),
      mOwnership(true),
      mpElementAttributes(nullptr)

{
    // Sanity checking
    assert(ELEMENT_DIM <= SPACE_DIM);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractElement<ELEMENT_DIM, SPACE_DIM>::AbstractElement(unsigned index)
    : mEdgeHelper(nullptr),
      mIndex(index),
      mIsDeleted(false),
      mOwnership(true),
      mpElementAttributes(nullptr)

{}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractElement<ELEMENT_DIM, SPACE_DIM>::~AbstractElement()
{

    delete mpElementAttributes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
{
    assert(pOldNode != pNewNode);
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->mNodes[i] == pOldNode)
        {
            UpdateNode(i, pNewNode);
            return;
        }
    }
    EXCEPTION("You didn't have that node to start with.");

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNodeLocation(unsigned localIndex, unsigned dimension) const
{
    assert(dimension < SPACE_DIM);
    assert((unsigned)localIndex < mNodes.size());
    return mNodes[localIndex]->rGetLocation()[dimension];
}

/*
 * Note for future reference: this used to return a reference to a c_vector, in which case a
 * weird error arose where it compiled, ran and passed on some machines but failed the tests
 * (bad_size errors) on another machine.  So be careful if you think about changing it!
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNodeLocation(unsigned localIndex) const
{
    assert((unsigned)localIndex < mNodes.size());
    return mNodes[localIndex]->rGetLocation();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNodeGlobalIndex(unsigned localIndex) const
{
    assert((unsigned)localIndex < mNodes.size());
    return mNodes[localIndex]->GetIndex();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned localIndex) const
{
    assert((unsigned)localIndex < mNodes.size());
    return mNodes[localIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode)
{
    mNodes.push_back(pNode);
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::SetEdgeHelper(EdgeHelper<SPACE_DIM> *edgeHelper)
{
    this->mEdgeHelper = edgeHelper;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::ClearEdges()
{
    for (auto edge: mEdges)
    {
        edge->RemoveElement(this->mIndex);
    }
    mEdges.clear();
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::BuildEdges()
{
    assert(mEdgeHelper != nullptr);

    // If SPACE_DIM == 2 then we can assume that the node layout
    // in the array corresponds to its connections
    // We can then infer the edge connection information
    if (SPACE_DIM == 2)
    {
        this->ClearEdges();
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            unsigned i_next = (i+1) % mNodes.size();
            mEdges.push_back(mEdgeHelper->GetEdgeFromNodes(this->mIndex, mNodes[i], mNodes[i_next]));
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned int SPACE_DIM>
unsigned AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetEdgeGlobalIndex(unsigned localIndex) const
{
    assert(localIndex < mEdges.size());
    return mEdges[localIndex]->GetIndex();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Edge<SPACE_DIM> *AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetEdge(unsigned localIndex) const
{
    assert(localIndex < mEdges.size());
    return mEdges[localIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNumEdges() const {
    return mEdges.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringElementAtEdgeIndex(unsigned localIndex)
{
    assert(localIndex < mEdges.size());
    return mEdges[localIndex]->GetOtherElements(this->mIndex);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractElement<ELEMENT_DIM, SPACE_DIM>::CheckEdgesAreValid()
{
    if (SPACE_DIM != 2)
    {
        return false;
    }

    if (mEdges.size() != mNodes.size())
    {
        return false;
    }

    bool edgesOk = true;
    for (unsigned i = 0; i < mEdges.size(); i++)
    {
        auto i_next = (i+1) % mEdges.size();
        auto edge = mEdges[i];

        if (!edge->IsEdgeValid())
        {
            printf("[Error] Element index: %i edge index: %i - edge is not valid \n", this->mIndex, i);
            edgesOk = false;
        }

        // Edge at index i must contain nodes at index i and i+1
        if (!(edge->ContainsNode(mNodes[i]) && edge->ContainsNode(mNodes[i_next])))
        {
            printf("[Error] Element index: %i edge index: %i - doesn't contain the correct nodes \n", this->mIndex, i);
            edgesOk = false;
        }
    }

    if (!edgesOk)
    {
        printf("[Debug] Element index: %i Node Indices: ", this->mIndex);
        for (auto nodeOut: this->mNodes)
        {
            printf(" %i ", nodeOut->GetIndex());
        }
        printf("\n Edge-node indices");
        for (auto edgeOut : this->mEdges)
        {
            printf(" %i (%i -> %i)", edgeOut->GetIndex(),
                    edgeOut->GetNode(0)->GetIndex(), edgeOut->GetNode(1)->GetIndex());
        }
        printf("\n");

        return false;
    }

    return true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractElement<ELEMENT_DIM, SPACE_DIM>::IsDeleted() const
{
    return mIsDeleted;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetIndex() const
{
    return mIndex;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::SetIndex(unsigned index)
{
    mIndex = index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetOwnership() const
{
    return mOwnership;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::SetOwnership(bool ownership)
{
    mOwnership = ownership;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::SetAttribute(double attribute)
{
    ConstructElementAttributes();

    mpElementAttributes->SetFirstAttribute(attribute);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetAttribute()
{
    ConstructElementAttributes();

    return mpElementAttributes->GetFirstAttribute();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetUnsignedAttribute()
{
    double double_attr = GetAttribute();
    unsigned unsigned_attr = (unsigned) (double_attr + 0.5);

    if (CompareDoubles::WithinAnyTolerance(double_attr, unsigned_attr) == false)
    {
        EXCEPTION("Element attribute '"<< double_attr <<"' cannot be converted to an unsigned.");
    }
    return unsigned_attr;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::ConstructElementAttributes()
{
    if (mpElementAttributes == nullptr)
    {
        mpElementAttributes = new ElementAttributes<ELEMENT_DIM, SPACE_DIM>();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::AddElementAttribute(double attribute)
{
    ConstructElementAttributes();

    mpElementAttributes->AddAttribute(attribute);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double>& AbstractElement<ELEMENT_DIM, SPACE_DIM>::rGetElementAttributes()
{
    if (mpElementAttributes == nullptr)
    {
        EXCEPTION("Element has no attributes associated with it. Construct attributes first");
    }

    return mpElementAttributes->rGetAttributes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes()
{
    return mpElementAttributes == nullptr ? 0 : mpElementAttributes->rGetAttributes().size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractElement<ELEMENT_DIM, SPACE_DIM>::ContainsEdge(const Edge<SPACE_DIM> *edge) const
{
    for (unsigned int i=0; i<mEdges.size(); ++i)
    {
        if ((*mEdges[i])==(*edge))
            return true;
    }
    return false;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
long AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetLocalEdgeIndex(const Edge<SPACE_DIM> *edge) const
{
    long result = -1;
    for (unsigned int i=0; i<mEdges.size(); ++i)
    {
        if ((*mEdges[i])==(*edge))
        {
            result = i;
        }
    }
    return result;
}


// Explicit instantiation
template class AbstractElement<0,1>;
template class AbstractElement<1,1>;
template class AbstractElement<0,2>;
template class AbstractElement<1,2>;
template class AbstractElement<2,2>;
template class AbstractElement<0,3>;
template class AbstractElement<1,3>;
template class AbstractElement<2,3>;
template class AbstractElement<3,3>;
