/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "AbstractElement.hpp"

#include "Exception.hpp"

#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractElement<ELEMENT_DIM, SPACE_DIM>::AbstractElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : mNodes(rNodes),
      mIndex(index),
      mRegion(0.0),
      mIsDeleted(false),
      mOwnership(true),
      mFlag(false)
{
    // Sanity checking
    assert(ELEMENT_DIM <= SPACE_DIM);

    // Flags must be initialised before the Jacobian calculations, or assertions trip
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractElement<ELEMENT_DIM, SPACE_DIM>::AbstractElement(unsigned index)
    : mIndex(index),
      mRegion(0.0),
      mIsDeleted(false),
      mOwnership(true),
      mFlag(false)
{}

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
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::Flag()
{
    mFlag = true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::Unflag()
{
    mFlag = false;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractElement<ELEMENT_DIM, SPACE_DIM>::IsFlagged() const
{
    return mFlag;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractElement<ELEMENT_DIM, SPACE_DIM>::SetRegion(double region)
{
    mRegion = region;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractElement<ELEMENT_DIM, SPACE_DIM>::GetRegion()
{
    return mRegion;
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractElement<0,1>;
template class AbstractElement<1,1>;
template class AbstractElement<0,2>;
template class AbstractElement<1,2>;
template class AbstractElement<2,2>;
template class AbstractElement<0,3>;
template class AbstractElement<1,3>;
template class AbstractElement<2,3>;
template class AbstractElement<3,3>;
