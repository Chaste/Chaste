/*

Copyright (C) University of Oxford, 2005-2012

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

#include <cassert>

#include "Node.hpp"
#include "Exception.hpp"

//////////////////////////////////////////////////////////////////////////
// Constructors
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::CommonConstructor(unsigned index, bool isBoundaryNode)
{
    mIndex = index;
    mIsBoundaryNode = isBoundaryNode;
    mIsInternal = false;
    mIsDeleted = false;
    mRegion = 0;
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, ChastePoint<SPACE_DIM> point, bool isBoundaryNode)
{
    mLocation = point.rGetLocation();
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, std::vector<double> coords, bool isBoundaryNode)
{
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        mLocation(i) = coords.at(i);
    }
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, c_vector<double, SPACE_DIM> location, bool isBoundaryNode)
{
    mLocation = location;
    CommonConstructor(index, isBoundaryNode);
}

template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, bool isBoundaryNode, double v1, double v2, double v3)
{
    mLocation[0] = v1;
    if (SPACE_DIM > 1)
    {
        mLocation[1] = v2;
        if (SPACE_DIM > 2)
        {
            mLocation[2] = v3;
        }
    }
    CommonConstructor(index, isBoundaryNode);
}
template<unsigned SPACE_DIM>
Node<SPACE_DIM>::Node(unsigned index, double *location, bool isBoundaryNode)
{
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        mLocation(i) = location[i];
    }
    CommonConstructor(index, isBoundaryNode);

}
//////////////////////////////////////////////////////////////////////////
// Methods dealing with node location
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetPoint(ChastePoint<SPACE_DIM> point)
{
    mLocation = point.rGetLocation();
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetIndex(unsigned index)
{
    mIndex = index;
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetAsBoundaryNode(bool value)
{
    mIsBoundaryNode = value;
}


template<unsigned SPACE_DIM>
ChastePoint<SPACE_DIM> Node<SPACE_DIM>::GetPoint() const
{
    return ChastePoint<SPACE_DIM>(mLocation);
}

template<unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& Node<SPACE_DIM>::rGetLocation() const
{
    assert(!mIsDeleted);
    return mLocation;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM>& Node<SPACE_DIM>::rGetModifiableLocation()
{
    assert(!mIsDeleted);
    return mLocation;
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetIndex() const
{
    return mIndex;
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsBoundaryNode() const
{
    return mIsBoundaryNode;
}


template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddNodeAttribute(double attribute)
{
    mNodeAttributes.push_back(attribute);
}

template<unsigned SPACE_DIM>
std::vector<double>& Node<SPACE_DIM>::rGetNodeAttributes()
{
    return mNodeAttributes;
}

//////////////////////////////////////////////////////////////////////////
// Tracking (boundary) elements which contain this node as a vertex
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddElement(unsigned index)
{
    mElementIndices.insert(index);
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::RemoveElement(unsigned index)
{
    unsigned count = mElementIndices.erase(index);
    if (count == 0)
    {
        EXCEPTION("Tried to remove an index which was not in the set");
    }
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::RemoveBoundaryElement(unsigned index)
{
    unsigned count = mBoundaryElementIndices.erase(index);
    if (count == 0)
    {
        EXCEPTION("Tried to remove an index which was not in the set");
    }
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::AddBoundaryElement(unsigned index)
{
    mBoundaryElementIndices.insert(index);
}

template<unsigned SPACE_DIM>
std::set<unsigned>& Node<SPACE_DIM>::rGetContainingElementIndices()
{
    return mElementIndices;
}

template<unsigned SPACE_DIM>
std::set<unsigned>& Node<SPACE_DIM>::rGetContainingBoundaryElementIndices()
{
    return mBoundaryElementIndices;
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetNumContainingElements() const
{
    return mElementIndices.size();
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetNumBoundaryElements() const
{
    return mBoundaryElementIndices.size();
}

//////////////////////////////////////////////////////////////////////////
// Methods dealing with some node flags (deleted, region)
//////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::MarkAsDeleted()
{
    mIsDeleted = true;
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsDeleted() const
{
    return mIsDeleted;
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::MarkAsInternal()
{
    mIsInternal = true;
}

template<unsigned SPACE_DIM>
bool Node<SPACE_DIM>::IsInternal() const
{
    return mIsInternal;
}

template<unsigned SPACE_DIM>
void Node<SPACE_DIM>::SetRegion(unsigned region)
{
    mRegion = region;
}

template<unsigned SPACE_DIM>
unsigned Node<SPACE_DIM>::GetRegion() const
{
    return mRegion;
}


//////////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////////

template class Node<1>;
template class Node<2>;
template class Node<3>;
