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

#include "BoundaryElement.hpp"

#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM, SPACE_DIM>::BoundaryElement()
    : AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>()
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM, SPACE_DIM>::BoundaryElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
    RegisterWithNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM, SPACE_DIM>::BoundaryElement(unsigned index, Node<SPACE_DIM>* pNode)
    : AbstractTetrahedralElement<ELEMENT_DIM,SPACE_DIM>(index)
{
    assert(ELEMENT_DIM == 0);

    // Store Node pointer
    this->mNodes.push_back(pNode);
    RegisterWithNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundaryElement<ELEMENT_DIM, SPACE_DIM>::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddBoundaryElement(this->mIndex);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundaryElement<ELEMENT_DIM, SPACE_DIM>::ResetIndex(unsigned index)
{
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveBoundaryElement(this->mIndex);
    }
    this->mIndex = index;
    RegisterWithNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundaryElement<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted()
{
    this->mIsDeleted = true;
//        this->mJacobianDeterminant = 0.0;
    // Update nodes in this element so they know they are not contained by us
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveBoundaryElement(this->mIndex);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundaryElement<ELEMENT_DIM, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[rIndex]->RemoveBoundaryElement(this->mIndex);

    // Update the node at this location
    this->mNodes[rIndex] = pNode;

    // Add element to this node
    this->mNodes[rIndex]->AddBoundaryElement(this->mIndex);
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class BoundaryElement<0,1>;
template class BoundaryElement<0,2>;
template class BoundaryElement<1,2>;
template class BoundaryElement<0,3>;
template class BoundaryElement<1,3>;
template class BoundaryElement<2,3>;
