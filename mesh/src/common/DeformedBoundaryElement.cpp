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

#include "DeformedBoundaryElement.hpp"

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
DeformedBoundaryElement<ELEM_DIM,SPACE_DIM>::DeformedBoundaryElement()
    : BoundaryElement<ELEM_DIM,SPACE_DIM>()
{
    assert(ELEM_DIM>0 && ELEM_DIM<3);
    assert(ELEM_DIM+1==SPACE_DIM);

    c_vector<double,SPACE_DIM> x0 = zero_vector<double>(SPACE_DIM);
    this->AddNode(new Node<SPACE_DIM>(UINT_MAX, x0));

    c_vector<double,SPACE_DIM> x1 = zero_vector<double>(SPACE_DIM);
    x1(0) = 1.0;
    this->AddNode(new Node<SPACE_DIM>(UINT_MAX, x1));

    if (ELEM_DIM==2)
    {
        c_vector<double,SPACE_DIM> x2 = zero_vector<double>(SPACE_DIM);
        x2(1) = 1.0;
        this->AddNode(new Node<SPACE_DIM>(UINT_MAX, x2));
    }
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
DeformedBoundaryElement<ELEM_DIM,SPACE_DIM>::~DeformedBoundaryElement()
{
    for(unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
}

template<unsigned ELEM_DIM, unsigned SPACE_DIM>
void DeformedBoundaryElement<ELEM_DIM,SPACE_DIM>::ApplyUndeformedElementAndDisplacement(BoundaryElement<ELEM_DIM,SPACE_DIM>* pUndeformedElement,
                                                                                        std::vector<c_vector<double,SPACE_DIM> >& rDisplacement)
{
    for (unsigned i=0; i<NUM_NODES; i++)
    {
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            this->GetNode(i)->rGetModifiableLocation()[j]
              = pUndeformedElement->GetNode(i)->rGetLocation()[j] + rDisplacement[i](j);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class DeformedBoundaryElement<1,2>;
template class DeformedBoundaryElement<2,3>;
