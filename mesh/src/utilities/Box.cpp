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
#include "Box.hpp"

template<unsigned DIM>
void Box<DIM>::AddNode(Node<DIM>* pNode)
{
    mNodesContained.insert(pNode);
}

template<unsigned DIM>
void Box<DIM>::RemoveNode(Node<DIM>* pNode)
{
    mNodesContained.erase(pNode);
}

template<unsigned DIM>
void Box<DIM>::ClearNodes()
{
    mNodesContained = std::set< Node<DIM>* >();
}

template<unsigned DIM>
std::set< Node<DIM>* >& Box<DIM>::rGetNodesContained()
{
    return mNodesContained;
}

template<unsigned DIM>
void Box<DIM>::AddElement(Element<DIM,DIM>* pElement)
{
    mElementsContained.insert(pElement);
}

///\todo #2308 there are no methods to remove or clear elements

template<unsigned DIM>
std::set< Element<DIM,DIM>* >& Box<DIM>::rGetElementsContained()
{
    return mElementsContained;
}

///////// Explicit instantiation///////

template class Box<1>;
template class Box<2>;
template class Box<3>;
