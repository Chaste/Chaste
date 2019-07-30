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


#ifndef _BOUNDARYELEMENT_HPP_
#define _BOUNDARYELEMENT_HPP_

#include <vector>

#include "AbstractTetrahedralElement.hpp"
#include "Node.hpp"

/**
 * Concrete boundary element class which inherits from AbstractTetrahedralElement.
 *
 * A 'face' in Chaste is shorthand for BoundaryElement<2,3>
 * i.e. a 2D boundary surface element on the edge of a 3D mesh.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BoundaryElement : public AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>
{
protected:
    /**
     * (Protected) constructor that does take in nodes. Only available
     * to subclasses. Calling code should use one of the other constructors.
     */
    BoundaryElement();

public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes  the nodes owned by the element
     */
    BoundaryElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Create a new boundary element from a Node.
     *
     * The element has ELEMENT_DIM=0 and SPACE_DIM identical to
     * that of the node from which it is constructed.
     *
     * @param index  the index of the element in the mesh
     * @param pNode  a pointer to the node
     */
    BoundaryElement(unsigned index, Node<SPACE_DIM>* pNode);

    /**
     * Inform all nodes forming this element that they are in this element.
     */
    void RegisterWithNodes();

    /**
     * Reset the index of this boundary element in the mesh.
     *
     * @param index the new index of the boundary element
     */
    void ResetIndex(unsigned index);

    /**
     * Mark the element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed.
     */
    void MarkAsDeleted();

    /**
     * Update node at the given index.
     *
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);
};

#endif //_BOUNDARYELEMENT_HPP_
