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


#ifndef _BOUNDARYELEMENT_HPP_
#define _BOUNDARYELEMENT_HPP_

#include <vector>

#include "AbstractTetrahedralElement.hpp"
#include "Node.hpp"

/**
 * Concrete boundary element class which inherits from AbstractTetrahedralElement.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BoundaryElement : public AbstractTetrahedralElement<ELEMENT_DIM, SPACE_DIM>
{
protected:
    /**
     *  (Protected) constructor that does take in nodes. Only available
     *  to subclasses. Calling code should use one of the other constructors.
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
