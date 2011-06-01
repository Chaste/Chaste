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
#ifndef POTTSELEMENT_HPP_
#define POTTSELEMENT_HPP_

#include "AbstractElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * An element class for use in the PottsMesh class.
 *
 * The main difference between this and the Element class is that a
 * PottsElement can have a variable number of nodes associated with
 * it and these represent the lattice sites contained in the element.
 */
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PottsElement : public AbstractElement<2, 2>
{
public:
    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     */
    PottsElement(unsigned index, const std::vector<Node<2>*>& rNodes);

    /**
     * Destructor.
     */
    ~PottsElement();

    /**
     * Overridden RegisterWithNodes() method.
     *
     * Informs all nodes forming this element that they are in this element.
     */
    void RegisterWithNodes();

    /**
     * Overridden MarkAsDeleted() method.
     *
     * Mark an element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed.
     */
    void MarkAsDeleted();

    /**
     * Reset the global index of the element and update its nodes.
     *
     * @param index the new global index
     */
    void ResetIndex(unsigned index);

    /**
     * Update node at the given index.
     *
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<2>* pNode);

    /**
     * Delete a node with given local index.
     *
     * @param rIndex is the local index of the node to remove
     */
    void DeleteNode(const unsigned& rIndex);

    /**
     * Add node to element. Note that we dont care about ordering in a potts
     * element so just add it to the end of the mNodes vector.
     *
     * @param pNode is a pointer to the new node
     */
    void AddNode(Node<2>* pNode);

    /**
     * Calculate the local index of a node given a global index
     * if node is not contained in element return UINT_MAX.
     *
     * @param globalIndex the global index of the node in the mesh
     * @return local_index.
     */
    unsigned GetNodeLocalIndex(unsigned globalIndex) const;

    /**
     * Get whether or not the element is on the boundary by seeing if contains boundary nodes.
     *
     * @return whether or not the element is on the boundary.
     */
    bool IsElementOnBoundary() const;
};

#endif /*POTTSELEMENT_HPP_*/
