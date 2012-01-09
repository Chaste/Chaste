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

#ifndef ABSTRACTELEMENT_HPP_
#define ABSTRACTELEMENT_HPP_

#include <vector>

#include "UblasVectorInclude.hpp"
#include "Node.hpp"

/*
 * When creating an element within a mesh one needs to specify its global index.
 * If the element is not used within a mesh the following constant is used instead.
 */
const unsigned INDEX_IS_NOT_USED=0;

/**
 * An abstract element class for use in finite element meshes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractElement
{
protected:

    /** The nodes forming this element. */
    std::vector<Node<SPACE_DIM>*> mNodes;

    /** The index of this element within the mesh */
    unsigned mIndex;

    /** A region ID. */
    double mRegion;

    /**
     * Whether this element has been deleted, and hence
     * whether its location in the mesh can be re-used.
     */
    bool mIsDeleted;

    /** Whether the current process owns this element. */
    bool mOwnership;

    /** A flag for the use of higher level algorithms. */
    bool mFlag;

public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes  the nodes owned by the element
     */
    AbstractElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Default constructor, which doesn't add any nodes: they must be added later.
     *
     * @param index  the index of the element in the mesh (defaults to INDEX_IS_NOT_USED)
     */
    AbstractElement(unsigned index=INDEX_IS_NOT_USED);

    /**
     * Virtual destructor, since this class has virtual methods.
     * Does nothing special.
     */
    virtual ~AbstractElement()
    {}

    /**
     * Update node at the given index.
     *
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     */
    virtual void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)=0;

    /**
     * Replace one of the nodes in this element with another.
     *
     * @param pOldNode  pointer to the current node
     * @param pNewNode  pointer to the replacement node
     */
    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode);

    /**
     * Mark the element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed.
     */
    virtual void MarkAsDeleted()=0;

    /**
     * Inform all nodes forming this element that they are in this element.
     */
    virtual void RegisterWithNodes()=0;

    /**
     * Get a single component of the location in space of one of the nodes
     * in this element.
     *
     * @param localIndex  the index of the node to query, in [0,N) where N
     *   is the number of nodes in this element.
     * @param dimension  the spatial dimension to query.
     */
    double GetNodeLocation(unsigned localIndex, unsigned dimension) const;

    /**
     * Get the location in space of one of the nodes in this element.
     *
     * @param localIndex  the index of the node to query, in [0,N) where N
     *   is the number of nodes in this element.
     */
    c_vector<double, SPACE_DIM> GetNodeLocation(unsigned localIndex) const;

    /**
     * Given the local index of a node owned by this element, return the
     * global index of the node in the mesh.
     *
     * @param localIndex the node's local index in this element
     * @return the global index
     */
    unsigned GetNodeGlobalIndex(unsigned localIndex) const;

    /**
     * Get the node with a given local index in this element.
     *
     * @param localIndex local index of the node in this element
     * @return a pointer to the node.
     */
    Node<SPACE_DIM>* GetNode(unsigned localIndex) const;

    /**
     * Get the number of nodes owned by this element.
     */
    unsigned GetNumNodes() const;

    /**
     * Add a node to this element.
     *
     * @param pNode pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode);

    /**
     * Get whether the element is marked as deleted.
     *
     * @return mIsDeleted
     */
    bool IsDeleted() const;

    /**
     *  Get the index of this element
     */
    unsigned GetIndex() const;

    /**
     * Set the index of this element in the mesh.
     *
     * @param index the new index
     */
    void SetIndex(unsigned index);

    /**
     * Get whether the current process owns this element.
     */
    bool GetOwnership() const;

    /**
     * Set whether the current process owns this element.
     *
     * @param ownership whether the current process now owns this element
     */
    void SetOwnership(bool ownership);

    /**
     * Mark the element as flagged.
     */
    void Flag();

    /**
     * Mark the element as not flagged.
     */
    void Unflag();

    /**
     * Get whether the element is flagged.
     */
    bool IsFlagged() const;

    /**
     * Set the element's region ID.
     *
     * @param region the element's new region ID
     */
    void SetRegion(double region);

    /**
     * Get the element's region ID.
     */
    double GetRegion();
};

#endif /*ABSTRACTELEMENT_HPP_*/
