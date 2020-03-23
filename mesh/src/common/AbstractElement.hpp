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

#ifndef ABSTRACTELEMENT_HPP_
#define ABSTRACTELEMENT_HPP_

#include <vector>

#include "UblasVectorInclude.hpp"
#include "Node.hpp"
#include "Edge.hpp"
#include "EdgeHelper.hpp"
#include "ElementAttributes.hpp"

template<unsigned int SPACE_DIM>
class EdgeHelper;
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

    /** The edges forming this element **/
    std::vector<Edge<SPACE_DIM>*> mEdges;

    /** EdgeHelper class to keep track of edges */
    EdgeHelper<SPACE_DIM>* mEdgeHelper;

    /** The index of this element within the mesh */
    unsigned mIndex;

    /**
     * Whether this element has been deleted, and hence
     * whether its location in the mesh can be re-used.
     */
    bool mIsDeleted;

    /** Whether the current process owns this element. */
    bool mOwnership;

    /** A pointer to an ElementAttributes object associated with this element. */
    ElementAttributes<ELEMENT_DIM, SPACE_DIM>* mpElementAttributes;

    /**
     * Construct an empty ElementAttributes container.
     */
    void ConstructElementAttributes();

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
     * Deletes added attributes (when they exist)
     */
    virtual ~AbstractElement();

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
     * @return a single component of the location in space of one of the nodes
     * in this element.
     *
     * @param localIndex  the index of the node to query, in [0,N) where N
     *   is the number of nodes in this element.
     * @param dimension  the spatial dimension to query.
     */
    double GetNodeLocation(unsigned localIndex, unsigned dimension) const;

    /**
     * @return the location in space of one of the nodes in this element.
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
     * @return the number of nodes owned by this element.
     */
    unsigned GetNumNodes() const;

    /**
     * Add a node to this element.
     *
     * @param pNode pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode);

    /**
     * TODO: Proper description
     * @return
     */
    bool CheckEdgesAreValid();


    /**
     * Get whether the element is marked as deleted.
     *
     * @return mIsDeleted
     */
    bool IsDeleted() const;

    /**
     *  @return the index of this element
     */
    unsigned GetIndex() const;

    /**
     * Set the index of this element in the mesh.
     *
     * @param index the new index
     */
    void SetIndex(unsigned index);

    /**
     * @return whether the current process owns this element.
     */
    bool GetOwnership() const;

    /**
     * Set whether the current process owns this element.
     *
     * @param ownership whether the current process now owns this element
     */
    void SetOwnership(bool ownership);

    /**
     * Set an attribute (a value associated with the element)
     *
     * @param attribute the value of an attribute
     */
    void SetAttribute(double attribute);

    /**
     * @return the element's attribute value
     */
    double GetAttribute();

    /**
     * This method converts the attribute (stored as a double) to an
     * unsigned, and throws an exception if this is not sensible.
     *
     * @return an unsigned attribute value.
     */
    unsigned GetUnsignedAttribute();

    /**
     * Add an attribute to the list of element attributes.
     *
     * @param attribute the value of the attribute to be added
     */
    void AddElementAttribute(double attribute);

    /**
     * @return reference to a vector containing the element attributes.
     */
    std::vector<double>& rGetElementAttributes();

    /**
     * @return the number of node attributes associated with this node.
     */
    unsigned GetNumElementAttributes();

    /**
     * Sets edge helper
     * @param edgeHelper
     */
    void SetEdgeHelper(EdgeHelper<SPACE_DIM>* edgeHelper);

    /**
     * Clear edges from element
     */
    void ClearEdges();

    /**
     * Builds edges from element nodes
     */
    void BuildEdges();

    /**
     * Gets the global index of the edge at localIndex
     * @param localIndex local index of the edge in this element
     * @return Global index of the edge
     */
    unsigned GetEdgeGlobalIndex(unsigned localIndex) const;

    /**
     * Gets the edge at localIndex
     * @param localIndex local index of the edge in this element
     * @return
     */
    Edge<SPACE_DIM>* GetEdge(unsigned localIndex) const;

    /**
     * @return Number of edges associated with this element
     */
    unsigned GetNumEdges() const;

    /**
     * Gets a set of element indices that neighours the element at the specified edge
     * @param localIndex Local index of the edge in this element
     * @return A set of element indices that neighbours this edge
     */
    std::set<unsigned> GetNeighbouringElementAtEdgeIndex(unsigned localIndex);

    /**
     * Checks if the element contains edge
     * @param edge
     * @return
     */
    bool ContainsEdge(const Edge<SPACE_DIM> *edge) const;

    /**
     * returns the local index of edge
     * @param edge
     * @return -1 if an edge was not found
     */
    long GetLocalEdgeIndex(const Edge<SPACE_DIM> *edge) const;
};

#endif /*ABSTRACTELEMENT_HPP_*/
