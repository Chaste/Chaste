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
#ifndef MUTABLEELEMENT_HPP_
#define MUTABLEELEMENT_HPP_

#include "AbstractElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 *  A mutable element containing functionality
 *  to add and remove nodes.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableElement : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{
private:


    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractElement<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     *
     * Alternative constructor.
     *
     * @param index global index of the element
     */
    MutableElement(unsigned index);

    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     */
    MutableElement(unsigned index,
                  const std::vector<Node<SPACE_DIM>*>& rNodes);


    /**
     * Destructor.
     */
    virtual ~MutableElement();

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
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);

    /**
     * Delete a node with given local index.
     *
     * @param rIndex is the local index of the node to remove
     */
    void DeleteNode(const unsigned& rIndex);

    /**
     * Add a node to the element between nodes at rIndex and rIndex+1.
     *
     * @param rIndex the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex);

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
    virtual bool IsElementOnBoundary() const;
};

//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
class MutableElement<1, SPACE_DIM> : public AbstractElement<1,SPACE_DIM>
{
public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes the nodes owned by the element
     */
    MutableElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     *
     * Alternative constructor.
     *
     * @param index global index of the element
     */
    MutableElement(unsigned index);

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~MutableElement();

    /**
     * Update node at the given index.
     *
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);

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
     * Delete a node with given local index.
     *
     * @param rIndex is the local index of the node to remove
     */
    void DeleteNode(const unsigned& rIndex);

    /**
     * Add a node to the element between nodes at rIndex and rIndex+1.
     *
     * @param rIndex the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex);

    /**
     * Calculate the local index of a node given a global index
     * if node is not contained in element return UINT_MAX
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
    virtual bool IsElementOnBoundary() const;
};

#endif /*MUTABLEELEMENT_HPP_*/
