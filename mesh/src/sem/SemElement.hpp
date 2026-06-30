/*

Copyright (c) 2005-2026, University of Oxford.
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
#ifndef SEMELEMENT_HPP_
#define SEMELEMENT_HPP_

#include "AbstractElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * An element class for use in the SemMesh class.
 *
 * An SemElement represents the subcellular node collection for one
 * biological cell. It records which mesh nodes belong to that cell and the
 * node membership bookkeeping needed for population lifecycle operations.
 * Each node carries a region label (0 = interior, 1 = boundary) set by the
 * mesh generator; forces such as SemRegionalForce read this via Node::GetRegion().
 */
template<unsigned DIM>
class SemElement : public AbstractElement<DIM, DIM>
{
private:

    /**
     * The id of the biological cell represented by this SEM element.
     */
    unsigned mCellId;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractElement<DIM, DIM> >(*this);
    }

public:

    /**
     * Default constructor, which doesn't add any nodes: they must be added later.
     *
     * This creates an empty SEM element shell with a mesh element index. Nodes
     * can then be added and registered once the element geometry is known.
     *
     * @param index global index of the element
     */
    SemElement(unsigned index);

    /**
     * Constructor which takes in a vector of nodes.
     *
     * This creates a complete SEM element around an existing set of
     * subcellular nodes and registers the element index with those nodes so
     * later force, damping, and deletion queries can find containing elements.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     */
    SemElement(unsigned index,
               const std::vector<Node<DIM>*>& rNodes);
    

    /**
     * Destructor.
     *
     * The element does not own its nodes; node memory is managed by SemMesh.
     */
    ~SemElement();
    
    /**
     * Set the id of the biological cell represented by this element.
     *
     * This records an external cell identifier on the SEM element for force,
     * output, or bookkeeping code that needs to associate element state with a
     * cell id.
     * 
     * @param id the new cell id
     */
    void SetCellId(unsigned int id);
    
    /**
     * Get the element's node vector.
     *
     * This provides direct access to the subcellular nodes owned by the SEM
     * element for code that needs to inspect or update the element geometry as
     * a node collection.
     *
     * @return reference to the vector of node pointers
     */
    std::vector<Node<DIM>*>& rGetNodes();

    /**
     * Replace one node in the element.
     *
     * This keeps node-to-element membership consistent while swapping the node
     * pointer at a local index: the old node is unregistered from this element
     * and the new node is registered.
     *
     * @param rIndex local node index to update
     * @param pNode replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<DIM>* pNode) override;

    /**
     * Mark this element as deleted.
     *
     * This removes the SEM element from active population use and unregisters
     * it from all of its nodes so later containing-element queries do not
     * include deleted cells.
     */
    void MarkAsDeleted() override;

    /**
     * Register this element with its nodes.
     *
     * This populates each node's containing-element set with this element
     * index, enabling node-based damping, neighbour, force, and
     * lifecycle queries to recover the SEM elements that contain a node.
     */
    void RegisterWithNodes() override;
};

#endif /*SEMELEMENT_HPP_*/
