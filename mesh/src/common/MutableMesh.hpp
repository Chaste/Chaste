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

#ifndef MUTABLEMESH_HPP_
#define MUTABLEMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "TetrahedralMesh.hpp"
#include "NodeMap.hpp"

/**
 * A concrete mutable mesh class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableMesh : public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save the mesh, along with attributes if nodes have attributes.
     * Saving of attributes is covered in TestNodesOnlyMesh.
     * @param archive the archive to save to.
     * @param version the version number.
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        // Assume that the first node is indicative of the rest.
        bool does_have_attributes = this->mNodes[0]->HasNodeAttributes();

        archive & does_have_attributes;

        if (does_have_attributes)
        {
            for (unsigned i=0; i<this->mNodes.size(); i++)
            {
                double radius;
                radius = this->mNodes[i]->GetRadius();
                archive & radius;

                bool is_particle;
                is_particle = this->mNodes[i]->IsParticle();
                archive & is_particle;
            }
        }
    }

    /**
     * Load the mesh, along with attributes if saved nodes have attributes.
     * Loading of attributes is covered in TestNodesOnlyMesh.
     * @param archive the archive to save to.
     * @param version the version number.
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        bool does_have_attributes;

        archive & does_have_attributes;

        if (does_have_attributes)
        {
            for (unsigned i=0; i<this->mNodes.size(); i++)
            {
                double radius;
                archive & radius;
                this->mNodes[i]->SetRadius(radius);

                bool is_particle;
                archive & is_particle;

                if (is_particle)
                {
                    this->mNodes[i]->SetIsParticle(true);
                }
            }
        }

        // If ELEMENT_DIM==SPACE_DIM do a remesh after archiving has finished to get right number of boundary nodes etc.
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        if (ELEMENT_DIM == SPACE_DIM)
        {
            NodeMap map(this->GetNumNodes());
            this->ReMesh(map);
            assert(map.IsIdentityMap()); // Otherwise the mesh will get VERY confused.
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected:

    /**
     * Indices of elements that have been marked as deleted.
     * These indices can be reused when adding new elements.
     */
    std::vector<unsigned> mDeletedElementIndices;

    /**
     * Indices of boundary elements that have been marked as deleted.
     * These indices can be reused when adding new boundary elements.
     */
    std::vector<unsigned> mDeletedBoundaryElementIndices;

    /**
     * Indices of nodes that have been marked as deleted.
     * These indices can be reused when adding new nodes.
     */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Whether any nodes have been added to the mesh. */
    bool mAddedNodes;

private:

    /**
     * @return true if the mesh is Voronoi local to the given element.
     * Check whether any neighbouring node is inside the circumsphere of this element.
     *
     * @param pElement pointer to an element
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of the element, as a proportion of the circumsphere radius.
     */
    bool CheckIsVoronoi(Element<ELEMENT_DIM, SPACE_DIM>* pElement, double maxPenetration);


public:

    /**
     * Constructor.
     */
    MutableMesh();

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param nodes  a vector of nodes
     */
    MutableMesh(std::vector<Node<SPACE_DIM> *> nodes);

    /**
     * Destructor.
     */
    virtual ~MutableMesh();

    /**
     * Clear all the data in the mesh.
     */
    void Clear();

    /**
     * @return the number of nodes that are actually in use.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of elements that are actually in use.
     */
    unsigned GetNumElements() const;

    /**
     * @return the number of boundary elements that are actually in use.
     */
    unsigned GetNumBoundaryElements() const;

    /// \todo should unsigned GetNumBoundaryNodes() be overloaded too??

    /**
     * Rescale the mesh from a boundary node.
     *
     * @param updatedPoint point determining the scale factor
     * @param boundaryNodeIndex index of the boundary node
     */
    void RescaleMeshFromBoundaryNode(ChastePoint<1> updatedPoint, unsigned boundaryNodeIndex);

    /**
     * Add a node to the mesh.
     *
     * NB. After calling this one or more times, you must then call ReMesh
     *
     * @param pNewNode  pointer to the new node
     * @return the index of the new node in the mesh
     */
    virtual unsigned AddNode(Node<SPACE_DIM>* pNewNode);

    /**
     * Helper method to add an element to the mesh.
     *
     * @param pNewElement pointer to the new element object. The object will be updated with the index assigned to the element.
     * @return new element index
     */
    unsigned AddElement(Element<ELEMENT_DIM,SPACE_DIM>* pNewElement);

    /**
     * Move the node with a particular index to a new point in space and
     * verifies that the signed areas of the supporting Elements are positive.
     *
     * @param index is the index of the node to be moved
     * @param point is the new target location of the node
     * @param concreteMove is set to false if we want to skip the signed area tests (defaults to true)
     */
    virtual void SetNode(unsigned index, ChastePoint<SPACE_DIM> point, bool concreteMove=true);

    /**
     * Move one node to another (i.e. merges the nodes), refreshing/deleting
     * elements as appropriate.
     *
     * @param index is the index of the node to be moved
     * @param targetIndex is the index of the node to move to
     * @param concreteMove can be set to false if you just want to check whether this will work (defaults to true).
     *     Set it to true if you're doing the merger for real, in order to do all the bookkeeping.
     */
    void MoveMergeNode(unsigned index, unsigned targetIndex, bool concreteMove=true);

    /**
     * Delete a node from the mesh by finding an appropriate neighbour node
     * to merge it with.
     *
     * @param index is the index of the node to be deleted
     */
    virtual void DeleteNode(unsigned index);

    /**
     * Delete an element from the mesh. Any remaining nodes not connected to an element are also removed.
     *
     * @param index The index of the element to be deleted
     */
    virtual void DeleteElement(unsigned index);

    /**
     * Mark a node as deleted. Note that this method DOES NOT deal with the
     * associated elements and therefore should only be called immediately prior
     * to a ReMesh() being called. (Thus saves work compared to DeleteNode()
     * function and does not MoveMerge the node and elements).
     *
     * @param index The index of the node to delete
     */
    void DeleteNodePriorToReMesh(unsigned index);

    /**
     * Refine an element at a given point.
     *
     * @param pElement  pointer to the element
     * @param point  a point located in the element
     * @return index of the new node which has been created at the given location
     */
    unsigned RefineElement(Element<ELEMENT_DIM,SPACE_DIM>* pElement, ChastePoint<SPACE_DIM> point);

    /**
     * Remove a boundary node, and update all the appropriate data structures.
     *
     * The deleted node is not removed from the list, merely marked as deleted,
     * and can be reused when a new node is added to the mesh.
     *
     * Any elements or boundary elements containing this node will be removed.
     * The boundary nodes information will be updated with new boundary node(s).
     * NB: New boundary elements WILL NOT be added.
     *
     * @param index  The index of the node to remove.
     */
    void DeleteBoundaryNodeAt(unsigned index);

    /**
     * Re-index a mesh so that it has no deleted elements or nodes.
     *
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    void ReIndex(NodeMap& map);


    /**
     * Re-mesh a mesh using triangle (via library calls) or tetgen
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    virtual void ReMesh(NodeMap& map);

    /**
     * Alternative version of remesh which takes no parameters, i.e. does not require a NodeMap.
     * It will create one and call the other ReMesh method.
     * Note: inherited classes should overload ReMesh(NodeMap&)
     */
    void ReMesh();


    /**
     * Find edges in the mesh longer than the given cutoff length and split them creating new elements as required.
     * @param cutoffLength cutoff length for edge splitting
     * @return returns a vector of triples with pointers to the new node followed by pointers to the nodes defining the bisected edge.
     */
    std::vector<c_vector<unsigned, 5> > SplitLongEdges(double cutoffLength);

    /**
     * Splits an edge in two and create all the relevant new elements and nodes
     *
     * @param pNodeA first point defining the edge
     * @param pNodeB second point defining the edge
     * @return the index of the new node
     */
    c_vector<unsigned, 3> SplitEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * @return true if Voronoi.  Checks the entire mesh element by element and checks whether any neighbouring node
     * is inside the circumsphere of this element.
     *
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of an element that it is not a member of, as a proportion of the
     * circumsphere radius.
     */
    bool CheckIsVoronoi(double maxPenetration=0.0);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableMesh)

#endif /*MUTABLEMESH_HPP_*/
