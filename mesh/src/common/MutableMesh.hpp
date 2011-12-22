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

#ifndef MUTABLEMESH_HPP_
#define MUTABLEMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

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
     * Serialize the mesh.
     *
     * Note that if you are calling this method (from subclasses) you should archive your
     * member variables FIRST. So that this method can call a ReMesh
     * (to convert from TrianglesMeshReader input format into your native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        // Do a remesh after archiving has finished to get right number of boundary nodes etc.
        // (strictly only relevant for load - but doesn't take long in the scheme of things.)
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        NodeMap map(this->GetNumNodes());
        this->ReMesh(map);
        assert(map.IsIdentityMap()); // Otherwise the mesh will get VERY confused.
    }

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

#define COVERAGE_IGNORE
    /**
     * Check whether any neighbouring node is inside the circumsphere of this element.
     *
     * @param pElement pointer to an element
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of the element, as a proportion of the circumsphere radius.
     */
    bool CheckIsVoronoi(Element<ELEMENT_DIM, SPACE_DIM>* pElement, double maxPenetration);
#undef COVERAGE_IGNORE

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
     * Get the number of nodes that are actually in use.
     */
    unsigned GetNumNodes() const;

    /**
     * Get the number of elements that are actually in use.
     */
    unsigned GetNumElements() const;

    /**
     * Get the number of boundary elements that are actually in use.
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
     */
    virtual unsigned AddNode(Node<SPACE_DIM>* pNewNode);

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

#define COVERAGE_IGNORE
    /**
     * Delete a node from the mesh by finding an appropriate neighbour node
     * to merge it with.
     *
     * @param index is the index of the node to be deleted
     */
    virtual void DeleteNode(unsigned index);
#undef COVERAGE_IGNORE

#define COVERAGE_IGNORE
    /**
     * Mark a node as deleted. Note that this method DOES NOT deal with the
     * associated elements and therefore should only be called immediately prior
     * to a ReMesh() being called. (Thus saves work compared to DeleteNode()
     * function and does not MoveMerge the node and elements).
     *
     * @param index The index of the node to delete
     */
    void DeleteNodePriorToReMesh(unsigned index);
#undef COVERAGE_IGNORE

    /**
     * Refine an element at a given point.
     *
     * @param pElement  pointer to the element
     * @param point  a point located in the element
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

#define COVERAGE_IGNORE
    /**
     * Re-index a mesh so that it has no deleted elements or nodes.
     *
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    void ReIndex(NodeMap& map);
#undef COVERAGE_IGNORE

#define COVERAGE_IGNORE
    /**
     * Re-mesh a mesh using triangle (via library calls) or tetgen
     * @param map is a NodeMap which associates the indices of nodes in the old mesh
     * with indices of nodes in the new mesh.  This should be created with the correct size (NumAllNodes)
     */
    virtual void ReMesh(NodeMap& map);
#undef COVERAGE_IGNORE


#define COVERAGE_IGNORE
    /**
     * Alternative version of remesh which takes no parameters does not require a NodeMap. Note: inherited
     * classes should overload ReMesh(NodeMap&)
     */
    void ReMesh();
#undef COVERAGE_IGNORE

#define COVERAGE_IGNORE
    /**
     * Checks the entire mesh element by element and checks whether any neighbouring node
     * is inside the circumsphere of this element.
     *
     * @param maxPenetration is the maximum distance a node is allowed to be inside the
     * circumsphere of an element that it is not a member of, as a proportion of the
     * circumsphere radius.
     */
    bool CheckIsVoronoi(double maxPenetration=0.0);
#undef COVERAGE_IGNORE
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableMesh)

#endif /*MUTABLEMESH_HPP_*/
