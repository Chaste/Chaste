/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef MUTABLEVERTEXMESH_HPP_
#define MUTABLEVERTEXMESH_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include <climits> // Work around Boost bug
#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "VertexMesh.hpp"
#include "ArchiveLocationInfo.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableVertexMesh : public VertexMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestMutableVertexMesh;
    friend class TestMutableVertexMeshReMesh;

protected:

    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangement. */
    double mCellRearrangementThreshold;

    /** The ratio between the minimum distance apart that two nodes in the mesh can be without causing element
     *  rearrangement and their seperation after remeshing.
     */
    double mCellRearrangementRatio;

    /** The maximum distance apart that neighbouring nodes in the mesh can be without the edge being divided. */
    double mEdgeDivisionThreshold;

    /** The area threshold at which T2 swaps occur in an apoptotic, triangular cell/element */
    double mT2Threshold;

    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Indices of elements that have been deleted. These indices can be reused when adding new elements. */
    std::vector<unsigned> mDeletedElementIndices;

    /**
     * Helper method for ReMesh to Identify the type of swap
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pNodeA one of the nodes to perform the swap with
     * @param pNodeB the other node to perform the swap
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    void IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, VertexElementMap& rElementMap);

    /**
     * Helper method for ReMesh to merge nodes when needed.
     * Replaces the node contained in the least number of elements with the other node.
     *
     * @param pNodeA one of the nodes to perform the merge with
     * @param pNodeB the other node to perform the merge with
     */
    void PerformNodeMerge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh to perform the T1 Swap
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pNodeA one of the nodes to perform the swap with
     * @param pNodeB the other node to perform the swap
     * @param rElementsContainingNodes set of common elements
     */
    void PerformT1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, std::set<unsigned>& rElementsContainingNodes);

    /**
     * Divide an element along the axis passing through two of its nodes.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement the element to divide
     * @param nodeAIndex the local index of one node within this element
     * @param nodeBIndex the local index of another node within this element
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                           unsigned nodeAIndex,
                           unsigned nodeBIndex,
                           bool placeOriginalElementBelow=false);

    /**
     * Called by ReMesh(). Moves a node, which has been found to overlap an element,
     * back onto the edge of that element and associates it with the element and adds
     * new nodes to maitain three elemets a node.
     *
     * @param pNode pointer to the node
     * @param elementIndex global index of the element in the mesh
     */
    void PerformT3Swap(Node<SPACE_DIM>* pNode, unsigned elementIndex);

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
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        archive & mCellRearrangementThreshold;
        archive & mCellRearrangementRatio;
        archive & mEdgeDivisionThreshold;
        archive & mT2Threshold;
        archive & mDeletedNodeIndices;
        archive & mDeletedElementIndices;

        archive & boost::serialization::base_object<VertexMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment (defaults to 0.01)
     * @param edgeDivisionThreshold the maximum threshold distance for edge division (defaults to DBL_MAX)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     */
    MutableVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
               double cellRearrangementThreshold=0.01,
               double edgeDivisionThreshold=DBL_MAX,
               double t2Threshold=0.001);

    /**
     * Default constructor for use by serializer.
     */
    MutableVertexMesh();

    /**
     * Destructor.
     */
    virtual ~MutableVertexMesh();

    /**
     * Set method for mCellRearrangementThreshold.
     *
     * @param cellRearrangementThreshold
     */
    void SetCellRearrangementThreshold(double cellRearrangementThreshold);

    /**
     * Set method for mEdgeDivisionThreshold.
     *
     * @param edgeDivisionThreshold
     */
    void SetEdgeDivisionThreshold(double edgeDivisionThreshold);

    /**
     * Set method for mT2Threshold.
     *
     * @param t2Threshold
     */
    void SetT2Threshold(double t2Threshold);

    /**
     *  Move the node with a particular index to a new point in space.
     *
      * @param nodeIndex the index of the node to be moved
      * @param point the new target location of the node
      */
    virtual void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);

    /**
     * @return mCellRearrangementThreshold
     */
    double GetCellRearrangementThreshold() const;

    /**
     * @return mEdgeDivisionThreshold
     */
    double GetEdgeDivisionThreshold() const;

    /**
     * @return mT2Threshold
     */
    double GetT2Threshold() const;

    /**
     * @return the number of Nodes in the mesh.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    unsigned GetNumElements() const;

    /**
     * Add a node to the mesh.
     *
     * Note: After calling this one or more times, you must then call ReMesh.
     *
     * @param pNewNode pointer to the new node
     * @return the global index of the new node in the mesh.
     */
    unsigned AddNode(Node<SPACE_DIM>* pNewNode);

    /**
     * Helper method for ReMesh to perform the T2 Swap
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param rElement the element to remove
     */
    void PerformT2Swap(VertexElement<ELEMENT_DIM,SPACE_DIM>& rElement);

    /**
     * Mark an element as deleted. Note that it DOES NOT deal with the associated
     * nodes and therefore should only be called immediately prior to a ReMesh()
     * being called.
     *
     * @param index  the global index of a specified vertex element
     */
    void DeleteElementPriorToReMesh(unsigned index);

    /**
     * Divide an element along its short axis.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement the element to divide
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongShortAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                         bool placeOriginalElementBelow=false);

    /**
     * Divide an element along a specified axis.
     *
     * If the new nodes (intersections of axis with element) are within
     * mCellRearrangementThreshold of existing nodes then they are
     * moved 2*mCellRearrangementThreshold away.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement the element to divide
     * @param axisOfDivision axis to divide the element by
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongGivenAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                         c_vector<double, SPACE_DIM> axisOfDivision,
                                         bool placeOriginalElementBelow=false);

    /**
     * Add an element to the mesh.
     *
     * @param pNewElement the new element
     *
     * @return the index of the new element in the mesh
     */
    unsigned AddElement(VertexElement<ELEMENT_DIM, SPACE_DIM>* pNewElement);

    /**
     * Delete mNodes and mElements.
     */
    void Clear();

    /**
     * Add a node on the edge between two nodes.
     *
     * @param pNodeA a pointer to one node
     * @param pNodeB a pointer to the other nodes
     */
    void DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh(). Removes the deleted nodes and elements from the mesh and updates the
     * rElementMap accordingly.
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    void RemoveDeletedNodesAndElements(VertexElementMap& rElementMap);

    /**
     * Helper method for ReMesh(). Removes the deleted nodes from the mesh and relabels the node indices.
     */
    void RemoveDeletedNodes();

    /**
     * Re-mesh the mesh.
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    void ReMesh(VertexElementMap& rElementMap);

    /**
     * Alternative version of remesh which takes no parameters does not require a VertexElementMap.
     * Note: inherited classes should overload ReMesh(VertexElementMap&).
     */
    void ReMesh();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableVertexMesh);

#endif /*MUTABLEVERTEXMESH_HPP_*/
