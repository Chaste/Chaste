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
#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>
#include <climits>


#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexElement.hpp"
#include "VertexElementMap.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestVertexMesh;

protected:

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangement. */
    double mCellRearrangementThreshold;

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
     *
     * @return the index of the new element
     */
    unsigned DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned nodeAIndex, unsigned nodeBIndex);

    /**
     * Test whether a given point lies inside a given element.
     *
     * We use a ray-casting algorithm, which relies on the following result:
     * if the point in question is not on the boundary of the element, then
     * the number of intersections is an even number if the point is outside,
     * and it is odd if inside.
     *
     * Currently the method is coded 'strictly', such that points lying on
     * an edge or at a vertex are considered to lie outside the element.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return if the point is included in the element.
     */
    bool ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /**
     * Get the local index of a given element which is the start vertex of the edge
     * of the element that the overlapping point rTestPoint is closest to.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return the local index
     */
    unsigned GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /**
     * Called by ReMesh(). Moves a node, which has been found to overlap an element,
     * back onto the edge of that element and associates it with the element.
     *
     * @param pNode pointer to the node
     * @param elementIndex global index of the element in the mesh
     */
    void MoveOverlappingNodeOntoEdgeOfElement(Node<SPACE_DIM>* pNode, unsigned elementIndex);

    /**
     * Solve node mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the VertexMesh and its member variables. Note that this will
     * write out a VertexMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mCellRearrangementThreshold;
        archive & mEdgeDivisionThreshold;
        archive & mT2Threshold;
        archive & mDeletedNodeIndices;
        archive & mDeletedElementIndices;

        // Create a mesh writer pointing to the correct file and directory
        VertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                             ArchiveLocationInfo::GetMeshFilename(),
                                                             false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<VertexMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
    }

    /**
     * Loads a mesh by using VertexMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mCellRearrangementThreshold;
        archive & mEdgeDivisionThreshold;
        archive & mT2Threshold;
        archive & mDeletedNodeIndices;
        archive & mDeletedElementIndices;

        VertexMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Forward declaration */
    class VertexElementIterator;

    /**
     * Get an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline VertexElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * Get an iterator to one past the last element in the mesh.
     */
    inline VertexElementIterator GetElementIteratorEnd();

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
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
               double cellRearrangementThreshold=0.01,
               double edgeDivisionThreshold=DBL_MAX,
               double t2Threshold=0.001);

    /**
     * Default constructor for use by serializer.
     */
    VertexMesh();

    /**
     * Destructor.
     */
    virtual ~VertexMesh();

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
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @param index  the global index of a specified vertex element
     *
     * @return a pointer to the vertex element
     */
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Compute the area of an element.
     *
     * This needs to be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the area of the element
     */
    virtual double GetAreaOfElement(unsigned index);

    /**
     * Compute the perimeter of an element.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the perimeter of the element
     */
    double GetPerimeterOfElement(unsigned index);

    /**
     * Compute the centroid of an element.
     *
     * This needs to be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Compute the area gradient of an element at one of its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the area of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of an element ending at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that ends at this node.
     */
    c_vector<double, SPACE_DIM> GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of an element starting at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that starts at this node.
     */
    c_vector<double, SPACE_DIM> GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the perimeter of an element at its nodes.
     * This returns the sum of GetPreviousEdgeGradientAtNode() and GetNextEdgeGradientAtNode().
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the perimeter of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the second moments of area of a given (polygonal) element.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfElement(unsigned index);

    /**
     * Calculate the vector of the shortest axis of a given element.
     * This is the eigenvector associated with the largest eigenvalue
     * of the inertial tensor. If the polygon is regular then the
     * eigenvalues are the same, so we return a random unit vector.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * \todo This method is only called inside DivideElementAlongShortAxis() -
     *       get rid of it and move the code into that method?
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (short_axis_x, short_axis_y).
     */
    c_vector<double, SPACE_DIM> GetShortAxisOfElement(unsigned index);

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Given a node and one of its containing elements, find a set containing
     * the indices of those neighbouring node(s) that are NOT also in the element.
     *
     * Note that we allow for more than one such index, since there is no reason
     * a priori to assume that each node is contained by exactly three elements.
     *
     * @param nodeIndex global index of the node
     * @param elemIndex global index of the element
     *
     * @return its neighbouring nodes that are not in the element
     */
    std::set<unsigned> GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);

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
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongShortAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement);

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
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongGivenAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, c_vector<double, SPACE_DIM> axisOfDivision);

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

    /**
     * Translate the mesh given the displacement vector.
     * This is the translation method that actually does the work.
     *
     * @param rDisplacement is a translation vector of the correct size
     */
    void Translate(c_vector<double, SPACE_DIM>& rDisplacement);

    /**
     * Translate the mesh given the coordinate displacements separately.
     *
     * @param xMovement is the x-displacement (defaults to 0.0)
     * @param yMovement is the y-displacement (defaults to 0.0)
     * @param zMovement is the z-displacement (defaults to 0.0)
     */
    void Translate(const double xMovement=0.0, const double yMovement=0.0, const double zMovement=0.0);

    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////

    /**
     * A smart iterator over the elements in the mesh.\todo This is the same as in AbstractTetrahedralMesh
     */
    class VertexElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         *
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline VertexElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * AbstractTetrahedralMesh::GetElementIteratorBegin and AbstractTetrahedralMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        VertexElementIterator(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
                        bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        VertexMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         */
        inline bool IsAllowedElement();
    };
};

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh);


//////////////////////////////////////////////////////////////////////////////
// VertexElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return VertexElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return VertexElementIterator(*this, mElements.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator!=(const VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::VertexElementIterator(
        VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
        bool skipDeletedElements)
    : mrMesh(rMesh),
      mElementIter(elementIter),
      mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mElements.size() == 0)
    {
        // Cope with empty meshes
        mElementIter = mrMesh.mElements.end();
    }
    else
    {
        // Make sure we start at an allowed element
        if (mElementIter == mrMesh.mElements.begin() && !IsAllowedElement())
        {
            ++(*this);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}


#endif /*VERTEXMESH_HPP_*/
