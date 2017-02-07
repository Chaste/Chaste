/*

Copyright (c) 2005-2017, University of Oxford.
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
#ifndef VERTEXELEMENT_HPP_
#define VERTEXELEMENT_HPP_

#include "MutableElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * An element class for use in the VertexMesh class. The main
 * difference between this and the Element class is that a
 * VertexElement can have a variable number of nodes associated
 * with it.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement : public MutableElement<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * Faces of the VertexElement, which should be distinct.
     */
    std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> mFaces;

    /**
     * How each face is oriented. From the perspective of the centre
     * of the element, the vertices of each face should be ordered
     * anti clockwise. If and only if this is false, the order of vertices
     * in the corresponding face should be reversed.
     *
     * N.B. Most faces belong to two elements, but with opposite
     * orientations. This allows us to reuse the face data across the
     * two cells.
     */
    std::vector<bool> mOrientations;

    /**
     * Set of indices of elements containing this face.
     */
    std::set<unsigned> mFaceContainingElementIndices;

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
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        // This needs to be first so that MeshBasedCellPopulation::Validate() doesn't go mental.
        archive& mFaces;
        archive& mOrientations;
        archive& mFaceContainingElementIndices;
        archive& boost::serialization::base_object<MutableElement<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:
    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rFaces vector of faces associated with the element
     * @param rOrientations vector of orientations of the faces associated with the element
     */
    VertexElement(unsigned index,
                  const std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                  const std::vector<bool>& rOrientations);

    /**
     *
     * Alternative constructor.
     *
     * When constructing a VertexMesh as the Voronoi dual to a Delaunay mesh,
     * each VertexElement is initially constructed without nodes.
     *
     * @param index global index of the element
     */
    VertexElement(unsigned index);

    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     */
    VertexElement(unsigned index,
                  const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Constructor used to specify the element completely. This ensures that
     * the nodes and faces are owned by the element *in a specified order*.
     * See #1076 and #1377 for more details.
     *
     * @param index global index of the element
     * @param rFaces vector of faces associated with the element
     * @param rOrientations vector of orientations of the faces associated with the element
     * @param rNodes vector of Nodes associated with the element
     */
    VertexElement(unsigned index,
                  const std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                  const std::vector<bool>& rOrientations,
                  const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Destructor.
     */
    ~VertexElement();

    /**
     * Replace one of the nodes in this element with another. Override method in AbstractElement
     * to replace nodes of the faces as well.
     *
     * @param pOldNode  pointer to the current node
     * @param pNewNode  pointer to the replacement node
     */
    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode);

    /**
     * Compute the centroid of an element. Exact same function as VertexMesh::GetCentroidOfElement(const unsigned&)
     * but it has been added so that it is accessible for faces as well
     *
     * @return centroid of the VertexElement.
     */
    c_vector<double, SPACE_DIM> GetCentroid() const;

    /**
     * Overridden MarkAsDeleted() method.
     *
     * Mark an element as having been removed from the mesh.
     * Also notify nodes and faces in the element that it has been removed.
     */
    void MarkAsDeleted();

    /**
     * Reset the global index of the element and update its nodes and faces.
     *
     * @param newIndex the new global index
     */
    void ResetIndex(unsigned newIndex);

    /**
     * @return the number of faces owned by this element.
     */
    unsigned GetNumFaces() const;

    /**
     * @param index the local index of a specified face
     *
     * @return a pointer to the face
     */
    VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* GetFace(const unsigned localIndex) const;

    /** Calculate the local index of a face given a global index
     * if face is not contained in element return UINT_MAX.
     *
     * @param globalIndex the global index of the face
     * @return local_index.
     */
    unsigned GetFaceLocalIndex(const unsigned globalIndex) const;

    /**
     * @return whether the face with a given index is oriented anti clockwise.
     *
     * @param localIndex the local index of the face
     */
    bool FaceIsOrientatedAntiClockwise(const unsigned localIndex) const;

    /**
     * @param globalIndex  the global index of the face
     * @return  whether the face is oriented anti clockwise
     */
    bool GetFaceOrientationWithGlobalIndex(const unsigned globalIndex) const;

    /**
     * Add a face to the element.
     *
     * @param pFace a pointer to the new face
     */
    void AddFace(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace);

    /**
     * Add a face of an element with given local index. This method will
     * add the corresponding face orientation as well.
     *
     * @param pFace is the face to add
     * @param Orientation is the orientation of pFace
     * @param rIndex is the local index of the face to add
     */
    void AddFace(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, bool Orientation, const unsigned& rIndex);

    /**
     * Delete a face of an element with given local index. This method will
     * remove the corresponding face orientation as well.
     *
     * @param rIndex is the local index of the face to remove
     */
    void DeleteFace(const unsigned& rIndex);

    /**
     * Replace old face with the new face and its orientation
     */
    void ReplaceFace(const unsigned oldFaceLocalIndex, VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pNewFace,
                     const bool newFaceOrientation);

    /**
     * Method for monolayer element. Delete ApicalNode1 and BasalNode1.
     * This method will delete nodes from the element, apical & basal face,
     * update lateral face and return the index of the face which contains all
     * four nodes (does NOT mark it as deleted).
     *
     * @param pApicalNodeDelete is the pointer of apical node which will be removed.
     * @param pBasalNodeDelete is the pointer of basal node which will be removed.
     * @param pApicalNodeStay is the pointer of apical node which serve as reference to identify the lateral face.
     * @param pBasalNodeStay is the pointer of basal node similar to pApicalNode2.
     *
     * @return the global index of the face which contains all 4 nodes and its orientation
     */
    std::vector<unsigned> MonolayerElementDeleteNodes(const Node<SPACE_DIM>* pApicalNodeDelete, const Node<SPACE_DIM>* pBasalNodeDelete,
                                                      Node<SPACE_DIM>* pApicalNodeStay, Node<SPACE_DIM>* pBasalNodeStay);

    /**
     * Method for monolayer element to rearrange faces and nodes so that
     * other operations can be done more efficiently (e.g. no need to loop
     * through to find the next face of a node etc.)
     *
     * 'Correct order' is defined as such:
     * (1) basal face will be at the first place;
     * (2) apical face will come second;
     * (3) basal nodes come first before apical nodes, and in the same order as basal and apical face;
     * (4) lateral faces are adjusted by LateralFaceRearrangeNodes;
     * (5) lateral faces (after basal and apical faces) have the same order as the nodes, where
     *     lateral_faces[i] contains nodes[i] and nodes[i+1].
     * Assumption:  nodes within the faces are already in cyclic order;
     *              the nodes in apical and basal faces are synchronised.
     */
    void MonolayerElementRearrangeFacesNodes();

    /**
     * Check the orientation of individual face.
     * This function is called in MonolayerElementRearrangeFacesNodes
     * @param faceLocalIndex local index of face
     */
    void CheckFaceOrientationOfElement(const unsigned faceLocalIndex);

    /**
     * Method for lateral face of monolayer element to rearrange its nodes such that the
     * basal nodes are always the at the beginning.
     * This function is called in MonolayerElementRearrangeFacesNodes
     */
    void LateralFaceRearrangeNodes();

    /**
     * Reset the global index of the face and update its nodes.
     *
     * @param index the new global index
     */
    void FaceResetIndex(unsigned index);

    /**
     * Method for faces. Add a node to the element between nodes at Index and Index+1.
     *
     * @param Index the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     */
    void FaceAddNode(Node<SPACE_DIM>* pNode, const unsigned Index = UINT_MAX);

    /**
     * Method for faces. Delete a node with given local index.
     *
     * @param Index is the local index of the node to remove
     */
    void FaceDeleteNode(const unsigned Index);

    /**
     * Method for faces. Delete a node.
     *
     * @param rIndex is the local index of the node to remove
     */
    void FaceDeleteNode(const Node<SPACE_DIM>* pNode);

    /**
     * Method for faces. Update node at the given index.
     *
     * @param Index is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     */
    void FaceUpdateNode(const unsigned Index, Node<SPACE_DIM>* pNode);

    /**
     * Method for faces. Update oldNode with newNode.
     *
     * @param pOldNode is a pointer which will be changed
     * @param pNewNode is a pointer to the replacement node
     */
    void FaceUpdateNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode);

    /**
     * Informs all nodes forming this face that they are in this face.
     * Called normally together with FaceAddElement to do the registry together.
     */
    void RegisterFaceWithNodes();

    /**
     * Mark a face as having been removed from the mesh.
     * Doesn't notify nodes and element that it has been removed.
     */
    void MarkFaceAsDeleted();

    /**
     * For VertexElement<2, 2> only.
     * Arrange nodes of the faces so that they will be in correct order (CCW from top).
     * @return true if there are changes to the order of nodes.
     * This function does not update the element about the changes #2850
     */
    bool FaceRearrangeNodes();

    /**
     * For VertexElement<2, 3> only.
     * Arrange nodes of the faces so that they will be in correct order (CCW from PointOfView).
     * @param PointOfView is the position from which the face is observed.
     * @return true if there are changes to the order of nodes.
     * This function does not update the element about the changes #2850
     */
    bool FaceRearrangeNodes(const c_vector<double, SPACE_DIM> &PointOfView);

    /**
     * face will add element index into its registry.
     * @param elementIndex  the index of the element
     */
    void FaceAddElement(const unsigned elementIndex);

    /**
     * Remove a element that contains this face.
     *
     * @param elementIndex
     */
    void FaceRemoveElement(const unsigned elementIndex);

    /**
     *  @return the number of elements in the mesh that contain this face.
     */
    unsigned FaceGetNumContainingElements() const;

    /**
     * @return a set of indices of elements containing this face.
     */
    std::set<unsigned>& rFaceGetContainingElementIndices();
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
template <unsigned SPACE_DIM>
class VertexElement<1, SPACE_DIM> : public MutableElement<1, SPACE_DIM>
{

private:
    /**
     * Set of indices of elements containing this face.
     */
    std::set<unsigned> mFaceContainingElementIndices;

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
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        // This needs to be first so that MeshBasedCellPopulation::Validate() doesn't go mental.
        archive& mFaceContainingElementIndices;
        archive& boost::serialization::base_object<MutableElement<1, SPACE_DIM> >(*this);
    }

public:
    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes the nodes owned by the element
     */
    VertexElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * @return the number of faces owned by this element.
     */
    unsigned GetNumFaces() const;

    /**
     * @param index the global index of a specified face
     *
     * @return a pointer to the face
     */
    VertexElement<0, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * @param index the global index of a specified face
     *
     * @return the local index of the face
     */
    unsigned GetFaceLocalIndex(unsigned globalIndex) const;

    /**
     * @return whether the face with a given index is oriented anti clockwise.
     *
     * @param index the index of the face
     */
    bool FaceIsOrientatedAntiClockwise(const unsigned index) const;

    /**
     * Compute the centroid of an element. Exact same function as VertexMesh::GetCentroidOfElement(const unsigned&)
     * but it has been added so that it is accessible for faces as well
     *
     * @return centroid of the VertexElement.
     */
    c_vector<double, SPACE_DIM> GetCentroid() const;

    /**
     * Dummy function declarations to satisfy compiler.
     */
    void AddFace(VertexElement<0, SPACE_DIM>* pFace);
    void AddFace(VertexElement<0, SPACE_DIM>* pFace, bool Orientation, const unsigned& rIndex);
    void DeleteFace(const unsigned& rIndex);
    void ReplaceFace(const unsigned oldFaceLocalIndex, VertexElement<0, SPACE_DIM>* pNewFace,
                     const bool newFaceOrientation);
    std::vector<unsigned> MonolayerElementDeleteNodes(const Node<SPACE_DIM>* pApicalNodeDelete, const Node<SPACE_DIM>* pBasalNodeDelete,
                                                      Node<SPACE_DIM>* pApicalNodeStay, Node<SPACE_DIM>* pBasalNodeStay);
    void MonolayerElementRearrangeFacesNodes();
    void FaceResetIndex(unsigned index);
    void FaceAddNode(Node<SPACE_DIM>* pNode, const unsigned Index);
    void FaceDeleteNode(const unsigned Index);
    void FaceDeleteNode(const Node<SPACE_DIM>* pNode);
    void FaceUpdateNode(const unsigned Index, Node<SPACE_DIM>* pNode);
    void RegisterFaceWithNodes();
    void MarkFaceAsDeleted();
    void FaceAddElement(const unsigned elementIndex);
    void FaceRemoveElement(const unsigned elementIndex);
    unsigned FaceGetNumContainingElements() const;
    std::set<unsigned>& rFaceGetContainingElementIndices();
};

#endif /*VERTEXELEMENT_HPP_*/
