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
#ifndef POTTSELEMENT_HPP_
#define POTTSELEMENT_HPP_

#include "AbstractElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * An element class for use in the VertexMesh class. The main
 * difference between this and the Element class is that a
 * VertexElement can have a variable number of nodes associated
 * with it.
 */
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PottsElement //: public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{
public:
    //** Default Constructor can remove once other constructors written */
    PottsElement();

//private:
//
//    /**
//     * Faces of the VertexElement, which should be distinct.
//     */
//    std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*> mFaces;
//
//    /**
//     * How each face is oriented. From the perspective of the centre
//     * of the element, the vertices of each face should be ordered
//     * anti clockwise. If and only if this is false, the order of vertices
//     * in the corresponding face should be reversed.
//     *
//     * N.B. Most faces belong to two elements, but with opposite
//     * orientations. This allows us to reuse the face data across the
//     * two cells.
//     */
//    std::vector<bool> mOrientations;
//
//    /** Needed for serialization. */
//    friend class boost::serialization::access;
//    /**
//     * Serialize the object and its member variables.
//     *
//     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
//     *
//     * Note also that member data related to writers is not saved - output must
//     * be set up again by the caller after a restart.
//     *
//     * @param archive the archive
//     * @param version the current version of this class
//     */
//    template<class Archive>
//    void serialize(Archive & archive, const unsigned int version)
//    {
//        // This needs to be first so that MeshBasedCellPopulation::Validate() doesn't go mental.
//     //   archive & mFaces;
//     //   archive & mOrientations;
//     //   archive & boost::serialization::base_object<AbstractElement<ELEMENT_DIM, SPACE_DIM> >(*this);
//    }
//
//public:
//
//    /**
//     * Constructor.
//     *
//     * @param index global index of the element
//     * @param rFaces vector of faces associated with the element
//     * @param rOrientations vector of orientations of the faces associated with the element
//     */
//    VertexElement(unsigned index,
//                  const std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*>& rFaces,
//                  const std::vector<bool>& rOrientations);
//
//    /**
//     *
//     * Alternative constructor.
//     *
//     * When constructing a VertexMesh as the Voronoi dual to a Delaunay mesh,
//     * each VertexElement is initially constructed without nodes.
//     *
//     * @param index global index of the element
//     */
//    VertexElement(unsigned index);
//
//    /**
//     * Constructor.
//     *
//     * @param index global index of the element
//     * @param rNodes vector of Nodes associated with the element
//     */
//    VertexElement(unsigned index,
//                  const std::vector<Node<SPACE_DIM>*>& rNodes);
//
//    /**
//     * Constructor used to specify the element completely. This ensures that
//     * the nodes and faces are owned by the element *in a specified order*.
//     * See #1076 and #1377 for more details.
//     *
//     * @param index global index of the element
//     * @param rFaces vector of faces associated with the element
//     * @param rOrientations vector of orientations of the faces associated with the element
//     * @param rNodes vector of Nodes associated with the element
//     */
//    VertexElement(unsigned index,
//                  const std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*>& rFaces,
//                  const std::vector<bool>& rOrientations,
//                  const std::vector<Node<SPACE_DIM>*>& rNodes);
//
//    /**
//     * Destructor.
//     */
//    ~VertexElement();
//
//    /**
//     * Get the number of faces owned by this element.
//     */
//    unsigned GetNumFaces() const;
//
//    /**
//     * Overridden RegisterWithNodes() method.
//     *
//     * Informs all nodes forming this element that they are in this element.
//     */
//    void RegisterWithNodes();
//
//    /**
//     * Overridden MarkAsDeleted() method.
//     *
//     * Mark an element as having been removed from the mesh.
//     * Also notify nodes in the element that it has been removed.
//     */
//    void MarkAsDeleted();
//
//    /**
//     * Reset the global index of the element and update its nodes.
//     *
//     * @param index the new global index
//     */
//    void ResetIndex(unsigned index);
//
//    /**
//     * Update node at the given index.
//     *
//     * @param rIndex is an local index to which node to change
//     * @param pNode is a pointer to the replacement node
//     */
//    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);
//
//    /**
//     * Delete a node with given local index.
//     *
//     * @param rIndex is the local index of the node to remove
//     */
//    void DeleteNode(const unsigned& rIndex);
//
//    /**
//     * Add a node to the element between nodes at rIndex and rIndex+1.
//     *
//     * @param rIndex the local index of the node after which the new node is added
//     * @param pNode a pointer to the new node
//     */
//    void AddNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);
//
//    /**
//     * Add a face to the element.
//     *
//     * @param pFace a pointer to the new face
//     */
//    void AddFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace);
//
//    /**
//     * Calculate the local index of a node given a global index
//     * if node is not contained in element return UINT_MAX.
//     *
//     * @param globalIndex the global index of the node in the mesh
//     * @return local_index.
//     */
//    unsigned GetNodeLocalIndex(unsigned globalIndex) const;
//
//    /**
//     * @param index the global index of a specified face
//     *
//     * @return a pointer to the face
//     */
//    VertexElement<ELEMENT_DIM-1, SPACE_DIM>* GetFace(unsigned index) const;
//
//    /**
//     * Get whether the face with a given index is oriented clockwise.
//     *
//     * @param index the index of the face
//     */
//    bool FaceIsOrientatedClockwise(unsigned index) const;
//
//    /**
//     * Get whether or not the element is on the boundary by seeing if contains boundary nodes.
//     *
//     * @return whether or not the element is on the boundary.
//     */
//    bool IsElementOnBoundary() const;
//
//};
//
//
////////////////////////////////////////////////////////////////////////
////                  Specialization for 1d elements                  //
////                                                                  //
////                 1d elements are just edges (lines)               //
////////////////////////////////////////////////////////////////////////
//
///**
// * Specialization for 1d elements so we don't get errors from Boost on some
// * compilers.
// */
//template<unsigned SPACE_DIM>
//class VertexElement<1, SPACE_DIM> : public AbstractElement<1,SPACE_DIM>
//{
//public:
//
//    /**
//     * Constructor which takes in a vector of nodes.
//     *
//     * @param index  the index of the element in the mesh
//     * @param rNodes the nodes owned by the element
//     */
//    VertexElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);
//
//    /**
//     * Virtual destructor, since this class has virtual methods.
//     */
//    virtual ~VertexElement();
//
//    /**
//     * Get the number of faces owned by this element.
//     */
//    unsigned GetNumFaces() const;
//
//    /**
//     * Update node at the given index.
//     *
//     * @param rIndex is an local index to which node to change
//     * @param pNode is a pointer to the replacement node
//     */
//    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);
//
//    /**
//     * Overridden RegisterWithNodes() method.
//     *
//     * Informs all nodes forming this element that they are in this element.
//     */
//    void RegisterWithNodes();
//
//    /**
//     * Overridden MarkAsDeleted() method.
//     *
//     * Mark an element as having been removed from the mesh.
//     * Also notify nodes in the element that it has been removed.
//     */
//    void MarkAsDeleted();
//
//    /**
//     * Reset the global index of the element and update its nodes.
//     *
//     * @param index the new global index
//     */
//    void ResetIndex(unsigned index);
//
//    /**
//     * Delete a node with given local index.
//     *
//     * @param rIndex is the local index of the node to remove
//     */
//    void DeleteNode(const unsigned& rIndex);
//
//    /**
//     * Add a node to the element between nodes at rIndex and rIndex+1.
//     *
//     * @param rIndex the local index of the node after which the new node is added
//     * @param pNode a pointer to the new node
//     */
//    void AddNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode);
//
//    /**
//     * Calculate the local index of a node given a global index
//     * if node is not contained in element return UINT_MAX
//     *
//     * @param globalIndex the global index of the node in the mesh
//     * @return local_index.
//     */
//    unsigned GetNodeLocalIndex(unsigned globalIndex) const;
//
//    /**
//     * @param index the global index of a specified face
//     *
//     * @return a pointer to the face
//     */
//    VertexElement<0, SPACE_DIM>* GetFace(unsigned index) const;
//
//    /**
//     * Get whether the face with a given index is oriented clockwise.
//     *
//     * @param index the index of the face
//     */
//    bool FaceIsOrientatedClockwise(unsigned index) const;
//
//    /**
//     * Get whether or not the element is on the boundary by seeing if contains boundary nodes.
//     *
//     * @return whether or not the element is on the boundary.
//     */
//    bool IsElementOnBoundary() const;
};

#endif /*POTTSELEMENT_HPP_*/
