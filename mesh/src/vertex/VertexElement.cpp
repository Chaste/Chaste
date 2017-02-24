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
#include "VertexElement.hpp"

#include <algorithm> // std::swap, std::min, std::max, std::rotate, std::copy
#include "UblasCustomFunctions.hpp" // VectorProduct
#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "Debug.hpp"

//////////////////////////////////////////////////
///       Constructors and destructors         ///
//////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
          mFaces(rFaces),
          mOrientations(rOrientations)
{
    // This constructor should be used on <2,3> and <3,3> mesh
    // <2,3> mesh is used currently for geodesic sphere
    assert(ELEMENT_DIM >= 2 && SPACE_DIM == 3); // LCOV_EXCL_LINE

    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    // Register element&nodes with faces
    assert(mFaceContainingElementIndices.size() == 0);
    for (unsigned face_index = 0; face_index < mFaces.size(); ++face_index)
    {
        mFaces[face_index]->FaceAddElement(index);
        mFaces[face_index]->RegisterFaceWithNodes();
    }

    // Register element with nodes
    this->RegisterWithNodes();

    // A simple sanity check
    if (ELEMENT_DIM == 3)
    {
        unsigned edge_2x = 0;
        for (unsigned face_index = 0; face_index < mFaces.size(); ++face_index)
        {
            edge_2x += mFaces[face_index]->GetNumNodes();
        }
        assert(edge_2x % 2 == 0); // LCOV_EXCL_LINE
        assert(mFaces.size() + this->GetNumNodes() - edge_2x / 2 == 2); // LCOV_EXCL_LINE
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index),
          mFaces(rFaces),
          mOrientations(rOrientations)
{
    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    // Register element&nodes with faces
    assert(mFaceContainingElementIndices.size() == 0);
    for (unsigned face_index = 0; face_index < mFaces.size(); ++face_index)
    {
        mFaces[face_index]->FaceAddElement(index);
        mFaces[face_index]->RegisterFaceWithNodes();
    }

    // Make a set of nodes with mFaces
    std::set<Node<SPACE_DIM>*> nodes_set;
    for (unsigned face_index = 0; face_index < mFaces.size(); ++face_index)
    {
        for (unsigned node_index = 0; node_index < mFaces[face_index]->GetNumNodes(); node_index++)
        {
            nodes_set.insert(mFaces[face_index]->GetNode(node_index));
        }
    }

    // Populate mNodes
    for (typename std::set<Node<SPACE_DIM>*>::iterator node_iter = nodes_set.begin();
         node_iter != nodes_set.end();
         ++node_iter)
    {
        this->mNodes.push_back(*node_iter);
    }

    // Register element with nodes
    this->RegisterWithNodes();

    // A simple sanity check
    if (ELEMENT_DIM == 3)
    {
        unsigned edge_2x = 0;
        for (unsigned face_index = 0; face_index < mFaces.size(); ++face_index)
        {
            edge_2x += mFaces[face_index]->GetNumNodes();
        }
        assert(edge_2x % 2 == 0); // LCOV_EXCL_LINE
        assert(mFaces.size() + this->GetNumNodes() - edge_2x / 2 == 2); // LCOV_EXCL_LINE
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::~VertexElement()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
{
    // first call the function implementated in abstract Element to rnt
    AbstractElement<ELEMENT_DIM, SPACE_DIM>::ReplaceNode(pOldNode, pNewNode);

    for (unsigned face_index = 0; face_index < mFaces.size(); ++face_index)
    {
        VertexElement<ELEMENT_DIM - 1, SPACE_DIM>& r_face = *(mFaces[face_index]);
        const unsigned node_local_index = r_face.GetNodeLocalIndex(pOldNode->GetIndex());
        if (node_local_index != UINT_MAX)
            r_face.FaceUpdateNode(node_local_index, pNewNode);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<ELEMENT_DIM, SPACE_DIM>::GetCentroid() const
{
    unsigned num_nodes = this->GetNumNodes();

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

    switch (SPACE_DIM)
    {
        case 2:
        {
            double centroid_x = 0;
            double centroid_y = 0;

            // Note that we cannot use GetVolumeOfElement() below as it returns the absolute, rather than signed, area
            double element_signed_area = 0.0;

            // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
            c_vector<double, SPACE_DIM> first_node_location = this->GetNodeLocation(0);
            c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

            // Loop over vertices
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                c_vector<double, SPACE_DIM> next_node_location = this->GetNodeLocation((local_index + 1) % num_nodes);
                c_vector<double, SPACE_DIM> pos_2 = next_node_location - first_node_location;

                double this_x = pos_1[0];
                double this_y = pos_1[1];
                double next_x = pos_2[0];
                double next_y = pos_2[1];

                double signed_area_term = this_x * next_y - this_y * next_x;

                centroid_x += (this_x + next_x) * signed_area_term;
                centroid_y += (this_y + next_y) * signed_area_term;
                element_signed_area += 0.5 * signed_area_term;

                pos_1 = pos_2;
            }

            assert(element_signed_area != 0.0);

            // Finally, map back and employ GetVectorFromAtoB() to allow for periodicity
            centroid = first_node_location;
            centroid(0) += centroid_x / (6.0 * element_signed_area);
            centroid(1) += centroid_y / (6.0 * element_signed_area);
        }
        break;
        case 3:
        {
            ///\todo compute centroid rather than centre of mass (see #1422)
            for (unsigned local_index = 0; local_index < num_nodes; ++local_index)
            {
                centroid += this->GetNodeLocation(local_index);
            }
            centroid /= ((double)num_nodes);
        }
        break;
        default:
            NEVER_REACHED;
    }
    return centroid;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted()
{
    this->MutableElement<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted();

    // Update faces in the element so they know they are not contained by it
    for (unsigned i = 0; i < this->GetNumFaces(); i++)
    {
        this->mFaces[i]->FaceRemoveElement(this->mIndex);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::ResetIndex(const unsigned newIndex)
{
    for (unsigned i = 0; i < this->GetNumFaces(); ++i)
    {
        this->mFaces[i]->FaceRemoveElement(this->mIndex);
        this->mFaces[i]->FaceAddElement(newIndex);
    }

    this->MutableElement<ELEMENT_DIM, SPACE_DIM>::ResetIndex(newIndex);
}

//////////////////////////////////////////////////
///        Element manipulate faces            ///
//////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* VertexElement<ELEMENT_DIM, SPACE_DIM>::GetFace(const unsigned localIndex) const
{
    assert(localIndex < mFaces.size());
    return mFaces[localIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetFaceLocalIndex(const unsigned globalIndex) const
{
    unsigned local_index = UINT_MAX;
    for (unsigned i = 0; i < this->mFaces.size(); i++)
    {
        if (this->mFaces[i]->GetIndex() == globalIndex)
        {
            local_index = i;
            break;
        }
    }

    return local_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceIsOrientatedAntiClockwise(const unsigned localIndex) const
{
    assert(localIndex < mOrientations.size());
    return mOrientations[localIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::GetFaceOrientationWithGlobalIndex(const unsigned globalIndex) const
{
    const unsigned local_index = this->GetFaceLocalIndex(globalIndex);
    if (local_index == UINT_MAX)
    {
        NEVER_REACHED;
    }
    if (local_index >= this->mFaces.size())
    {
        NEVER_REACHED;
    }

    return mOrientations[local_index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddFace(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    ///\todo: duplicated code here because AddFace with the +1!!! #2850
    // // Add pFace to the end of mFaces
    // this->AddFace(pFace, false, this->GetNumFaces() - 1);
    // return;

    this->mFaces.push_back(pFace);
    this->mOrientations.push_back(false);

    pFace->FaceAddElement(this->mIndex);
    pFace->RegisterFaceWithNodes();

    // Create a set of indices of nodes currently owned by this element
    std::set<unsigned> node_indices;
    for (unsigned local_index = 0; local_index < this->GetNumNodes(); local_index++)
    {
        node_indices.insert(this->GetNodeGlobalIndex(local_index));
    }

    // Loop over nodes owned by pFace
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes(); local_index++)
    {
        // If this node is not already owned by this element...
        unsigned global_index = pFace->GetNodeGlobalIndex(local_index);
        if (node_indices.find(global_index) == node_indices.end())
        {
            // ... then add it to the element (and vice versa)
            this->AddNode(pFace->GetNode(local_index), this->GetNumNodes() - 1);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddFace(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                    bool Orientation, const unsigned& rIndex)
{
    assert(rIndex < this->mFaces.size());

    // Add pFace to rIndex+1 element of mFaces pushing the others up
    this->mFaces.insert(this->mFaces.begin() + rIndex + 1, pFace);
    this->mOrientations.insert(this->mOrientations.begin() + rIndex + 1, Orientation);

    assert(mFaces.size() == mOrientations.size());

    pFace->FaceAddElement(this->mIndex);
    pFace->RegisterFaceWithNodes();

    // Create a set of indices of nodes currently owned by this element
    std::set<unsigned> node_indices;
    for (unsigned local_index = 0; local_index < this->GetNumNodes(); local_index++)
    {
        node_indices.insert(this->GetNodeGlobalIndex(local_index));
    }

    // Loop over nodes owned by pFace
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes(); local_index++)
    {
        // If this node is not already owned by this element...
        unsigned global_index = pFace->GetNodeGlobalIndex(local_index);
        if (node_indices.find(global_index) == node_indices.end())
        {
            // ... then add it to the element (and vice versa)
            this->AddNode(pFace->GetNode(local_index), this->GetNumNodes() - 1);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteFace(const unsigned index)
{
    assert(index < this->mFaces.size());

    // Remove element from the registry of the face
    mFaces[index]->FaceRemoveElement(this->mIndex);

    // Remove the face and orientation at index (removes face from element)
    this->mFaces.erase(this->mFaces.begin() + index);
    this->mOrientations.erase(this->mOrientations.begin() + index);

    assert(mFaces.size() == mOrientations.size());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteFace(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    const unsigned face_local_index = this->GetFaceLocalIndex(pFace->GetIndex());
    if (face_local_index == UINT_MAX)
    {
        EXCEPTION("Face is not in element and cannot be deleted!");
    }
    this->DeleteFace(face_local_index);
    
    return face_local_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::ReplaceFace(const unsigned oldFaceLocalIndex,
                                                        VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pNewFace, const bool newFaceOrientation)
{
    NEVER_REACHED;
}

template <>
void VertexElement<3, 3>::ReplaceFace(const unsigned oldFaceLocalIndex,
                                      VertexElement<2, 3>* pNewFace, const bool newFaceOrientation)
{
    assert(oldFaceLocalIndex < this->mFaces.size());

    // Remove element from the registry of the face
    this->mFaces[oldFaceLocalIndex]->FaceRemoveElement(this->mIndex);

    // Update the face at this location
    this->mFaces[oldFaceLocalIndex] = pNewFace;
    this->mOrientations[oldFaceLocalIndex] = newFaceOrientation;

    // Add the element to the registry of the face
    pNewFace->FaceAddElement(this->mIndex);
}

//////////////////////////////////////////////////////////////////////////
///            Implementation of faces (similar to element)            ///
//////////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceResetIndex(unsigned index)
{
    for (unsigned i = 0; i < this->GetNumNodes(); ++i)
    {
        this->mNodes[i]->RemoveFace(this->mIndex);
    }
    this->mIndex = index;
    this->RegisterFaceWithNodes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceAddNode(Node<SPACE_DIM>* pNode, const unsigned Index)
{
    const unsigned real_index = (Index == UINT_MAX - 11) ? this->GetNumNodes() - 1 : Index;

    assert(real_index < this->mNodes.size());

    // Add pNode to real_index+1 face of mNodes pushing the others up
    this->mNodes.insert(this->mNodes.begin() + real_index + 1, pNode);

    // Add face to this node
    pNode->AddFace(this->mIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceDeleteNode(const unsigned Index)
{
    assert(Index < this->mNodes.size());

    // Remove face from the node at this location
    this->mNodes[Index]->RemoveFace(this->mIndex);

    // Remove the node at Index (removes node from face)
    this->mNodes.erase(this->mNodes.begin() + Index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceDeleteNode(const Node<SPACE_DIM>* pNode)
{
    const unsigned node_local_index = this->GetNodeLocalIndex(pNode->GetIndex());
    if (node_local_index == UINT_MAX)
    {
        EXCEPTION("Node is not in face and cannot be deleted!");
    }
    FaceDeleteNode(node_local_index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceUpdateNode(const unsigned Index, Node<SPACE_DIM>* pNode)
{
    assert(Index < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[Index]->RemoveFace(this->mIndex);

    // Update the node at this location
    this->mNodes[Index] = pNode;

    // Add face to this node
    pNode->AddFace(this->mIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceUpdateNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
{
    const unsigned old_local_index(this->GetNodeLocalIndex(pOldNode->GetIndex()));
    if (old_local_index == UINT_MAX)
    {
        EXCEPTION("Face (" << this->mIndex << ") does not have Node (" << pOldNode->GetIndex() << ")!");
    }
    this->FaceUpdateNode(old_local_index, pNewNode);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::RegisterFaceWithNodes()
{
    for (unsigned i = 0; i < this->mNodes.size(); ++i)
    {
        this->mNodes[i]->AddFace(this->mIndex);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::MarkFaceAsDeleted()
{
    // Mark face as deleted
    this->mIsDeleted = true;

    // Update nodes in the face so they know they are not contained by it
    for (unsigned i = 0; i < this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveFace(this->mIndex);
    }
    this->mNodes.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceRearrangeNodes(const c_vector<double, SPACE_DIM>& PointOfView)
{
    NEVER_REACHED;
}

template <>
bool VertexElement<2, 3>::FaceRearrangeNodes(const c_vector<double, 3>& PointOfView)
{
    const c_vector<double, 3> centroid = this->GetCentroid();

    c_vector<double, 3> normal = centroid - PointOfView;
    normal /= norm_2(normal);
    c_vector<double, 3> e1 = this->mNodes[0]->rGetLocation() - centroid;
    // Gramm-Schmidt Process
    e1 -= inner_prod(e1, normal) * normal;
    e1 /= norm_2(e1);
    c_vector<double, 3> e2 = VectorProduct(normal, e1);

    std::vector<std::pair<double, Node<3>*> > angles_and_nodes;
    for (unsigned i = 0; i < this->GetNumNodes(); ++i)
    {
        const c_vector<double, 3> vec_tmp = this->GetNode(i)->rGetLocation() - centroid;
        double tmp_angle = atan2(inner_prod(vec_tmp, e2), inner_prod(vec_tmp, e1));
        if (tmp_angle < 0)
        {
            tmp_angle += 2 * M_PI;
        }
        angles_and_nodes.push_back(std::make_pair(tmp_angle, this->GetNode(i)));
    }
    std::sort(angles_and_nodes.begin(), angles_and_nodes.end());

    std::vector<Node<3>*> nodes_tmp(this->GetNumNodes());
    for (unsigned i = 0; i < nodes_tmp.size(); ++i)
    {
        nodes_tmp[i] = angles_and_nodes[i].second;
    }
    
    const bool return_val = (nodes_tmp != this->mNodes);
    if (return_val)
    {
        std::swap(nodes_tmp, this->mNodes);
    }
    return return_val;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceRearrangeNodes()
{
    NEVER_REACHED;
}

template <>
bool VertexElement<2, 2>::FaceRearrangeNodes()
{
    const c_vector<double, 2> centroid = this->GetCentroid();

    std::vector<std::pair<double, Node<2>*> > angles_and_nodes;
    for (unsigned i = 0; i < this->GetNumNodes(); ++i)
    {
        const c_vector<double, 2>& loc_tmp = this->GetNode(i)->rGetLocation() - centroid;
        double tmp_angle = atan2(loc_tmp[1], loc_tmp[0]);
        if (tmp_angle < 0)
        {
            tmp_angle += 2 * M_PI;
        }
        angles_and_nodes.push_back(std::make_pair(tmp_angle, this->GetNode(i)));
    }
    std::sort(angles_and_nodes.begin(), angles_and_nodes.end());

    std::vector<Node<2>*> nodes_tmp(this->GetNumNodes());
    for (unsigned i = 0; i < nodes_tmp.size(); ++i)
    {
        nodes_tmp[i] = angles_and_nodes[i].second;
    }

    const bool return_val = (nodes_tmp != this->mNodes);
    if (return_val)
    {
        std::swap(nodes_tmp, this->mNodes);
    }
    return return_val;
}

//////////////////////////////////////////////////////////////////////////
///               Tracking element which contain this face             ///
//////////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceAddElement(unsigned elementIndex)
{
    mFaceContainingElementIndices.insert(elementIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceRemoveElement(unsigned elementIndex)
{
    unsigned count = mFaceContainingElementIndices.erase(elementIndex);
    if (count == 0)
    {
        EXCEPTION("Tried to remove an element index(" << elementIndex << ") from face(" << this->GetIndex() << ") which was not in the set");
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned>& VertexElement<ELEMENT_DIM, SPACE_DIM>::rFaceGetContainingElementIndices()
{
    return mFaceContainingElementIndices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceGetNumContainingElements() const
{
    return mFaceContainingElementIndices.size();
}

//////////////////////////////////////////////////////////////////////////
///                   Functions for monolayer elements                 ///
//////////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::CheckFaceOrientationOfElement(const unsigned faceLocalIndex)
{
    assert(ELEMENT_DIM == 3 && SPACE_DIM == 3);

    // Check Orientation
    const c_vector<double, SPACE_DIM> elem_centroid = this->GetCentroid();
    
    const VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = this->GetFace(faceLocalIndex);
    c_vector<double, SPACE_DIM> face_centroid = p_face->GetCentroid();
    c_vector<double, SPACE_DIM> face_normal = zero_vector<double>(SPACE_DIM);

    // Calculation w.r.t. centroid so that it will not be affected by uneven surface.
    c_vector<double, SPACE_DIM> v_minus_c0 = p_face->GetNode(p_face->GetNumNodes() - 1)->rGetLocation() - face_centroid;
    for (unsigned node_index = 0; node_index < p_face->GetNumNodes(); ++node_index)
    {
        c_vector<double, SPACE_DIM> vnext_minus_c0 = p_face->GetNode(node_index)->rGetLocation() - face_centroid;
        face_normal += VectorProduct(v_minus_c0, vnext_minus_c0);
        v_minus_c0 = vnext_minus_c0;
    }
    this->mOrientations[faceLocalIndex] = inner_prod(face_normal, face_centroid - elem_centroid) < 0;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::LateralFaceRearrangeNodes()
{
    NEVER_REACHED;
}

template <>
void VertexElement<2, 3>::LateralFaceRearrangeNodes()
{
    if (!IsLateralFace(this))
    {
        EXCEPTION("This function should be only called by lateral face.");
    }
    if (this->GetNumNodes() < 4)
    {
        return;
    }

    // Store the local indices of basal nodes for later use
    unsigned min_cyc_index = UINT_MAX;
    unsigned second_index = UINT_MAX;
    for (unsigned i = 0; i < this->GetNumNodes(); ++i)
    {
        if (IsBasalNode(this->GetNode(i)))
        {
            if (min_cyc_index == UINT_MAX)
            {
                min_cyc_index = i;
            }
            else
            {
                assert(second_index == UINT_MAX);
                second_index = i;
            }
        }
    }

    /*
     * if {0,num_nodes-1}, then the min should be actually num_nodes-1
     * other cases of {i, i+1}, just remain as they are.
     *
     * Then rotate mNodes accordingly.
     */
    if (min_cyc_index == 0u && second_index == this->GetNumNodes() - 1)
    {
        min_cyc_index = second_index;
    }
    else if (min_cyc_index + 1 != second_index)
    {
        NEVER_REACHED;
    }

    std::rotate(this->mNodes.begin(), this->mNodes.begin() + min_cyc_index, this->mNodes.end());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::MonolayerElementRearrangeFacesNodes()
{
    NEVER_REACHED;
}

template <>
void VertexElement<3, 3>::MonolayerElementRearrangeFacesNodes()
{
    if (!IsMonolayerElement(this))
    {
        NEVER_REACHED;
    }

    const VertexElement<2, 3>* p_basal = GetBasalFace(this);
    const VertexElement<2, 3>* p_apical = GetApicalFace(this);

    // Check and change basal and apical faces to 0th and 1st position respectively
    if (this->GetFace(0)->GetIndex() != p_basal->GetIndex())
    {
        const unsigned tmp_index = this->GetFaceLocalIndex(p_basal->GetIndex());
        assert(tmp_index != UINT_MAX);
        std::swap(this->mFaces[0], this->mFaces[tmp_index]);
        this->mOrientations.swap(this->mOrientations[0], this->mOrientations[tmp_index]);
    }
    if (this->GetFace(1)->GetIndex() != p_apical->GetIndex())
    {
        const unsigned tmp_index = this->GetFaceLocalIndex(p_apical->GetIndex());
        assert(tmp_index != UINT_MAX);
        std::swap(this->mFaces[1], this->mFaces[tmp_index]);
        this->mOrientations.swap(this->mOrientations[1], this->mOrientations[tmp_index]);
    }

    // Check and rearrange nodes if necessary
    std::vector<Node<3>*> elem_nodes(this->mNodes);
    std::vector<Node<3>*> face_nodes;
    face_nodes.reserve(this->GetNumNodes());
    {
        const std::vector<Node<3>*> lateral_nodes = GetNodesWithType(p_basal, Monolayer::BasalValue);
        face_nodes.insert(face_nodes.end(), lateral_nodes.begin(), lateral_nodes.end());
    }
    {
        const std::vector<Node<3>*> lateral_nodes = GetNodesWithType(p_apical, Monolayer::ApicalValue);
        face_nodes.insert(face_nodes.end(), lateral_nodes.begin(), lateral_nodes.end());
    }
    {
        const std::vector<Node<3>*> lateral_nodes = GetNodesWithType(this, Monolayer::LateralValue);
        face_nodes.insert(face_nodes.end(), lateral_nodes.begin(), lateral_nodes.end());
    }

    assert(elem_nodes.size() == face_nodes.size());
    // If nodes in element are not in proper order, they are rearranged
    if (elem_nodes != face_nodes)
    {
        const std::set<Node<3>*> tmp_elem_nodes(elem_nodes.begin(), elem_nodes.end());
        const std::set<Node<3>*> tmp_face_nodes(face_nodes.begin(), face_nodes.end());

        // Make sure the nodes are the same.
        if (tmp_elem_nodes != tmp_face_nodes)
        {
            NEVER_REACHED;
        }
        this->mNodes = face_nodes;
    }

    {
        std::set<Node<3>*> nodes_tmp;
        for (unsigned i = 0; i < this->GetNumFaces(); ++i)
        {
            VertexElement<2, 3>* p_face = this->GetFace(i);
            if (IsBasalFace(p_face) || IsApicalFace(p_face))
            {
                continue;
            }
            for (unsigned j = 0; j < p_face->GetNumNodes(); ++j)
            {
                nodes_tmp.insert(p_face->GetNode(j));
            }
        }

        const std::set<Node<3>*> actual_nodes(this->mNodes.begin(), this->mNodes.end());
        if (nodes_tmp != actual_nodes)
        {
            NEVER_REACHED;
        }
    }

    // With new functionalities, it is not necessary to arrange faces in "correct" order, as there 
    // isn't unambiguous way to arrange faces anymore.

    // Check Orientation
    for (unsigned face_index = 0; face_index < this->GetNumFaces(); ++face_index)
    {
        this->CheckFaceOrientationOfElement(face_index);
    }

    ///\todo: #2850 remove the both assumptions of this function (refer to hpp)
}

////////////////////////////////////////////////////////////////////////
///                  Specialization for 1d elements                  ///
///                                                                  ///
///                 1d elements are just edges (lines)               ///
////////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template <unsigned SPACE_DIM>
VertexElement<1, SPACE_DIM>::VertexElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
        : MutableElement<1, SPACE_DIM>(index, rNodes)
{
}

template <unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::GetNumFaces() const
{
    return 0;
}

template <unsigned SPACE_DIM>
VertexElement<0, SPACE_DIM>* VertexElement<1, SPACE_DIM>::GetFace(unsigned index) const
{
    return NULL;
}

template <unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::GetFaceLocalIndex(unsigned globalIndex) const
{
    NEVER_REACHED;
}

template <unsigned SPACE_DIM>
bool VertexElement<1, SPACE_DIM>::FaceIsOrientatedAntiClockwise(const unsigned index) const
{
    return false;
}

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<1, SPACE_DIM>::GetCentroid() const
{
    if (this->GetNumNodes() != 2u)
    {
        EXCEPTION("1D line should only has 2 nodes!");
    }
    return 0.5 * (this->GetNodeLocation(0) + this->GetNodeLocation(1));
}

template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::AddFace(VertexElement<0, SPACE_DIM>* pFace)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::AddFace(VertexElement<0, SPACE_DIM>* pFace, bool Orientation, const unsigned& rIndex)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::DeleteFace(const unsigned index)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::DeleteFace(VertexElement<0, SPACE_DIM>* pFace)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::ReplaceFace(const unsigned oldFaceLocalIndex,
                                              VertexElement<0, SPACE_DIM>* pNewFace, const bool newFaceOrientation)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::MonolayerElementRearrangeFacesNodes()
{
    NEVER_REACHED;
}

template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::FaceResetIndex(unsigned index)
{
    for (unsigned i = 0; i < this->GetNumNodes(); ++i)
    {
        this->mNodes[i]->RemoveFace(this->mIndex);
    }
    this->mIndex = index;
    this->RegisterFaceWithNodes();
}

template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::FaceAddNode(Node<SPACE_DIM>* pNode, const unsigned Index)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::FaceDeleteNode(const unsigned Index)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::FaceDeleteNode(const Node<SPACE_DIM>* pNode)
{
    NEVER_REACHED;
}
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::FaceUpdateNode(const unsigned Index, Node<SPACE_DIM>* pNode)
{
    NEVER_REACHED;
}

template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::RegisterFaceWithNodes()
{
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddFace(this->mIndex);
    }
}

template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::MarkFaceAsDeleted()
{
    // Mark face as deleted
    this->mIsDeleted = true;

    // Update nodes in the face so they know they are not contained by it
    for (unsigned i = 0; i < this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveFace(this->mIndex);
    }
}

//////////////////////////////////////////////////////////////////////////
///               Tracking element which contain this face             ///
//////////////////////////////////////////////////////////////////////////
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::FaceAddElement(const unsigned elementIndex)
{
    mFaceContainingElementIndices.insert(elementIndex);
}

template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::FaceRemoveElement(const unsigned elementIndex)
{
    unsigned count = mFaceContainingElementIndices.erase(elementIndex);
    if (count == 0)
    {
        EXCEPTION("Tried to remove an element index(" << elementIndex << ") from face(" << this->GetIndex() << ") which was not in the set");
    }
}

template <unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::FaceGetNumContainingElements() const
{
    return mFaceContainingElementIndices.size();
}

template <unsigned SPACE_DIM>
std::set<unsigned>& VertexElement<1, SPACE_DIM>::rFaceGetContainingElementIndices()
{
    return mFaceContainingElementIndices;
}

// Explicit instantiation
template class VertexElement<1, 1>;
template class VertexElement<1, 2>;
template class VertexElement<1, 3>;
template class VertexElement<2, 2>;
template class VertexElement<2, 3>;
template class VertexElement<3, 3>;
