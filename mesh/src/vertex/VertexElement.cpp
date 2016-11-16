/*

Copyright (c) 2005-2016, University of Oxford.
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
#include <cassert>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
      mFaces(rFaces),
      mOrientations(rOrientations)
{
    // This constructor should only be used in 3D
    assert(SPACE_DIM == 3);    // LCOV_EXCL_LINE - code will be removed at compile time

    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    if (SPACE_DIM == ELEMENT_DIM)
    {
        // Register element with nodes
        this->RegisterWithNodes();
    }
}

template<unsigned SPACE_DIM>
class NodeLessThanYxz
{
public:
    static const double DELTA_YXZ=1e-6;
    bool operator() (const Node<SPACE_DIM>* pNode1, const Node<SPACE_DIM>* pNode2) const
    {
        const c_vector<double, SPACE_DIM>& location1 = pNode1->rGetLocation();
        const c_vector<double, SPACE_DIM>& location2 = pNode2->rGetLocation();
        bool result = false;

        // As double doesn't save the exact value, direct comparison could lead to error
        // e.g. comparison of 0.50000000000000011 and 0.50000000000000044
        std::vector<double> is_same(SPACE_DIM, false);
        for (unsigned i=0; i < SPACE_DIM; i++)
        {
            if (fabs(location1[i] - location2[i]) < DELTA_YXZ)
            {
                is_same[i] = true;
            }
        }

        // Only if the two values aren't the same, the < comparison will make sense
        if (!is_same[1] && (location1[1] < location2[1]))
        {
            result = true;
        }
        else if (is_same[1])
        {
            if (!is_same[0] && location1[0] < location2[0])
            {
                result = true;
            }
            else if ((SPACE_DIM == 3) && is_same[0] && !is_same[2] && (location1[2] < location2[2]))
            {
                result = true;
            }
        }
        return result;
    }
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index),
      mFaces(rFaces),
      mOrientations(rOrientations)
{
//EXCEPTION("THIS CONSTRUCTOR SHOULD NOT BE USED YET"); ///\todo #2850
    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    // Make a set of nodes with mFaces
    std::set<Node<SPACE_DIM>*, NodeLessThanYxz<SPACE_DIM> > nodes_set;
    for (unsigned face_index=0; face_index<mFaces.size(); face_index++)
    {
        for (unsigned node_index=0; node_index<mFaces[face_index]->GetNumNodes(); node_index++)
        {
            nodes_set.insert(mFaces[face_index]->GetNode(node_index));
        }
    }

    // Populate mNodes
    for (typename std::set< Node<SPACE_DIM>*, NodeLessThanYxz<SPACE_DIM> >::iterator node_iter = nodes_set.begin();
         node_iter != nodes_set.end();
         ++node_iter)
    {
         this->mNodes.push_back(*node_iter);
    }

    // Register element with nodes
    this->RegisterWithNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::~VertexElement()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
{
    // first call the function implementated in abstract Element to rnt
    AbstractElement<ELEMENT_DIM, SPACE_DIM>::ReplaceNode(pOldNode, pNewNode);

    for (unsigned face_index=0 ; face_index<mFaces.size() ; ++face_index)
    {
        VertexElement<ELEMENT_DIM-1, SPACE_DIM>& r_face = *(mFaces[face_index]);
        const unsigned node_local_index = r_face.GetNodeLocalIndex(pOldNode->GetIndex());
        if ( node_local_index != UINT_MAX )
            r_face.FaceUpdateNode(node_local_index, pNewNode);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddFace(VertexElement<ELEMENT_DIM-1,SPACE_DIM>* pFace)
{
    // Add pFace to the end of mFaces
    this->mFaces.push_back(pFace);

    // Create a set of indices of nodes currently owned by this element
    std::set<unsigned> node_indices;
    for (unsigned local_index=0; local_index<this->GetNumNodes(); local_index++)
    {
        node_indices.insert(this->GetNodeGlobalIndex(local_index));
    }

    // Loop over nodes owned by pFace
    unsigned end_index = this->GetNumNodes()-1;
    for (unsigned local_index=0; local_index<pFace->GetNumNodes(); local_index++)
    {
        // If this node is not already owned by this element...
        unsigned global_index = pFace->GetNodeGlobalIndex(local_index);
        if (node_indices.find(global_index) == node_indices.end())
        {
            // ... then add it to the element (and vice versa)
            this->AddNode(pFace->GetNode(local_index), end_index);
            end_index++;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM-1,  SPACE_DIM>* VertexElement<ELEMENT_DIM, SPACE_DIM>::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
    return mFaces[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceIsOrientatedAntiClockwise(unsigned index) const
{
    assert(index < mOrientations.size());
    return mOrientations[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::MonolayerElementDeleteNodes(const Node<SPACE_DIM>* pApicalNode1,
        const Node<SPACE_DIM>* pBasalNode1, Node<SPACE_DIM>* pApicalNode2, Node<SPACE_DIM>* pBasalNode2)
{
    NEVER_REACHED;
}

template<>
unsigned VertexElement<3, 3>::MonolayerElementDeleteNodes(const Node<3>* pApicalNodeDelete,
        const Node<3>* pBasalNodeDelete, Node<3>* pApicalNodeStay, Node<3>* pBasalNodeStay)
{
    const unsigned apical_node_delete_index = pApicalNodeDelete->GetIndex();
    const unsigned basal_node_delete_index = pBasalNodeDelete->GetIndex();
    const unsigned apical_node_stay_index = pApicalNodeStay->GetIndex();

    // We will delete the node1's from the element and faces and also remove the face from the element.
    // Operation for faces: I apical & basal faces will just remove the node1's
    // II the lateral face not belong to the deleted element needs to update node1's to node2's
    // III the lateral face belong to the deleted element need to be remove from the registry of neighbouring element
    this->DeleteNode(this->GetNodeLocalIndex(apical_node_delete_index));
    this->DeleteNode(this->GetNodeLocalIndex(basal_node_delete_index));
    // Operation I
    this->GetFace(0)->FaceDeleteNode(pBasalNodeDelete);
    this->GetFace(1)->FaceDeleteNode(pApicalNodeDelete);

    unsigned delete_lateral_face_index(UINT_MAX);
    // Operation II and III To update another face and remove the extra face
    for (unsigned face_index=2 ; face_index<this->GetNumFaces() ; ++face_index)
    {
        VertexElement<2, 3>* p_face = this->GetFace(face_index);
        const unsigned apical_node1_local_index = p_face->GetNodeLocalIndex(apical_node_delete_index);
        if (apical_node1_local_index != UINT_MAX)
        {
            if (p_face->GetNodeLocalIndex(apical_node_stay_index) == UINT_MAX)
            {
                // Operation II
                p_face->FaceUpdateNode(apical_node1_local_index, pApicalNodeStay);
                const unsigned basal_node1_local_index = p_face->GetNodeLocalIndex(basal_node_delete_index);
                p_face->FaceUpdateNode(basal_node1_local_index, pBasalNodeStay);
            }
            else
            {
                // Operation III
                delete_lateral_face_index = p_face->GetIndex();
                this->DeleteFace(face_index);

                // To compensate the changes in NumFaces
                --face_index;
            }
        }
    }
    return delete_lateral_face_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceUpdateNode(const unsigned Index, Node<SPACE_DIM>* pNode)
{
    assert(Index < this->mNodes.size());

    // Update the node at this location
    this->mNodes[Index] = pNode;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceDeleteNode(const unsigned Index)
{
    assert(Index < this->mNodes.size());

    // Remove the node at rIndex (removes node from element)
    this->mNodes.erase(this->mNodes.begin() + Index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceDeleteNode(const Node<SPACE_DIM>* pNode)
{
    const unsigned node_local_index = this->GetNodeLocalIndex(pNode->GetIndex());
    if (node_local_index == UINT_MAX)
    {
        EXCEPTION("Node is not in face and cannot be deleted!");
    }
    FaceDeleteNode(node_local_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceAddNode(Node<SPACE_DIM>* pNode, const unsigned Index)
{
    assert(Index < this->mNodes.size());

    // Add pNode to rIndex+1 element of mNodes pushing the others up
    this->mNodes.insert(this->mNodes.begin() + Index+1,  pNode);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteFace(const unsigned& rIndex)
{
    assert(rIndex < this->mFaces.size());

    // Remove the face and orientation at rIndex (removes face from element)
    this->mFaces.erase(this->mFaces.begin() + rIndex);
    this->mOrientations.erase(this->mOrientations.begin() + rIndex);

    assert(mFaces.size() == mOrientations.size());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::MarkFaceAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace,
                                                    bool& Orientation,const unsigned& rIndex)
{
    assert(rIndex < this->mFaces.size());

    // Add pFace to rIndex+1 element of mFaces pushing the others up
    this->mFaces.insert(this->mFaces.begin() + rIndex+1,  pFace);
    this->mOrientations.insert(this->mOrientations.begin() +rIndex+1, Orientation);

    assert(mFaces.size() == mOrientations.size());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
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
            for (unsigned local_index=0; local_index<num_nodes; local_index++)
            {
                c_vector<double, SPACE_DIM> next_node_location = this->GetNodeLocation((local_index+1)%num_nodes);
                c_vector<double, SPACE_DIM> pos_2 = next_node_location - first_node_location;

                double this_x = pos_1[0];
                double this_y = pos_1[1];
                double next_x = pos_2[0];
                double next_y = pos_2[1];

                double signed_area_term = this_x*next_y - this_y*next_x;

                centroid_x += (this_x + next_x)*signed_area_term;
                centroid_y += (this_y + next_y)*signed_area_term;
                element_signed_area += 0.5*signed_area_term;

                pos_1 = pos_2;
            }

            assert(element_signed_area != 0.0);

            // Finally, map back and employ GetVectorFromAtoB() to allow for periodicity
            centroid = first_node_location;
            centroid(0) += centroid_x / (6.0*element_signed_area);
            centroid(1) += centroid_y / (6.0*element_signed_area);
        }
        break;
        case 3:
        {
            ///\todo compute centroid rather than centre of mass (see #1422)
            for (unsigned local_index=0; local_index<num_nodes; ++local_index)
            {
                centroid += this->GetNodeLocation(local_index);//PRINT_VARIABLE(p_element->GetNodeLocation(local_index));
            }
            centroid /= ((double) num_nodes);
        }
        break;
        default:
            NEVER_REACHED;
    }
    return centroid;
}

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
VertexElement<1, SPACE_DIM>::VertexElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<1, SPACE_DIM>(index, rNodes)
{
}

template<unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::GetNumFaces() const
{
    return 0;
}

template<unsigned SPACE_DIM>
VertexElement<0, SPACE_DIM>* VertexElement<1, SPACE_DIM>::GetFace(unsigned index) const
{
    return NULL;
}

template<unsigned SPACE_DIM>
bool VertexElement<1, SPACE_DIM>::FaceIsOrientatedAntiClockwise(unsigned index) const
{
    return false;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<1, SPACE_DIM>::GetCentroid() const
{
    return 0.5*(this->GetNodeLocation(0) + this->GetNodeLocation(1));
}

// Explicit instantiation
template class VertexElement<1,1>;
template class VertexElement<1,2>;
template class VertexElement<1,3>;
template class VertexElement<2,2>;
template class VertexElement<2,3>;
template class VertexElement<3,3>;
