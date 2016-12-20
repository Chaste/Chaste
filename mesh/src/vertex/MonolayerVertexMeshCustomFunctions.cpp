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

#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "Node.hpp"
#include "VertexElement.hpp"
#include <set>
#include <algorithm>


/// ===============================================================
/// Some function that can be added into trunk and relevant for all

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ElementHasNode(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, const unsigned nodeIndex)
{
    return pElement->GetNodeLocalIndex(nodeIndex)!=UINT_MAX;
}

std::set<unsigned> GetSharedElementIndices(const Node<3>* pNodeA, const Node<3>* pNodeB)
{
    Node<3>* p_node_a = const_cast<Node<3>*>(pNodeA);
    Node<3>* p_node_b = const_cast<Node<3>*>(pNodeB);
    std::set<unsigned> elems_with_node_a = p_node_a->rGetContainingElementIndices();
    std::set<unsigned> elems_with_node_b = p_node_b->rGetContainingElementIndices();

    std::set<unsigned> shared_elements;
    std::set_intersection(elems_with_node_a.begin(), elems_with_node_a.end(),
                          elems_with_node_b.begin(),elems_with_node_b.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    return shared_elements;
}

std::set<unsigned> GetSharedFaceIndices(const Node<3>* pNodeA, const Node<3>* pNodeB)
{
    Node<3>* p_node_a = const_cast<Node<3>*>(pNodeA);
    Node<3>* p_node_b = const_cast<Node<3>*>(pNodeB);
    std::set<unsigned> elems_with_node_a = p_node_a->rGetContainingFaceIndices();
    std::set<unsigned> elems_with_node_b = p_node_b->rGetContainingFaceIndices();

    std::set<unsigned> shared_faces;
    std::set_intersection(elems_with_node_a.begin(), elems_with_node_a.end(),
                          elems_with_node_b.begin(),elems_with_node_b.end(),
                          std::inserter(shared_faces, shared_faces.begin()));

    return shared_faces;
}


/**
 * Face is only considered a boundary face when all nodes are boundary nodes.
 * @return
 */
bool IsFaceOnBoundary(VertexElement<2, 3>* pFace)
{
    ///\todo: #2850 cannot handle such case
    /*
     * The new node is boundary node if the 2 nodes are boundary nodes and the elements don't look like
     *   ___A___
     *  |   |   |
     *  |___|___|
     *      B
     */
    bool is_element_on_boundary = true;
    for (unsigned i=0; i<pFace->GetNumNodes(); ++i)
    {
        if (!pFace->GetNode(i)->IsBoundaryNode())
        {
            is_element_on_boundary = false;
            break;
        }
    }

    return is_element_on_boundary;
}


/// ===============================================================
/// Functions for monolayer




void SetNodeAsApical(Node<3>* pNode)
{
    assert(pNode->GetNumNodeAttributes() == 0u);    // LCOV_EXCL_LINE
    pNode->AddNodeAttribute(Monolayer::SetApicalValue);
}

void SetNodeAsBasal(Node<3>* pNode)
{
    assert(pNode->GetNumNodeAttributes() == 0u);    // LCOV_EXCL_LINE
    pNode->AddNodeAttribute(Monolayer::SetBasalValue);
}

Monolayer::v_type GetNodeType(const Node<3>* pNode)
{
    /* implemented const node as there is no modication for this
     * function. However, there isn't suitable const function for
     * NodeAttribute and hence required a const_cast
     */
    Node<3>* p_non_const_node = const_cast<Node<3>*>(pNode);
    assert(p_non_const_node->GetNumNodeAttributes() == 1u);    // LCOV_EXCL_LINE
    return Monolayer::v_type(p_non_const_node->rGetNodeAttributes()[0]);
}

bool IsApicalNode(const Node<3>* pNode)
{
    return GetNodeType(pNode)==Monolayer::ApicalValue;
}

bool IsBasalNode(const Node<3>* pNode)
{
    return GetNodeType(pNode)==Monolayer::BasalValue;
}



// VertexElement
void SetFaceAsApical(VertexElement<2, 3>* pFace)
{
    assert(pFace->GetNumElementAttributes() == 0u);    // LCOV_EXCL_LINE
    pFace->AddElementAttribute(Monolayer::SetApicalValue);
}

void SetFaceAsBasal(VertexElement<2, 3>* pFace)
{
    assert(pFace->GetNumElementAttributes() == 0u);    // LCOV_EXCL_LINE
    pFace->AddElementAttribute(Monolayer::SetBasalValue);
}

void SetFaceAsLateral(VertexElement<2, 3>* pFace)
{
    assert(pFace->GetNumElementAttributes() == 0u);    // LCOV_EXCL_LINE
    pFace->AddElementAttribute(Monolayer::SetLateralValue);
}

Monolayer::v_type GetFaceType(const VertexElement<2, 3>* pFace)
{
    /* implemented const node as there is no modication for this
     * function. However, there isn't suitable const function for
     * ElementAttribute and hence required a const_cast
     */
    VertexElement<2, 3>* p_non_const_face = const_cast<VertexElement<2, 3>*>(pFace);
    assert(p_non_const_face->GetNumElementAttributes() == 1u);    // LCOV_EXCL_LINE
    return Monolayer::v_type(p_non_const_face->rGetElementAttributes()[0]);
}



bool IsApicalFace(const VertexElement<2, 3>* pFace)
{
    return GetFaceType(pFace)==Monolayer::ApicalValue;
}

bool IsBasalFace(const VertexElement<2, 3>* pFace)
{
    return GetFaceType(pFace)==Monolayer::BasalValue;
}

bool IsLateralFace(const VertexElement<2, 3>* pFace)
{
    return GetFaceType(pFace)==Monolayer::LateralValue;
}

void SetElementAsMonolayer(VertexElement<3, 3>* pElement)
{
    assert(pElement->GetNumElementAttributes() == 0u);    // LCOV_EXCL_LINE
    pElement->AddElementAttribute(Monolayer::SetElementValue);
    pElement->MonolayerElementRearrangeFacesNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool IsMonolayerElement(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    if (ELEMENT_DIM!=3 || SPACE_DIM!=3)
    {
        return false;
    }
    else
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_non_const_elem = const_cast<VertexElement<ELEMENT_DIM, SPACE_DIM>*>(pElement);
        return (p_non_const_elem->GetNumElementAttributes() == 1u) &&
               (Monolayer::v_type(p_non_const_elem->rGetElementAttributes()[0]) == Monolayer::ElementValue);
    }
}
// Template instantiation is for the compiler
template bool IsMonolayerElement<2, 2>(const VertexElement<2, 2>* pElement);
template bool IsMonolayerElement<2, 3>(const VertexElement<2, 3>* pElement);
template bool IsMonolayerElement<3, 3>(const VertexElement<3, 3>* pElement);


VertexElement<2, 3>* GetApicalFace(const VertexElement<3, 3>* pElement)
{
    VertexElement<2, 3>* p_face = pElement->GetFace(1);
    assert(IsApicalFace(p_face));   // LCOV_EXCL_LINE
    return p_face;
}

VertexElement<2, 3>* GetBasalFace(const VertexElement<3, 3>* pElement)
{
    VertexElement<2, 3>* p_face = pElement->GetFace(0);
    assert(IsBasalFace(p_face));   // LCOV_EXCL_LINE
    return p_face;
}

/**
 * Get half the number of nodes of a monolayer vertex element
 * to eliminate /2 everywhere in the code.
 * @param pElement
 */
unsigned MonolayerGetNumNodes(const VertexElement<3, 3>* pElement)
{
    unsigned num_nodes = pElement->GetNumNodes();
    assert(num_nodes%2 == 0);   // LCOV_EXCL_LINE
    num_nodes /= 2;
    assert(num_nodes == GetApicalFace(pElement)->GetNumNodes());  // LCOV_EXCL_LINE
    assert(num_nodes == GetBasalFace(pElement)->GetNumNodes());  // LCOV_EXCL_LINE
    return num_nodes;
}

std::vector<unsigned> GetLateralFace(const VertexElement<3, 3>* pElement, const unsigned nodeIndexA, const unsigned nodeIndexB)
{
    std::vector<unsigned> return_vector;
    for (unsigned face_local_index=2; face_local_index<pElement->GetNumFaces(); ++face_local_index)
    {
        VertexElement<2, 3>* p_tmp_face = pElement->GetFace(face_local_index);

        if (ElementHasNode(p_tmp_face, nodeIndexA) && ElementHasNode(p_tmp_face, nodeIndexB))
        {
            return_vector.push_back(p_tmp_face->GetIndex());
            return_vector.push_back(pElement->FaceIsOrientatedAntiClockwise(face_local_index));
            return_vector.push_back(face_local_index);
            break;
        }
    }

    if (return_vector.size() == 0)
    {
        return_vector.push_back(UINT_MAX);
    }
    else
    {
        assert(return_vector.size() == 3);  // LCOV_EXCL_LINE
    }
    return return_vector;
}

bool GetFaceOrientation(const VertexElement<3, 3>* pElement, const unsigned faceIndex)
{
    bool face_orientation;
    bool face_found (false);
    for (unsigned face_index = 0; face_index<pElement->GetNumFaces(); ++face_index)
    {
        if (faceIndex == pElement->GetFace(face_index)->GetIndex())
        {
            face_orientation = pElement->FaceIsOrientatedAntiClockwise(face_index);
            face_found = true;
            break;
        }
    }
    assert(face_found);
    return face_orientation;
}




//void AddPairNode(VertexElement<3, 3>* pElement, const unsigned index, Node<3>* pBasalNode, Node<3>* pApicalNode)
//{
//    const unsigned num_nodes = MonolayerGetNumNodes(pElement);
//    unsigned smaller_index = index;
//    if (index>num_nodes)
//    {
//        smaller_index -= num_nodes;
//        assert(smaller_index < num_nodes);
//    }
//
//    GetBasalFace(pElement)->FaceAddNode(pBasalNode, smaller_index);
//    GetApicalFace(pElement)->FaceAddNode(pApicalNode, smaller_index);
//
//    // Add the apical node first so that I don't changes its indices.
//    pElement->AddNode(pApicalNode, smaller_index+num_nodes);
//    pElement->AddNode(pBasalNode, smaller_index);
//
//    // NOTE: I doesn't rearrange here because there is still a missing face.
//}
