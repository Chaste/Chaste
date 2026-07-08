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

#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "MutableVertexMesh.hpp"
#include "Node.hpp"
#include "VertexElement.hpp"

#include <algorithm>
#include <map>
#include <string>

#include <iomanip> // PrintElement
#include <iostream> // PrintMesh and PrintElement

#include "Debug.hpp"

/////////////////////////////////////////////////////////
///      Some function that are relevant for all      ///
/////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ElementHasNode(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, const unsigned nodeIndex)
{
    return pElement->GetNodeLocalIndex(nodeIndex) != UINT_MAX;
}
// Template instantiation is for the compiler
template bool ElementHasNode(const VertexElement<1, 1>* pElement, const unsigned nodeIndex);
template bool ElementHasNode(const VertexElement<1, 2>* pElement, const unsigned nodeIndex);
template bool ElementHasNode(const VertexElement<1, 3>* pElement, const unsigned nodeIndex);
template bool ElementHasNode(const VertexElement<2, 2>* pElement, const unsigned nodeIndex);
template bool ElementHasNode(const VertexElement<2, 3>* pElement, const unsigned nodeIndex);
template bool ElementHasNode(const VertexElement<3, 3>* pElement, const unsigned nodeIndex);

std::set<unsigned> GetSharedElementIndices(const Node<3>* pNodeA, const Node<3>* pNodeB)
{
    std::set<unsigned> shared_elements;

    if (pNodeA != pNodeB)
    {
        Node<3>* p_node_a = const_cast<Node<3>*>(pNodeA);
        Node<3>* p_node_b = const_cast<Node<3>*>(pNodeB);
        const std::set<unsigned>& elems_with_node_a = p_node_a->rGetContainingElementIndices();
        const std::set<unsigned>& elems_with_node_b = p_node_b->rGetContainingElementIndices();

        std::set_intersection(elems_with_node_a.begin(), elems_with_node_a.end(),
                              elems_with_node_b.begin(), elems_with_node_b.end(),
                              std::inserter(shared_elements, shared_elements.begin()));
    }

    return shared_elements;
}

std::set<unsigned> GetSharedFaceIndices(const Node<3>* pNodeA, const Node<3>* pNodeB)
{
    std::set<unsigned> shared_faces;

    if (pNodeA != pNodeB)
    {
        Node<3>* p_node_a = const_cast<Node<3>*>(pNodeA);
        Node<3>* p_node_b = const_cast<Node<3>*>(pNodeB);
        const std::set<unsigned>& elems_with_node_a = p_node_a->rGetContainingFaceIndices();
        const std::set<unsigned>& elems_with_node_b = p_node_b->rGetContainingFaceIndices();

        std::set_intersection(elems_with_node_a.begin(), elems_with_node_a.end(),
                              elems_with_node_b.begin(), elems_with_node_b.end(),
                              std::inserter(shared_faces, shared_faces.begin()));
    }

    return shared_faces;
}

VertexElement<2, 3>* GetSharedLateralFace(const VertexElement<3, 3>* pElemA,
                                          const VertexElement<3, 3>* pElemB)
{
    if (pElemA == nullptr || pElemB == nullptr)
    {
        EXCEPTION("Two elements do not share lateral face."); //LCOV_EXCL_LINE
    }

    std::set<unsigned> s1, s2, s_return;
    for (unsigned i = 0; i < pElemA->GetNumFaces(); ++i)
    {
        if (IsLateralFace(pElemA->GetFace(i)))
        {
            s1.insert(pElemA->GetFace(i)->GetIndex());
        }
    }
    for (unsigned i = 0; i < pElemB->GetNumFaces(); ++i)
    {
        if (IsLateralFace(pElemB->GetFace(i)))
        {
            s2.insert(pElemB->GetFace(i)->GetIndex());
        }
    }
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                          std::inserter(s_return, s_return.begin()));

    switch (s_return.size())
    {
        case 0:
            EXCEPTION("Two elements do not share lateral face."); //LCOV_EXCL_LINE
        case 1:
            return pElemA->GetFace(pElemA->GetFaceLocalIndex(*(s_return.begin())));
        default:
            EXCEPTION("Probably some errors occur. Two elements share more than 1 lateral face."); //LCOV_EXCL_LINE
    }
}

template <typename VertexObject>
VertexElement<2, 3>* GetSharedLateralFace(const Node<3>* pNodeA, const Node<3>* pNodeB, const VertexObject* pObject)
{
    const std::set<unsigned> shared_face_ids = GetSharedFaceIndices(pNodeA, pNodeB);
    const std::set<VertexElement<2, 3>*> shared_lateral_faces = GetFacesWithIndices(shared_face_ids, pObject,
                                                                                    Monolayer::LateralValue);

    switch (shared_lateral_faces.size())
    {
        case 0:
            EXCEPTION("Two nodes do not share lateral face."); //LCOV_EXCL_LINE
        case 1:
            return no1(shared_lateral_faces);
        default:
            EXCEPTION("Probably some errors occur. Two nodes share more than 1 lateral face in monolayer mesh."); //LCOV_EXCL_LINE
    }
}
template VertexElement<2, 3>* GetSharedLateralFace(const Node<3>*, const Node<3>*, const VertexMesh<3, 3>*);
template VertexElement<2, 3>* GetSharedLateralFace(const Node<3>*, const Node<3>*, const MutableVertexMesh<3, 3>*);
template VertexElement<2, 3>* GetSharedLateralFace(const Node<3>*, const Node<3>*, const VertexElement<3, 3>*);

std::set<VertexElement<2, 3>*> GetFacesWithIndices(const std::set<unsigned>& face_indices, const VertexElement<3, 3>* pElement,
                                                   const Monolayer::v_type faceType)
{
    std::set<VertexElement<2, 3>*> s_result;
    for (unsigned i = 0; i < pElement->GetNumFaces(); ++i)
    {
        VertexElement<2, 3>* p_tmp_face = pElement->GetFace(i);
        if (face_indices.count(p_tmp_face->GetIndex()) != 0)
        {
            s_result.insert(p_tmp_face);
        }
    }

    if (faceType != Monolayer::AllTypes)
    {
        std::set<VertexElement<2, 3>*> s_matched;
        for (VertexElement<2, 3>* p_result_face : s_result)
        {
            if (GetFaceType(p_result_face) == faceType)
            {
                s_matched.insert(p_result_face);
            }
        }

        s_matched.swap(s_result);
    }

    return s_result;
}

std::set<VertexElement<2, 3>*> GetFacesWithIndices(const std::set<unsigned>& face_indices, const VertexMesh<3, 3>* pMesh,
                                                   const Monolayer::v_type faceType)
{
    std::set<VertexElement<2, 3>*> s_result;
    for (unsigned face_index : face_indices)
    {
        s_result.insert(pMesh->GetFace(face_index));
    }

    if (faceType != Monolayer::AllTypes)
    {
        std::set<VertexElement<2, 3>*> s_matched;
        for (VertexElement<2, 3>* p_result_face : s_result)
        {
            if (GetFaceType(p_result_face) == faceType)
            {
                s_matched.insert(p_result_face);
            }
        }

        s_matched.swap(s_result);
    }

    return s_result;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PrintElement(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    NEVER_REACHED;
}

// LCOV_EXCL_START
template <>
void PrintElement(const VertexElement<3, 3>* pElement)
{
    const std::string TAB = "    ";

    std::cout << "=================================================================================" << std::endl;
    // Printing out each elements
    std::cout << "ELEMENT: " << pElement->GetIndex() << (pElement->IsDeleted() ? " (DELETED)" : "") << std::endl;

    std::cout << TAB << "number of Faces : " << pElement->GetNumFaces() << " {";
    for (unsigned j = 0; j < pElement->GetNumFaces(); ++j)
    {
        std::cout << std::setw(3) << pElement->GetFace(j)->GetIndex() << "  ";
    }
    std::cout << "}" << std::endl;

    std::cout << TAB << "Face oriented.. : " << pElement->GetNumFaces() << " {";
    for (unsigned j = 0; j < pElement->GetNumFaces(); ++j)
    {
        std::cout << std::setw(3) << pElement->FaceIsOrientatedAntiClockwise(j) << "  ";
    }
    std::cout << "}" << std::endl;

    std::cout << TAB << "number of Nodes : " << pElement->GetNumNodes() << " {  ";
    for (unsigned j = 0; j < pElement->GetNumNodes(); ++j)
    {
        std::cout << pElement->GetNode(j)->GetIndex() << "  ";
    }
    std::cout << "}" << std::endl;

    VertexElement<2, 3>& basal = *(pElement->GetFace(0));
    std::cout << TAB << "Nodes for basal face " << basal.GetIndex() << " {  ";
    for (unsigned j = 0; j < basal.GetNumNodes(); ++j)
    {
        std::cout << basal.GetNode(j)->GetIndex() << "  ";
    }
    std::cout << "}" << std::endl
              << "---------------------------------------------------------" << std::endl;

    std::cout << "***************************************************************" << std::endl;
    // Now printing all faces
    for (unsigned i = 0; i < pElement->GetNumFaces(); ++i)
    {
        VertexElement<2, 3>* p_face = pElement->GetFace(i);
        std::cout << "FACE (" << i << ") : " << p_face->GetIndex() << (p_face->IsDeleted() ? " (DELETED)" : "") << std::endl;
        std::cout << TAB << "Face Attribute : " << Monolayer::ValueToString[GetFaceType(p_face)] << (IsFaceOnBoundary(p_face) ? " (BOUNDARY)" : "") << std::endl;

        std::set<unsigned> set_tmp = p_face->rFaceGetContainingElementIndices();
        std::cout << TAB << "number of Elements : " << set_tmp.size() << " {  ";
        for (unsigned elem_index : set_tmp)
        {
            std::cout << elem_index << "  ";
        }
        std::cout << "}" << std::endl;

        std::cout << TAB << "number of Nodes : " << p_face->GetNumNodes() << " {  ";
        for (unsigned j = 0; j < p_face->GetNumNodes(); ++j)
        {
            std::cout << p_face->GetNode(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl
                  << "---------------------------------------------------------" << std::endl;
    }

    std::cout << "***************************************************************" << std::endl;
    //Now printing all the nodes
    for (unsigned i = 0; i < pElement->GetNumNodes(); ++i)
    {
        Node<3>* p_node = pElement->GetNode(i);
        std::cout << "NODE (" << i << ") : " << p_node->GetIndex() << (p_node->IsDeleted() ? " (DELETED)" : "") << std::endl;
        std::cout << TAB << "Node Attribute : " << Monolayer::ValueToString[GetNodeType(p_node)] << (p_node->IsBoundaryNode() ? " (BOUNDARY)" : "") << std::endl;

        std::set<unsigned> set_tmp = p_node->rGetContainingElementIndices();
        std::cout << TAB << "number of Elements : " << set_tmp.size() << " {  ";
        for (unsigned elem_index : set_tmp)
        {
            std::cout << elem_index << "  ";
        }
        std::cout << "}" << std::endl;

        std::set<unsigned> face_tmp = p_node->rGetContainingFaceIndices();
        std::cout << TAB << "number of Faces : " << face_tmp.size() << " {  ";
        for (unsigned face_index : face_tmp)
        {
            std::cout << face_index << "  ";
        }
        std::cout << "}" << std::endl
                  << "---------------------------------------------------------" << std::endl;
    }
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
template <>
void PrintElement(const VertexElement<2, 3>* pFace)
{
    const std::string TAB = "    ";

    std::cout << "=================================================================================" << std::endl;

    std::cout << "FACE: " << pFace->GetIndex() << " " << Monolayer::ValueToString[GetFaceType(pFace)]
              << (pFace->IsDeleted() ? " (DELETED)" : "") << std::endl;

    const std::set<unsigned>& s_elems = (const_cast<VertexElement<2, 3>*>(pFace))->rFaceGetContainingElementIndices();
    std::cout << TAB << "number of Elements : " << pFace->FaceGetNumContainingElements() << " {";
    for (unsigned elem_index : s_elems)
    {
        std::cout << std::setw(3) << elem_index << "  ";
    }
    std::cout << "}" << std::endl;

    std::cout << TAB << "number of Nodes : " << pFace->GetNumNodes() << " {  ";
    for (unsigned j = 0; j < pFace->GetNumNodes(); ++j)
    {
        std::cout << pFace->GetNode(j)->GetIndex() << "  ";
    }
    std::cout << "}" << std::endl;

    std::cout << "=================================================================================" << std::endl;
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
void PrintMesh(const VertexMesh<3, 3>* pMesh, const bool printDeletedObjects)
{
    const std::string TAB = "    ";

    std::cout << "=================================================================================" << std::endl;
    // Printing out each elements
    const unsigned num_elems = printDeletedObjects ? pMesh->GetNumAllElements() : pMesh->GetNumElements();
    for (unsigned i = 0; i < num_elems; ++i)
    {
        VertexElement<3, 3>& elem = *(pMesh->GetElement(i));
        std::cout << "ELEMENT (" << i << ") : " << elem.GetIndex() << (elem.IsDeleted() ? " (DELETED)" : "") << std::endl;
        std::cout << TAB << "number of Faces : " << elem.GetNumFaces() << " {";
        for (unsigned j = 0; j < elem.GetNumFaces(); ++j)
        {
            std::cout << std::setw(3) << elem.GetFace(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl;

        std::cout << TAB << "Face oriented.. : " << elem.GetNumFaces() << " {";
        for (unsigned j = 0; j < elem.GetNumFaces(); ++j)
        {
            std::cout << std::setw(3) << elem.FaceIsOrientatedAntiClockwise(j) << "  ";
        }
        std::cout << "}" << std::endl;

        std::cout << TAB << "number of Nodes : " << elem.GetNumNodes() << " {  ";
        for (unsigned j = 0; j < elem.GetNumNodes(); ++j)
        {
            std::cout << elem.GetNode(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl
                  << "---------------------------------------------------------" << std::endl;
    }

    std::cout << "***************************************************************" << std::endl;
    // Now printing all faces
    const unsigned num_faces = printDeletedObjects ? pMesh->GetNumAllFaces() : pMesh->GetNumFaces();
    for (unsigned i = 0; i < num_faces; ++i)
    {
        VertexElement<2, 3>& face = *(pMesh->GetFace(i));
        std::cout << "FACE (" << i << ") : " << face.GetIndex() << (face.IsDeleted() ? " (DELETED)" : "") << std::endl;
        std::cout << TAB << "Face Attribute : " << Monolayer::ValueToString[GetFaceType(&face)] << (IsFaceOnBoundary(&face) ? " (BOUNDARY)" : "") << std::endl;

        std::set<unsigned> set_tmp = face.rFaceGetContainingElementIndices();
        std::cout << TAB << "number of Elements : " << set_tmp.size() << " {  ";
        for (unsigned elem_index : set_tmp)
        {
            std::cout << elem_index << "  ";
        }
        std::cout << "}" << std::endl;

        std::cout << TAB << "number of Nodes : " << face.GetNumNodes() << " {  ";
        for (unsigned j = 0; j < face.GetNumNodes(); ++j)
        {
            std::cout << face.GetNode(j)->GetIndex() << "  ";
        }
        std::cout << "}" << std::endl
                  << "---------------------------------------------------------" << std::endl;
    }

    std::cout << "***************************************************************" << std::endl;
    //Now printing all the nodes
    const unsigned num_nodes = printDeletedObjects ? pMesh->GetNumAllNodes() : pMesh->GetNumNodes();
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        Node<3>* p_node = pMesh->GetNode(i);
        std::cout << "NODE (" << i << ") : " << p_node->GetIndex() << (p_node->IsDeleted() ? " (DELETED)" : "") << std::endl;
        std::cout << TAB << "Node Attribute : " << Monolayer::ValueToString[GetNodeType(p_node)] << (p_node->IsBoundaryNode() ? " (BOUNDARY)" : "") << std::endl;

        std::set<unsigned> set_tmp = p_node->rGetContainingElementIndices();
        std::cout << TAB << "number of Elements : " << set_tmp.size() << " {  ";
        for (unsigned elem_index : set_tmp)
        {
            std::cout << elem_index << "  ";
        }
        std::cout << "}" << std::endl;

        std::set<unsigned> face_tmp = p_node->rGetContainingFaceIndices();
        std::cout << TAB << "number of Faces : " << face_tmp.size() << " {  ";
        for (unsigned face_index : face_tmp)
        {
            std::cout << face_index << "  ";
        }
        std::cout << "}" << std::endl
                  << "---------------------------------------------------------" << std::endl;
    }
}
// LCOV_EXCL_STOP

bool IsFaceOnBoundary(const VertexElement<2, 3>* pFace)
{
    return pFace->FaceGetNumContainingElements() == 1;
}

/**
 * Compute the cyclic node order of a monolayer element's apical or basal face from
 * the element's lateral faces, instead of sorting the nodes geometrically by angle.
 *
 * Each lateral face of a prism is a quad sharing exactly two nodes with the apical/basal face (its
 * two nodes of the same type); those two nodes are adjacent on the face. Collecting these
 * adjacencies over all lateral faces yields the face's node ring, independent of geometry. This is
 * robust on curved (e.g. cylindrical) meshes and after division makes a face non-convex, where the
 * geometric angular sort in VertexElement<2,3>::FaceRearrangeNodes() can transpose two nodes and so
 * break the prism invariant that every apical/basal edge has a matching lateral face.
 *
 * Returns false (so the caller can fall back to the geometric sort) if the lateral faces do not
 * form a single clean ring over this face's nodes - e.g. when an asynchronous T1 swap has left a
 * scutoid interface node, so a lateral face is a triangle sharing only one node with the face.
 *
 * @param pElement  the (single) element owning this apical/basal face
 * @param pFace  the apical or basal face to be ordered
 * @param rOrderedNodes  filled with the topological node order on success
 * @return whether a clean ring was found
 */
static bool ComputeApicalBasalNodeOrderFromLateralFaces(const VertexElement<3, 3>* pElement,
                                                        const VertexElement<2, 3>* pFace,
                                                        std::vector<Node<3>*>& rOrderedNodes)
{
    const unsigned num_nodes = pFace->GetNumNodes();
    if (num_nodes < 3)
    {
        return false;
    }

    std::vector<Node<3>*> face_node_vec(num_nodes);
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        face_node_vec[i] = pFace->GetNode(i);
    }
    const std::set<Node<3>*> face_node_set(face_node_vec.begin(), face_node_vec.end());

    // Build the adjacency among this face's nodes from the element's lateral faces.
    std::map<Node<3>*, std::set<Node<3>*> > adjacency;
    for (unsigned f = 0; f < pElement->GetNumFaces(); ++f)
    {
        const VertexElement<2, 3>* p_lateral = pElement->GetFace(f);
        if (!IsLateralFace(p_lateral))
        {
            continue;
        }
        std::vector<Node<3>*> shared;
        for (unsigned j = 0; j < p_lateral->GetNumNodes(); ++j)
        {
            Node<3>* p_node = p_lateral->GetNode(j);
            if (face_node_set.count(p_node) == 1u)
            {
                shared.push_back(p_node);
            }
        }
        if (shared.size() != 2u)
        {
            // Not a clean prism (e.g. a scutoid interface node); let the caller fall back.
            return false;
        }
        adjacency[shared[0]].insert(shared[1]);
        adjacency[shared[1]].insert(shared[0]);
    }

    // For a single closed ring every node must have exactly two neighbours.
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        if (adjacency[face_node_vec[i]].size() != 2u)
        {
            return false;
        }
    }

    // Walk the ring, starting from node 0 and never stepping back to the previous node.
    rOrderedNodes.clear();
    rOrderedNodes.reserve(num_nodes);
    Node<3>* p_start = face_node_vec[0];
    Node<3>* p_prev = nullptr;
    Node<3>* p_curr = p_start;
    do
    {
        rOrderedNodes.push_back(p_curr);
        Node<3>* p_next = nullptr;
        for (Node<3>* p_candidate : adjacency[p_curr])
        {
            if (p_candidate != p_prev)
            {
                p_next = p_candidate;
                break;
            }
        }
        p_prev = p_curr;
        p_curr = p_next;
    } while (p_curr != p_start && p_curr != nullptr && rOrderedNodes.size() <= num_nodes);

    // Every node must have been visited exactly once and the loop must have closed.
    if (rOrderedNodes.size() != num_nodes || p_curr != p_start)
    {
        return false;
    }

    /*
     * The walk direction is arbitrary, so fix the winding to match the historical convention used by
     * the geometric sort (VertexElement<2,3>::FaceRearrangeNodes): nodes ordered counterclockwise as
     * viewed from outside the element, i.e. the polygon's Newell normal points away from the element
     * centroid. This keeps the node order (and thus serialised meshes and any code that relies on the
     * winding) unchanged apart from repairing transposed nodes.
     */
    const c_vector<double, 3> face_centroid = pFace->GetCentroid();
    const c_vector<double, 3> outward_direction = face_centroid - pElement->GetCentroid();
    c_vector<double, 3> newell_normal = zero_vector<double>(3);
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        const c_vector<double, 3> a = rOrderedNodes[i]->rGetLocation() - face_centroid;
        const c_vector<double, 3> b = rOrderedNodes[(i + 1) % num_nodes]->rGetLocation() - face_centroid;
        // Accumulate the cross product a x b (Newell's method) without depending on VectorProduct,
        // which is only included later in this translation unit.
        newell_normal[0] += a[1] * b[2] - a[2] * b[1];
        newell_normal[1] += a[2] * b[0] - a[0] * b[2];
        newell_normal[2] += a[0] * b[1] - a[1] * b[0];
    }
    if (inner_prod(newell_normal, outward_direction) < 0.0)
    {
        std::reverse(rOrderedNodes.begin(), rOrderedNodes.end());
    }

    return true;
}

/**
 * Whether a lateral quad's two basal nodes sit at cyclically consecutive positions
 * (so the face reads basal-basal-apical-apical). The geometric angular sort can instead interleave
 * them (basal-apical-basal-apical) on non-planar quads, which would trip the assertion in
 * VertexElement<2,3>::LateralFaceRearrangeNodes(). Returns true for anything that is not a simple
 * two-basal/two-apical quad, leaving the existing handling to deal with it.
 */
static bool AreLateralBasalNodesAdjacent(const VertexElement<2, 3>* pFace)
{
    const unsigned num_nodes = pFace->GetNumNodes();
    std::vector<unsigned> basal_local_indices;
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        if (IsBasalNode(pFace->GetNode(i)))
        {
            basal_local_indices.push_back(i);
        }
    }
    if (basal_local_indices.size() != 2u)
    {
        return true;
    }
    const unsigned gap = basal_local_indices[1] - basal_local_indices[0];
    return (gap == 1u || gap == num_nodes - 1u);
}

/**
 * Compute a lateral quad's node order topologically from the opposite-node (vertical
 * edge) relationship, independent of geometry: basal_a, basal_b, opposite(basal_b), opposite(basal_a).
 * This keeps the two basal nodes adjacent and each basal node next to its apical partner, which is the
 * invariant LateralFaceRearrangeNodes() and the lateral-face consumers expect. Returns false if the
 * face is not a simple two-basal/two-apical quad, or if the opposite nodes are not this face's apical
 * nodes, so the caller can leave the face alone.
 *
 * @param pMesh  the mesh owning the face (needed to resolve opposite nodes)
 * @param pFace  the lateral face
 * @param rOrderedNodes  filled with the topological node order on success
 * @return whether a valid quad order was produced
 */
static bool ComputeLateralFaceNodeOrder(VertexMesh<3, 3>* pMesh, const VertexElement<2, 3>* pFace, std::vector<Node<3>*>& rOrderedNodes)
{
    const std::vector<Node<3>*> basal_nodes = GetNodesWithType(pFace, Monolayer::BasalValue);
    const std::vector<Node<3>*> apical_nodes = GetNodesWithType(pFace, Monolayer::ApicalValue);
    if (basal_nodes.size() != 2u || apical_nodes.size() != 2u)
    {
        return false;
    }
    Node<3>* p_basal_a = basal_nodes[0];
    Node<3>* p_basal_b = basal_nodes[1];
    Node<3>* p_apical_a = GetOppositeNode(p_basal_a, pMesh);
    Node<3>* p_apical_b = GetOppositeNode(p_basal_b, pMesh);
    const bool matches = (p_apical_a == apical_nodes[0] && p_apical_b == apical_nodes[1])
        || (p_apical_a == apical_nodes[1] && p_apical_b == apical_nodes[0]);
    if (!matches)
    {
        return false;
    }
    rOrderedNodes.clear();
    rOrderedNodes.push_back(p_basal_a);
    rOrderedNodes.push_back(p_basal_b);
    rOrderedNodes.push_back(p_apical_b);
    rOrderedNodes.push_back(p_apical_a);
    return true;
}

void FaceRearrangeNodesInMesh(VertexMesh<3, 3>* pMesh, VertexElement<2, 3>* pFace)
{
    const std::set<unsigned>& set_tmp = pFace->rFaceGetContainingElementIndices();
    const c_vector<double, 3> centroid_tmp = pMesh->GetElement(no1(set_tmp))->GetCentroid();

    bool has_changes = pFace->FaceRearrangeNodes(centroid_tmp);

    /*
     * The geometric angular sort can interleave a lateral quad's basal and apical nodes on
     * non-planar (e.g. cylindrical) faces, which would trip the assertion in LateralFaceRearrangeNodes.
     * When that has happened, re-derive the node order topologically from the opposite-node relationship.
     * Faces the geometric sort ordered consistently (e.g. on flat meshes) are left untouched.
     */
    if (IsLateralFace(pFace) && !AreLateralBasalNodesAdjacent(pFace))
    {
        std::vector<Node<3>*> ordered_nodes;
        if (ComputeLateralFaceNodeOrder(pMesh, pFace, ordered_nodes))
        {
            has_changes = pFace->FaceSetNodeOrder(ordered_nodes) || has_changes;
        }
    }

    if (has_changes)
    {
        for (unsigned elem_index : set_tmp)
        {
            VertexElement<3, 3>* p_elem = pMesh->GetElement(elem_index);
            p_elem->CheckFaceOrientationOfElement(p_elem->GetFaceLocalIndex(pFace->GetIndex()));
        }

        if (IsLateralFace(pFace))
        {
            pFace->LateralFaceRearrangeNodes();
        }
    }
}

/**
 * Whether an apical/basal face's stored node order is consistent with the element's
 * lateral faces: every consecutive pair of nodes must share exactly one lateral face. This is the
 * prism invariant that DivideElementAlongGivenAxis() (via GetSharedLateralFace()) relies on.
 */
static bool IsApicalBasalFaceConsistent(const VertexElement<3, 3>* pElement, const VertexElement<2, 3>* pFace)
{
    const unsigned num_nodes = pFace->GetNumNodes();
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        const std::set<VertexElement<2, 3>*> shared_lateral = GetFacesWithIndices(
            GetSharedFaceIndices(pFace->GetNode(i), pFace->GetNode((i + 1) % num_nodes)),
            pElement, Monolayer::LateralValue);
        if (shared_lateral.size() != 1u)
        {
            return false;
        }
    }
    return true;
}

void RepairApicalBasalFaceOrdering(VertexMesh<3, 3>* pMesh, VertexElement<2, 3>* pFace)
{
    if (!(IsApicalFace(pFace) || IsBasalFace(pFace)))
    {
        return;
    }
    const std::set<unsigned>& set_tmp = pFace->rFaceGetContainingElementIndices();
    if (set_tmp.empty())
    {
        return;
    }
    VertexElement<3, 3>* p_elem = pMesh->GetElement(no1(set_tmp));

    /*
     * The geometric angular sort in VertexElement<2,3>::FaceRearrangeNodes() is correct for planar,
     * convex faces but can transpose nodes on curved or non-convex apical/basal faces, breaking the
     * prism invariant. Only in that case do we re-derive the node order topologically from the
     * lateral faces; consistent faces are left exactly as they are, so flat/convex meshes (and their
     * stored reference outputs) are unaffected.
     */
    if (IsApicalBasalFaceConsistent(p_elem, pFace))
    {
        return;
    }
    std::vector<Node<3>*> ordered_nodes;
    if (ComputeApicalBasalNodeOrderFromLateralFaces(p_elem, pFace, ordered_nodes))
    {
        if (pFace->FaceSetNodeOrder(ordered_nodes))
        {
            p_elem->CheckFaceOrientationOfElement(p_elem->GetFaceLocalIndex(pFace->GetIndex()));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////
///                       Functions for monolayer classes                       ///
///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
///     Functions for nodes     ///
///////////////////////////////////

void SetNodeAsApical(Node<3>* pNode)
{
    assert(pNode->GetNumNodeAttributes() == 0u); // LCOV_EXCL_LINE
    pNode->AddNodeAttribute(Monolayer::SetApicalValue);
}

void SetNodeAsBasal(Node<3>* pNode)
{
    assert(pNode->GetNumNodeAttributes() == 0u); // LCOV_EXCL_LINE
    pNode->AddNodeAttribute(Monolayer::SetBasalValue);
}

void SetNodeAsLateral(Node<3>* pNode)
{
    assert(pNode->GetNumNodeAttributes() == 0u); // LCOV_EXCL_LINE
    pNode->AddNodeAttribute(Monolayer::SetLateralValue);
}

Monolayer::v_type GetNodeType(const Node<3>* pNode)
{
    /* implemented const node as there is no modication for this
     * function. However, there isn't suitable const function for
     * NodeAttribute and hence required a const_cast
     */
    Node<3>* p_non_const_node = const_cast<Node<3>*>(pNode);
    switch (p_non_const_node->GetNumNodeAttributes())
    {
        case 0:
        {
            return 0;
        }
        case 1:
        {
            return static_cast<Monolayer::v_type>(p_non_const_node->rGetNodeAttributes()[0]);
        }
        default:
            NEVER_REACHED;
    }
}

bool IsApicalNode(const Node<3>* pNode)
{
    return GetNodeType(pNode) == Monolayer::ApicalValue;
}

bool IsBasalNode(const Node<3>* pNode)
{
    return GetNodeType(pNode) == Monolayer::BasalValue;
}

bool IsLateralNode(const Node<3>* pNode)
{
    return GetNodeType(pNode) == Monolayer::LateralValue;
}

//////////////////////////////////
///     Functions for face     ///
//////////////////////////////////

void SetFaceAsApical(VertexElement<2, 3>* pFace)
{
    assert(pFace->GetNumElementAttributes() == 0u); // LCOV_EXCL_LINE
    pFace->AddElementAttribute(Monolayer::SetApicalValue);
}

void SetFaceAsBasal(VertexElement<2, 3>* pFace)
{
    assert(pFace->GetNumElementAttributes() == 0u); // LCOV_EXCL_LINE
    pFace->AddElementAttribute(Monolayer::SetBasalValue);
}

void SetFaceAsLateral(VertexElement<2, 3>* pFace)
{
    assert(pFace->GetNumElementAttributes() == 0u); // LCOV_EXCL_LINE
    pFace->AddElementAttribute(Monolayer::SetLateralValue);
}

Monolayer::v_type GetFaceType(const VertexElement<2, 3>* pFace)
{
    /* implemented const node as there is no modication for this
     * function. However, there isn't suitable const function for
     * ElementAttribute and hence required a const_cast
     */
    VertexElement<2, 3>* p_non_const_face = const_cast<VertexElement<2, 3>*>(pFace);
    switch (p_non_const_face->GetNumElementAttributes())
    {
        case 0:
        {
            return 0;
        }
        case 1:
        {
            return static_cast<Monolayer::v_type>(p_non_const_face->rGetElementAttributes()[0]);
        }
        default:
            NEVER_REACHED;
    }
}

bool IsApicalFace(const VertexElement<2, 3>* pFace)
{
    return GetFaceType(pFace) == Monolayer::ApicalValue;
}

bool IsBasalFace(const VertexElement<2, 3>* pFace)
{
    return GetFaceType(pFace) == Monolayer::BasalValue;
}

bool IsLateralFace(const VertexElement<2, 3>* pFace)
{
    return GetFaceType(pFace) == Monolayer::LateralValue;
}

////////////////////////////////////////////
///     Functions for Vertex Element     ///
////////////////////////////////////////////

void SetElementAsMonolayer(VertexElement<3, 3>* pElement)
{
    assert(pElement->GetNumElementAttributes() == 0u); // LCOV_EXCL_LINE
    pElement->AddElementAttribute(Monolayer::SetElementValue);

    for (unsigned i = 0; i < pElement->GetNumNodes(); ++i)
    {
        if (pElement->GetNode(i)->GetNumNodeAttributes() == 0u)
        {
            NEVER_REACHED;
        }
    }

    for (unsigned i = 0; i < pElement->GetNumFaces(); ++i)
    {
        if (pElement->GetFace(i)->GetNumElementAttributes() == 0u)
        {
            VertexElement<2, 3>* p_face = pElement->GetFace(i);
            std::set<Monolayer::v_type> node_vals;
            for (unsigned j = 0; j < p_face->GetNumNodes(); ++j)
            {
                node_vals.insert(GetNodeType(p_face->GetNode(j)));
            }

            if (node_vals.size() == 2)
            {
                SetFaceAsLateral(p_face);
            }
            else
            {
                assert(node_vals.size() == 1);
                switch (*(node_vals.begin()))
                {
                    case Monolayer::ApicalValue:
                        SetFaceAsApical(p_face);
                        break;
                    case Monolayer::BasalValue:
                        SetFaceAsBasal(p_face);
                        break;
                    default:
                        NEVER_REACHED;
                }
            }
        }
    }

    pElement->MonolayerElementRearrangeFacesNodes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool IsMonolayerElement(const VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    if (ELEMENT_DIM != 3 || SPACE_DIM != 3)
    {
        return false;
    }
    else
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_non_const_elem = const_cast<VertexElement<ELEMENT_DIM, SPACE_DIM>*>(pElement);
        return (p_non_const_elem->GetNumElementAttributes() == 1u) && (Monolayer::v_type(p_non_const_elem->rGetElementAttributes()[0]) == Monolayer::ElementValue);
    }
}
// Template instantiation is for the compiler
template bool IsMonolayerElement<1, 1>(const VertexElement<1, 1>* pElement);
template bool IsMonolayerElement<1, 2>(const VertexElement<1, 2>* pElement);
template bool IsMonolayerElement<1, 3>(const VertexElement<1, 3>* pElement);
template bool IsMonolayerElement<2, 2>(const VertexElement<2, 2>* pElement);
template bool IsMonolayerElement<2, 3>(const VertexElement<2, 3>* pElement);
template bool IsMonolayerElement<3, 3>(const VertexElement<3, 3>* pElement);

VertexElement<2, 3>* GetApicalFace(const VertexElement<3, 3>* pElement)
{
    VertexElement<2, 3>* p_face = pElement->GetFace(1);
    if (!IsApicalFace(p_face))
    {
        for (unsigned i = 0; i < pElement->GetNumFaces(); ++i)
        {
            if (IsApicalFace(pElement->GetFace(i)))
            {
                p_face = pElement->GetFace(i);
                break;
            }
        }
    }

    assert(p_face != nullptr && IsApicalFace(p_face)); // LCOV_EXCL_LINE
    return p_face;
}

VertexElement<2, 3>* GetBasalFace(const VertexElement<3, 3>* pElement)
{
    VertexElement<2, 3>* p_face(nullptr);
    for (unsigned i = 0; i < pElement->GetNumFaces(); ++i)
    {
        if (IsBasalFace(pElement->GetFace(i)))
        {
            p_face = pElement->GetFace(i);
            break;
        }
    }

    assert(p_face != nullptr && IsBasalFace(p_face)); // LCOV_EXCL_LINE
    return p_face;
}

template <unsigned ELEMENT_DIM>
std::vector<Node<3>*> GetNodesWithType(const VertexElement<ELEMENT_DIM, 3>* pElement, const Monolayer::v_type nodeType)
{
    std::vector<Node<3>*> return_v;
    for (unsigned i = 0; i < pElement->GetNumNodes(); ++i)
    {
        if (GetNodeType(pElement->GetNode(i)) == nodeType)
        {
            return_v.push_back(pElement->GetNode(i));
        }
    }

    return return_v;
}
template std::vector<Node<3>*> GetNodesWithType(const VertexElement<3, 3>* pElement, const Monolayer::v_type nodeType);
template std::vector<Node<3>*> GetNodesWithType(const VertexElement<2, 3>* pFace, const Monolayer::v_type nodeType);

std::vector<VertexElement<2, 3>*> GetFacesWithType(const VertexElement<3, 3>* pElement, const Monolayer::v_type faceType)
{
    std::vector<VertexElement<2, 3>*> return_v;
    for (unsigned i = 0; i < pElement->GetNumFaces(); ++i)
    {
        if (GetFaceType(pElement->GetFace(i)) == faceType)
        {
            return_v.push_back(pElement->GetFace(i));
        }
    }

    return return_v;
}

// Having overloaded functions so that Opposite node can be obtained in almost any scope.
Node<3>* GetOppositeNode(const Node<3>* pNode, const VertexElement<2, 3>* pFace)
{
    if (!(IsApicalNode(pNode) || IsBasalNode(pNode)))
    {
        EXCEPTION("No Opposite Node"); // LCOV_EXCL_LINE
    }

    const Monolayer::v_type opposite_type = IsApicalNode(pNode) ? Monolayer::BasalValue : Monolayer::ApicalValue;

    /*
     * #480 The opposite node is the node of the opposite type (apical<->basal) joined to pNode by a
     * vertical (lateral) edge. A basal node and an apical node can share only lateral faces, and the
     * opposite node lies on *every* lateral face containing pNode, whereas any other node of the
     * opposite type lies on strictly fewer. So the opposite node is the opposite-type node of pFace
     * that shares the most faces with pNode.
     *
     * This replaces requiring an exact match against a lateral-face count inferred from
     * (num containing faces - num containing elements): that count assumed exactly one basal/apical
     * face per containing element and so broke near interface (scutoid) nodes, raising "No Opposite
     * Node". Taking the maximum needs no such assumption and degrades gracefully.
     */
    Node<3>* p_opposite_node = nullptr;
    unsigned max_num_shared_faces = 0;
    for (unsigned i = 0; i < pFace->GetNumNodes(); ++i)
    {
        Node<3>* p_candidate_node = pFace->GetNode(i);
        if (GetNodeType(p_candidate_node) != opposite_type)
        {
            continue;
        }

        const unsigned num_shared_faces = GetSharedFaceIndices(pNode, p_candidate_node).size();
        if (num_shared_faces > max_num_shared_faces)
        {
            max_num_shared_faces = num_shared_faces;
            p_opposite_node = p_candidate_node;
        }
    }

    if (p_opposite_node == nullptr)
    {
        EXCEPTION("No Opposite Node"); // LCOV_EXCL_LINE
    }
    return p_opposite_node;
}

Node<3>* GetOppositeNode(const Node<3>* pNode, const VertexElement<3, 3>* pElement)
{
    VertexElement<2, 3>* p_other_face = nullptr;
    switch (GetNodeType(pNode))
    {
        case Monolayer::LateralValue:
            EXCEPTION("No Opposite Node"); // LCOV_EXCL_LINE
        case Monolayer::BasalValue:
            p_other_face = GetApicalFace(pElement);
            break;
        case Monolayer::ApicalValue:
            p_other_face = GetBasalFace(pElement);
            break;
        default:
            NEVER_REACHED;
    }

    return GetOppositeNode(pNode, p_other_face);
}

Node<3>* GetOppositeNode(const Node<3>* pNode, const VertexMesh<3, 3>* pMesh)
{
    const std::set<unsigned>& s_containing_elems = const_cast<Node<3>*>(pNode)->rGetContainingElementIndices();
    return GetOppositeNode(pNode, pMesh->GetElement(no1(s_containing_elems)));
}

void MeshUpdateNode(Node<3>* pOldNode, Node<3>* pNewNode, MutableVertexMesh<3, 3>* pMesh)
{
    const std::set<unsigned> containing_faces = pOldNode->rGetContainingFaceIndices();
    const std::set<unsigned> containing_elems = pOldNode->rGetContainingElementIndices();
    for (unsigned face_index : containing_faces)
    {
        if (pMesh->GetFace(face_index)->GetNodeLocalIndex(pNewNode->GetIndex()) != UINT_MAX)
        {
            pMesh->GetFace(face_index)->FaceDeleteNode(pOldNode);
        }
        else
        {
            pMesh->GetFace(face_index)->FaceUpdateNode(pOldNode, pNewNode);
        }
    }

    for (unsigned elem_index : containing_elems)
    {
        if (pMesh->GetElement(elem_index)->GetNodeLocalIndex(pNewNode->GetIndex()) != UINT_MAX)
        {
            pMesh->GetElement(elem_index)->DeleteNode(pOldNode);
        }
        else
        {
            pMesh->GetElement(elem_index)->UpdateNode(pOldNode, pNewNode);
        }
    }

    pMesh->DeleteNodePriorToReMesh(pOldNode->GetIndex());
}

template <unsigned ELEMENT_DIM>
Node<3>* GetNextNode(const unsigned localIndex, const VertexElement<ELEMENT_DIM, 3>* pElement)
{
    const unsigned num_nodes = pElement->GetNumNodes();
    return pElement->GetNode(plus1(localIndex, num_nodes));
}
template Node<3>* GetNextNode(const unsigned, const VertexElement<2, 3>*);
template Node<3>* GetNextNode(const unsigned, const VertexElement<3, 3>*);

template <unsigned ELEMENT_DIM>
Node<3>* GetNextNode(const Node<3>* pNode, const VertexElement<ELEMENT_DIM, 3>* pElement)
{
    const unsigned current_index = pElement->GetNodeLocalIndex(pNode->GetIndex());
    return GetNextNode(current_index, pElement);
}
template Node<3>* GetNextNode(const Node<3>*, const VertexElement<2, 3>*);
template Node<3>* GetNextNode(const Node<3>*, const VertexElement<3, 3>*);

#include "UblasCustomFunctions.hpp"

c_vector<double, 3> CalculateUnitNormalToFace(const VertexElement<2, 3>* pFace)
{
    const unsigned num_nodes = pFace->GetNumNodes();
    // A face in 3D needs at least three vertices to define a plane.
    assert(num_nodes >= 3u); // LCOV_EXCL_LINE

    /*
     * #480 Newell's method: accumulate the cross products of consecutive edge vectors, taken relative
     * to node 0 for numerical robustness. Unlike the previous least-squares plane fit, this is
     * well-conditioned for axis-aligned faces (e.g. every node at z = 0, which made that fit singular)
     * and is correct for non-planar and non-convex faces. A degenerate (zero-area) face yields the
     * zero vector.
     */
    const c_vector<double, 3> ref = pFace->GetNodeLocation(0);
    c_vector<double, 3> normal = zero_vector<double>(3);
    for (unsigned i = 0; i < num_nodes; ++i)
    {
        const c_vector<double, 3> v_this = pFace->GetNodeLocation(i) - ref;
        const c_vector<double, 3> v_next = pFace->GetNodeLocation((i + 1) % num_nodes) - ref;
        normal += VectorProduct(v_this, v_next);
    }

    const double magnitude = norm_2(normal);
    if (magnitude > 1e-12)
    {
        normal /= magnitude;
    }
    return normal;
}
