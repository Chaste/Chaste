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

#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "Node.hpp"
#include "VertexElement.hpp"
#include "MutableVertexMesh.hpp"

#include <algorithm>
#include <string>

#include <iostream> // PrintMesh
#include <iostream> // PrintElement
#include <iomanip> // PrintElement

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
    if (pElemA == NULL || pElemB == NULL)
    {
        EXCEPTION("Two elements do not share lateral face.");
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
            EXCEPTION("Two elements do not share lateral face.");
        case 1:
            return pElemA->GetFace(pElemA->GetFaceLocalIndex(*(s_return.begin())));
        default:
            EXCEPTION("Probably some errors occur. Two elements share more than 1 lateral face.");
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
            EXCEPTION("Two nodes do not share lateral face.");
        case 1:
            return no1(shared_lateral_faces);
        default:
            EXCEPTION("Probably some errors occur. Two nodes share more than 1 lateral face in monolayer mesh.");
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
        for (std::set<VertexElement<2, 3>*>::iterator it = s_result.begin(); it != s_result.end(); ++it)
        {
            if (GetFaceType(*it) == faceType)
            {
                s_matched.insert(*it);
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
    for (std::set<unsigned>::const_iterator it = face_indices.begin(); it != face_indices.end(); ++it)
    {
        s_result.insert(pMesh->GetFace(*it));
    }

    if (faceType != Monolayer::AllTypes)
    {
        std::set<VertexElement<2, 3>*> s_matched;
        for (std::set<VertexElement<2, 3>*>::iterator it = s_result.begin(); it != s_result.end(); ++it)
        {
            if (GetFaceType(*it) == faceType)
            {
                s_matched.insert(*it);
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
        for (std::set<unsigned>::iterator it = set_tmp.begin(); it != set_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
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
        for (std::set<unsigned>::iterator it = set_tmp.begin(); it != set_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
        }
        std::cout << "}" << std::endl;

        std::set<unsigned> face_tmp = p_node->rGetContainingFaceIndices();
        std::cout << TAB << "number of Faces : " << face_tmp.size() << " {  ";
        for (std::set<unsigned>::iterator it = face_tmp.begin(); it != face_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
        }
        std::cout << "}" << std::endl
                  << "---------------------------------------------------------" << std::endl;
    }
}

template <>
void PrintElement(const VertexElement<2, 3>* pFace)
{
    const std::string TAB = "    ";

    std::cout << "=================================================================================" << std::endl;

    std::cout << "FACE: " << pFace->GetIndex() << " " << Monolayer::ValueToString[GetFaceType(pFace)]
              << (pFace->IsDeleted() ? " (DELETED)" : "") << std::endl;

    const std::set<unsigned>& s_elems = (const_cast<VertexElement<2, 3>*>(pFace))->rFaceGetContainingElementIndices();
    std::cout << TAB << "number of Elements : " << pFace->FaceGetNumContainingElements() << " {";
    for (std::set<unsigned>::const_iterator it = s_elems.begin();
         it != s_elems.end(); ++it)
    {
        std::cout << std::setw(3) << *it << "  ";
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
        std::cout << TAB << "Face Attrbute : " << Monolayer::ValueToString[GetFaceType(&face)] << (IsFaceOnBoundary(&face) ? " (BOUNDARY)" : "") << std::endl;

        std::set<unsigned> set_tmp = face.rFaceGetContainingElementIndices();
        std::cout << TAB << "number of Elements : " << set_tmp.size() << " {  ";
        for (std::set<unsigned>::iterator it = set_tmp.begin(); it != set_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
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
        for (std::set<unsigned>::iterator it = set_tmp.begin(); it != set_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
        }
        std::cout << "}" << std::endl;

        std::set<unsigned> face_tmp = p_node->rGetContainingFaceIndices();
        std::cout << TAB << "number of Faces : " << face_tmp.size() << " {  ";
        for (std::set<unsigned>::iterator it = face_tmp.begin(); it != face_tmp.end(); ++it)
        {
            std::cout << *it << "  ";
        }
        std::cout << "}" << std::endl
                  << "---------------------------------------------------------" << std::endl;
    }
}

bool IsFaceOnBoundary(const VertexElement<2, 3>* pFace)
{
    return pFace->FaceGetNumContainingElements() == 1;
}

void FaceRearrangeNodesInMesh(VertexMesh<3, 3>* pMesh, VertexElement<2, 3>* pFace)
{
    const std::set<unsigned>& set_tmp = pFace->rFaceGetContainingElementIndices();
    const c_vector<double, 3> centroid_tmp = pMesh->GetElement(no1(set_tmp))->GetCentroid();

    const bool has_changes = pFace->FaceRearrangeNodes(centroid_tmp);

    if (has_changes)
    {
        for (std::set<unsigned>::const_iterator it = set_tmp.begin(); it != set_tmp.end(); ++it)
        {
            VertexElement<3, 3> *p_elem = pMesh->GetElement(*it);
            p_elem->CheckFaceOrientationOfElement(p_elem->GetFaceLocalIndex(pFace->GetIndex()));
        }

        if (IsLateralFace(pFace))
        {
            pFace->LateralFaceRearrangeNodes();
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

template <unsigned DIM>
Monolayer::v_type GetNodeType(const Node<DIM>* pNode)
{
    if (DIM != 3u)
    {
        NEVER_REACHED;
    }

    /* implemented const node as there is no modication for this
     * function. However, there isn't suitable const function for
     * NodeAttribute and hence required a const_cast
     */
    Node<DIM>* p_non_const_node = const_cast<Node<DIM>*>(pNode);
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
template Monolayer::v_type GetNodeType(const Node<1>* pNode);
template Monolayer::v_type GetNodeType(const Node<2>* pNode);
template Monolayer::v_type GetNodeType(const Node<3>* pNode);

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

    assert(p_face != NULL && IsApicalFace(p_face)); // LCOV_EXCL_LINE
    return p_face;
}

VertexElement<2, 3>* GetBasalFace(const VertexElement<3, 3>* pElement)
{
    VertexElement<2, 3>* p_face(NULL);
    for (unsigned i = 0; i < pElement->GetNumFaces(); ++i)
    {
        if (IsBasalFace(pElement->GetFace(i)))
        {
            p_face = pElement->GetFace(i);
            break;
        }
    }

    assert(p_face != NULL && IsBasalFace(p_face)); // LCOV_EXCL_LINE
    return p_face;
}

unsigned MonolayerGetHalfNumNodes(const VertexElement<3, 3>* pElement)
{
    if (!IsMonolayerElement(pElement))
    {
        NEVER_REACHED;
    }

    unsigned num_nodes = pElement->GetNumNodes();
    if (num_nodes % 2 != 0)
    {
        NEVER_REACHED;
    }
    num_nodes /= 2;

    if (num_nodes != GetApicalFace(pElement)->GetNumNodes()
        || num_nodes != GetBasalFace(pElement)->GetNumNodes())
    {
        NEVER_REACHED;
    }

    return num_nodes;
}

///\todo #2850 Remove this method?
std::vector<unsigned> GetLateralFace(const VertexElement<3, 3>* pElement, const unsigned nodeIndexA, const unsigned nodeIndexB)
{
    ///\todo #2850 think which will be the more effective way to find lateral face
    std::vector<unsigned> return_vector;
    for (unsigned face_local_index = 2; face_local_index < pElement->GetNumFaces(); ++face_local_index)
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
        assert(return_vector.size() == 3); // LCOV_EXCL_LINE
    }
    return return_vector;
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

///\todo: #2850 num_lateral_faces is obtained using heuristical method. If node were to have weak pointer of 
//          face&elements (using Boost or c++11)
// Having overloaded functions so that Opposite node can be obtained in almost any scope.
Node<3>* GetOppositeNode(const Node<3>* pNode, const VertexElement<2, 3>* pFace)
{
    if (!(IsApicalNode(pNode) || IsBasalNode(pNode)))
    {
        EXCEPTION("No Opposite Node");
    }

    const Monolayer::v_type opposite_type = IsApicalNode(pNode) ? Monolayer::BasalValue : Monolayer::ApicalValue;
    const unsigned num_lateral_faces = pNode->GetNumContainingFaces() - pNode->GetNumContainingElements();

    Node<3>* p_return_node = NULL;
    for (unsigned i = 0; i < pFace->GetNumNodes(); ++i)
    {
        const Node<3>* p_tmp_node = pFace->GetNode(i);
        if (GetNodeType(p_tmp_node) != opposite_type)
        {
            continue;
        }

        if (GetSharedFaceIndices(pNode, p_tmp_node).size() == num_lateral_faces)
        {
            return p_return_node = pFace->GetNode(i);
        }
    }

    EXCEPTION("No Opposite Node");
}

Node<3>* GetOppositeNode(const Node<3>* pNode, const VertexElement<3, 3>* pElement)
{
    VertexElement<2, 3>* p_other_face = NULL;
    switch (GetNodeType(pNode))
    {
        case Monolayer::LateralValue:
            EXCEPTION("No Opposite Node");
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