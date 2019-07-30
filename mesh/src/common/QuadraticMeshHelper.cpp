/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "QuadraticMeshHelper.hpp"

#define SEEK_TO_CONTENT(methNameDirect, methNameIncrement, index) \
    if (index > 0u) {                                             \
        if (rMeshReader.IsFileFormatBinary()) {                   \
            rMeshReader.methNameDirect(index - 1u);               \
        } else {                                                  \
            for (unsigned i=0; i<index-1u; ++i) {                 \
                rMeshReader.methNameIncrement();                  \
    } } }

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SeekToBoundaryElement(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                           unsigned boundaryElementIndex)
{
    SEEK_TO_CONTENT(GetFaceData, GetNextFaceData, boundaryElementIndex);
}

template<unsigned DIM>
void QuadraticMeshHelper<DIM>::AddInternalNodesToElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                                          AbstractMeshReader<DIM,DIM>* pMeshReader)
{
    assert(pMesh);
    assert(pMeshReader);

    if (pMesh->GetNumLocalElements() > 0u)
    {
        pMeshReader->Reset();

        // Create a set of element indices we own
        std::set<unsigned> owned_element_indices;
        for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter = pMesh->GetElementIteratorBegin();
             iter != pMesh->GetElementIteratorEnd();
             ++iter)
        {
            owned_element_indices.insert(iter->GetIndex());
        }

        const std::vector<unsigned>& r_node_perm = pMesh->rGetNodePermutation();

        // Add the extra nodes (1 extra node in 1D, 3 in 2D, 6 in 3D) to the element data
        for (typename AbstractMeshReader<DIM,DIM>::ElementIterator iter = pMeshReader->GetElementIteratorBegin(owned_element_indices);
             iter != pMeshReader->GetElementIteratorEnd();
             ++iter)
        {
            std::vector<unsigned> nodes = iter->NodeIndices;
            assert(nodes.size()==(DIM+1)*(DIM+2)/2);
            Element<DIM,DIM>* p_element = pMesh->GetElement(iter.GetIndex());
            assert(p_element->GetNumNodes()==DIM+1); // Element is initially linear

            // Add extra nodes to make it a quad element
            for (unsigned j=DIM+1; j<(DIM+1)*(DIM+2)/2; j++)
            {
                unsigned node_index = nodes[j];
                if (!r_node_perm.empty())
                {
                    node_index = r_node_perm[node_index];
                }
                Node<DIM>* p_node = pMesh->GetNodeOrHaloNode(node_index);
                p_element->AddNode(p_node);
                p_node->AddElement(p_element->GetIndex());
                p_node->MarkAsInternal();
            }
        }
    }
}

template<unsigned DIM>
void QuadraticMeshHelper<DIM>::AddInternalNodesToBoundaryElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                                                  AbstractMeshReader<DIM,DIM>* pMeshReader)
{
    assert(pMesh);
    assert(pMeshReader);
    // We only have boundary elements in 2d or 3d
    if (DIM > 1 && pMesh->GetNumLocalBoundaryElements() > 0u)
    {
        // If the data is on disk our job is easy
        if (pMeshReader->GetOrderOfBoundaryElements() == 2u)
        {
            // The work should have been done in the linear constructor, but let's check
            // that the first face has more than DIM nodes.
            assert((*pMesh->GetBoundaryElementIteratorBegin())->GetNumNodes()==DIM*(DIM+1)/2);
            return;
        }
        else
        {
            AddNodesToBoundaryElements(pMesh, pMeshReader);
        }
    }
}

template<unsigned DIM>
void QuadraticMeshHelper<DIM>::AddNodesToBoundaryElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                                          AbstractMeshReader<DIM,DIM>* pMeshReader)
 {
    // Loop over all boundary elements, find the equivalent face from all
    // the elements, and add the extra nodes to the boundary element
    bool boundary_element_file_has_containing_element_info = false;

    if (pMeshReader)
    {
        boundary_element_file_has_containing_element_info = pMeshReader->GetReadContainingElementOfBoundaryElement();
    }

    if (DIM > 1)
    {
        if (boundary_element_file_has_containing_element_info)
        {
            pMeshReader->Reset();
        }

        ///\todo #1930 Until there is an interator_facade over boundary elements,
        // we may need to skip through the boundary element file searching for containing elements hints.
        // This counter keeps track of our position in the file.
        unsigned next_face_on_file = 0u;

        for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
                 = pMesh->GetBoundaryElementIteratorBegin();
             iter != pMesh->GetBoundaryElementIteratorEnd();
             ++iter)
        {

            // collect the nodes of this boundary element in a set
            std::set<unsigned> boundary_element_node_indices;
            for (unsigned i=0; i<DIM; i++)
            {
                boundary_element_node_indices.insert( (*iter)->GetNodeGlobalIndex(i) );
            }

            bool found_this_boundary_element = false;

            // Loop over elements surrounding this face, then loop over each face of that element, and see if it matches
            // this boundary element.
            // Note, if we know what elem it should be in (boundary_element_file_has_containing_element_info==true)
            // we will reset elem_index immediately (below)
            Node<DIM>* p_representative_node = (*iter)->GetNode(0);
            for (typename Node<DIM>::ContainingElementIterator element_iter = p_representative_node->ContainingElementsBegin();
                 element_iter != p_representative_node->ContainingElementsEnd();
                 ++element_iter)
            {
                unsigned elem_index = *element_iter;

                // We know what elem it should be in (but we'll still check the node indices match in case)
                if (boundary_element_file_has_containing_element_info)
                {
                    unsigned face_index = (*iter)->GetIndex();
                    ///\todo #1930 Once there is an interator_facade then we can do: elem_index = pMeshReader->GetFaceData(face_index).ContainingElement;
                    do
                    {
                        elem_index = pMeshReader->GetNextFaceData().ContainingElement;
                        next_face_on_file++;
                    }
                    while (face_index >= next_face_on_file);
                }

                Element<DIM,DIM>* p_element = pMesh->GetElement(elem_index);

                // For each element, loop over faces (the opposites to a node)
                for (unsigned face=0; face<DIM+1; face++)
                {
                    // Collect the node indices for this face
                    std::set<unsigned> node_indices;
                    for (unsigned local_node_index=0; local_node_index<DIM+1; local_node_index++)
                    {
                        if (local_node_index != face)
                        {
                            node_indices.insert( p_element->GetNodeGlobalIndex(local_node_index) );
                        }
                    }

                    assert(node_indices.size()==DIM);

                    // See if this face matches the boundary element, and add internal nodes if so
                    if (node_indices == boundary_element_node_indices)
                    {
                        QuadraticMeshHelper<DIM>::AddExtraBoundaryNodes(pMesh, *iter, p_element, face);

                        found_this_boundary_element = true;
                        break;
                    }
                }

                // If the containing element info was given, we must have found the face first time
                if (boundary_element_file_has_containing_element_info && !found_this_boundary_element)
                {
                    // LCOV_EXCL_START
                    //std::cout << (*iter)->GetIndex() << " " <<  pMeshReader->GetNextFaceData().ContainingElement << "\n";
                    EXCEPTION("Boundary element " << (*iter)->GetIndex()
                              << "wasn't found in the containing element given for it "
                              << elem_index);
                    // LCOV_EXCL_STOP
                }

                if (found_this_boundary_element)
                {
                    break;
                }
            }

            if (!found_this_boundary_element)
            {
                // LCOV_EXCL_START
                EXCEPTION("Unable to find a face of an element which matches one of the boundary elements");
                // LCOV_EXCL_STOP
            }
        }
    }
}

template<unsigned DIM>
void QuadraticMeshHelper<DIM>::CheckBoundaryElements(AbstractTetrahedralMesh<DIM, DIM>* pMesh)
{
#ifndef NDEBUG
    unsigned expected_num_nodes = DIM*(DIM+1)/2;
    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
             = pMesh->GetBoundaryElementIteratorBegin();
         iter != pMesh->GetBoundaryElementIteratorEnd();
         ++iter)
    {
        assert((*iter)->GetNumNodes() == expected_num_nodes);
    }
#endif
}

template<unsigned DIM>
void QuadraticMeshHelper<DIM>::AddNodeToBoundaryElement(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                                        BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                                        Node<DIM>* pNode)
{
    assert(DIM > 1); // LCOV_EXCL_LINE

    // Add node to the boundary node list
    if (!pNode->IsBoundaryNode())
    {
        pNode->SetAsBoundaryNode();
        pMesh->mBoundaryNodes.push_back(pNode);
    }
    // Add it to the boundary element
    pBoundaryElement->AddNode(pNode);
}

template<unsigned DIM>
void QuadraticMeshHelper<DIM>::AddNodeToBoundaryElement(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                                        BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                                        Element<DIM,DIM>* pElement,
                                                        unsigned internalNode)
{
    assert(DIM > 1); // LCOV_EXCL_LINE
    assert(internalNode >= DIM+1);
    assert(internalNode < (DIM+1)*(DIM+2)/2);
    Node<DIM>* p_internal_node = pElement->GetNode(internalNode);
    AddNodeToBoundaryElement(pMesh, pBoundaryElement, p_internal_node);
}

template<unsigned DIM>
void QuadraticMeshHelper<DIM>::AddExtraBoundaryNodes(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                                     BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                                     Element<DIM,DIM>* pElement,
                                                     unsigned nodeIndexOppositeToFace)
{
    assert(DIM!=1); // LCOV_EXCL_LINE
    if (DIM==2)
    {
        assert(nodeIndexOppositeToFace<3);
        // the single internal node of the element's face will be numbered 'face+3'
        AddNodeToBoundaryElement(pMesh, pBoundaryElement, pElement, nodeIndexOppositeToFace+3);
    }
    else
    {
        assert(DIM==3);

        unsigned b_elem_n0 = pBoundaryElement->GetNodeGlobalIndex(0);
        unsigned b_elem_n1 = pBoundaryElement->GetNodeGlobalIndex(1);

        unsigned offset;
        bool reverse;

        if (nodeIndexOppositeToFace==0)
        {
            // face opposite to node 0 = {1,2,3}, with corresponding internals {9,8,5}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 1, 2, 3, offset, reverse);
            HelperMethod2(pMesh, pBoundaryElement, pElement, 9, 8, 5, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==1)
        {
            // face opposite to node 1 = {2,0,3}, with corresponding internals {7,9,6}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 2, 0, 3, offset, reverse);
            HelperMethod2(pMesh, pBoundaryElement, pElement, 7, 9, 6, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==2)
        {
            // face opposite to node 2 = {0,1,3}, with corresponding internals {8,7,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 3, offset, reverse);
            HelperMethod2(pMesh, pBoundaryElement, pElement, 8, 7, 4, offset, reverse);
        }
        else
        {
            assert(nodeIndexOppositeToFace==3);
            // face opposite to node 3 = {0,1,2}, with corresponding internals {5,6,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 2, offset, reverse);
            HelperMethod2(pMesh, pBoundaryElement, pElement, 5, 6, 4, offset, reverse);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// two unpleasant helper methods for AddExtraBoundaryNodes()
///////////////////////////////////////////////////////////////////////////////

// LCOV_EXCL_START /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void QuadraticMeshHelper<DIM>::HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
                                             Element<DIM,DIM>* pElement,
                                             unsigned node0, unsigned node1, unsigned node2,
                                             unsigned& rOffset,
                                             bool& rReverse)
{
    if (pElement->GetNodeGlobalIndex(node0)==boundaryElemNode0)
    {
        rOffset = 0;
        if (pElement->GetNodeGlobalIndex(node1)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node2)==boundaryElemNode1);
            rReverse = true;
        }
    }
    else if (pElement->GetNodeGlobalIndex(node1)==boundaryElemNode0)
    {
        rOffset = 1;
        if (pElement->GetNodeGlobalIndex(node2)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node0)==boundaryElemNode1);
            rReverse = true;
        }
    }
    else
    {
        assert(pElement->GetNodeGlobalIndex(node2)==boundaryElemNode0);
        rOffset = 2;
        if (pElement->GetNodeGlobalIndex(node0)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node1)==boundaryElemNode1);
            rReverse = true;
        }
    }
}
// LCOV_EXCL_STOP /// \todo These helper methods aren't properly covered


// LCOV_EXCL_START /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void QuadraticMeshHelper<DIM>::HelperMethod2(AbstractTetrahedralMesh<DIM, DIM>* pMesh,
                                             BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                             Element<DIM,DIM>* pElement,
                                             unsigned internalNode0, unsigned internalNode1, unsigned internalNode2,
                                             unsigned offset,
                                             bool reverse)
{
    if (offset==1)
    {
        unsigned temp = internalNode0;
        internalNode0 = internalNode1;
        internalNode1 = internalNode2;
        internalNode2 = temp;
    }
    else if (offset == 2)
    {
        unsigned temp = internalNode0;
        internalNode0 = internalNode2;
        internalNode2 = internalNode1;
        internalNode1 = temp;
    }

    if (reverse)
    {
        unsigned temp = internalNode1;
        internalNode1 = internalNode2;
        internalNode2 = temp;
    }

    AddNodeToBoundaryElement(pMesh, pBoundaryElement, pElement, internalNode0);
    AddNodeToBoundaryElement(pMesh, pBoundaryElement, pElement, internalNode1);
    AddNodeToBoundaryElement(pMesh, pBoundaryElement, pElement, internalNode2);
}
// LCOV_EXCL_STOP /// \todo These helper methods aren't properly covered

// Explicit instantiation
template class QuadraticMeshHelper<1>;
template class QuadraticMeshHelper<2>;
template class QuadraticMeshHelper<3>;
