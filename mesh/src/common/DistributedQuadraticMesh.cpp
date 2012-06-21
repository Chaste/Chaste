/*

Copyright (c) 2005-2012, University of Oxford.
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


#include "DistributedQuadraticMesh.hpp"

template<unsigned DIM>
DistributedQuadraticMesh<DIM>::DistributedQuadraticMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod)
    : DistributedTetrahedralMesh<DIM, DIM>(partitioningMethod)
{

}

template<unsigned DIM>
DistributedQuadraticMesh<DIM>::~DistributedQuadraticMesh()
{
}

template<unsigned DIM>
void DistributedQuadraticMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rAbsMeshReader)
{
    TrianglesMeshReader<DIM, DIM>* p_mesh_reader=dynamic_cast<TrianglesMeshReader<DIM, DIM>*>(&rAbsMeshReader);
    assert(p_mesh_reader != NULL);

    unsigned order_of_elements = 1;
    if (p_mesh_reader)
    {
        //A triangles mesh reader will let you read with non-linear elements
        order_of_elements = p_mesh_reader->GetOrderOfElements();
    }

    assert(order_of_elements == 2);
    // If it is a linear TrianglesMeshReader or any other reader (which are all linear)
//    if (order_of_elements == 1)
//    {
//        WARNING("Reading a (linear) tetrahedral mesh and converting it to a DistributedQuadraticMesh.  This involves making an external library call to Triangle/Tetgen in order to compute internal nodes");
//        ConstructFromLinearMeshReader(rAbsMeshReader);
//        return;
//    }

    DistributedTetrahedralMesh<DIM,DIM>::ConstructFromMeshReader(*p_mesh_reader);
    assert(this->GetNumBoundaryElements() > 0);

    p_mesh_reader->Reset();

    // Add the extra nodes (1 extra node in 1D, 3 in 2D, 6 in 3D) to the element data
    for (unsigned i=0; i<this->GetNumElements(); i++)
    {
        std::vector<unsigned> nodes = p_mesh_reader->GetNextElementData().NodeIndices;
        assert(nodes.size()==(DIM+1)*(DIM+2)/2);
        assert(this->GetElement(i)->GetNumNodes()==DIM+1); // element is initially linear

        // Add extra nodes to make it a quad element
        for (unsigned j=DIM+1; j<(DIM+1)*(DIM+2)/2; j++)
        {
            this->GetElement(i)->AddNode( this->GetNode(nodes[j]) );
            this->GetNode(nodes[j])->AddElement(this->GetElement(i)->GetIndex());
            this->GetNode(nodes[j])->MarkAsInternal();
        }
    }

    //Add the extra nodes to boundary elements
    if (DIM > 1)
    {
        // if  OrderOfBoundaryElements is 2 it can read in the extra nodes for each boundary element, other have to compute them.
        if (p_mesh_reader->GetOrderOfBoundaryElements() == 2u)
        {
            p_mesh_reader->Reset();
            for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
                  = this->GetBoundaryElementIteratorBegin();
                 iter != this->GetBoundaryElementIteratorEnd();
                 ++iter)
            {
                std::vector<unsigned> nodes = p_mesh_reader->GetNextFaceData().NodeIndices;

                assert((*iter)->GetNumNodes()==DIM); // so far just the vertices
                assert(nodes.size()==DIM*(DIM+1)/2); // the reader should have got 6 nodes (3d) for each face

                for (unsigned j=DIM; j<DIM*(DIM+1)/2; j++)
                {
                    Node <DIM> *p_internal_node = this->GetNode(nodes[j]);
                    (*iter)->AddNode( p_internal_node );
                    // add node to the boundary node list
                    if (!p_internal_node->IsBoundaryNode())
                    {
                        p_internal_node->SetAsBoundaryNode();
                        this->mBoundaryNodes.push_back(p_internal_node);
                    }
                }
            }
        }
        else
        {
            AddNodesToBoundaryElements(p_mesh_reader);
        }
    }

    // Check each boundary element has a quadratic number of nodes
#ifndef NDEBUG
    unsigned expected_num_nodes = DIM*(DIM+1)/2;
    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
          = this->GetBoundaryElementIteratorBegin();
          iter != this->GetBoundaryElementIteratorEnd();
          ++iter)
    {
        assert((*iter)->GetNumNodes()==expected_num_nodes);
    }
#endif
}


template<unsigned DIM>
void DistributedQuadraticMesh<DIM>::AddNodesToBoundaryElements(TrianglesMeshReader<DIM,DIM>* pMeshReader)
 {
    // Loop over all boundary elements, find the equivalent face from all
    // the elements, and add the extra nodes to the boundary element
    bool boundary_element_file_has_containing_element_info=false;

    if (pMeshReader)
    {
        boundary_element_file_has_containing_element_info=pMeshReader->GetReadContainingElementOfBoundaryElement();
    }

    if (DIM > 1)
    {
        if (boundary_element_file_has_containing_element_info)
        {
            pMeshReader->Reset();
        }

        for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
               = this->GetBoundaryElementIteratorBegin();
             iter != this->GetBoundaryElementIteratorEnd();
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

                // We know what elem it should be in
                if (boundary_element_file_has_containing_element_info)
                {
                    elem_index = pMeshReader->GetNextFaceData().ContainingElement;
                }

                Element<DIM,DIM>* p_element = this->GetElement(elem_index);

                // For each element, loop over faces (the opposites to a node)
                for (unsigned face=0; face<DIM+1; face++)
                {
                    // Collect the node indices for this face
                    std::set<unsigned> node_indices;
                    for (unsigned local_node_index=0; local_node_index<DIM+1; local_node_index++)
                    {
                        if (local_node_index!=face)
                        {
                            node_indices.insert( p_element->GetNodeGlobalIndex(local_node_index) );
                        }
                    }

                    assert(node_indices.size()==DIM);

                    // see if this face matches the boundary element,
                    // and call AddExtraBoundaryNodes() if so
                    if (node_indices==boundary_element_node_indices)
                    {
                        AddExtraBoundaryNodes(*iter, p_element, face);

                        found_this_boundary_element = true;
                        break;
                    }
                }

                // If the containing element info was given, we must have found the face first time
                if (boundary_element_file_has_containing_element_info && !found_this_boundary_element)
                {
                    #define COVERAGE_IGNORE
                    //std::cout << (*iter)->GetIndex() << " " <<  pMeshReader->GetNextFaceData().ContainingElement << "\n";
                    EXCEPTION("Boundary element " << (*iter)->GetIndex()
                       << "wasn't found in the containing element given for it "
                       << elem_index);
                    #undef COVERAGE_IGNORE
                }

                if (found_this_boundary_element)
                {
                    break;
                }
            }

            if (!found_this_boundary_element)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("Unable to find a face of an element which matches one of the boundary elements");
                #undef COVERAGE_IGNORE
            }
        }
    }
}


template<unsigned DIM>
void DistributedQuadraticMesh<DIM>::AddNodeToBoundaryElement(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                                  Element<DIM,DIM>* pElement,
                                                  unsigned internalNode)
{
    assert(DIM > 1);
    assert(internalNode >= DIM+1);
    assert(internalNode < (DIM+1)*(DIM+2)/2);
    Node<DIM>* p_internal_node = pElement->GetNode(internalNode);

    // Add node to the boundary node list
    if (!p_internal_node->IsBoundaryNode())
    {
        p_internal_node->SetAsBoundaryNode();
        this->mBoundaryNodes.push_back(p_internal_node);
    }

    pBoundaryElement->AddNode( p_internal_node );
}


template<unsigned DIM>
void DistributedQuadraticMesh<DIM>::AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                               Element<DIM,DIM>* pElement,
                                               unsigned nodeIndexOppositeToFace)
{
    assert(DIM!=1);
    if (DIM==2)
    {
        assert(nodeIndexOppositeToFace<3);
        // the single internal node of the elements face will be numbered 'face+3'
        AddNodeToBoundaryElement(pBoundaryElement, pElement, nodeIndexOppositeToFace+3);
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
            HelperMethod2(pBoundaryElement, pElement, 9, 8, 5, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==1)
        {
            // face opposite to node 1 = {2,0,3}, with corresponding internals {7,9,6}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 2, 0, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 7, 9, 6, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==2)
        {
            // face opposite to node 2 = {0,1,3}, with corresponding internals {8,7,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 8, 7, 4, offset, reverse);
        }
        else
        {
            assert(nodeIndexOppositeToFace==3);
            // face opposite to node 3 = {0,1,2}, with corresponding internals {5,6,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 2, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 5, 6, 4, offset, reverse);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// two unpleasant helper methods for AddExtraBoundaryNodes()
///////////////////////////////////////////////////////////////////////////////

#define COVERAGE_IGNORE /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void DistributedQuadraticMesh<DIM>::HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
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
#undef COVERAGE_IGNORE /// \todo These helper methods aren't properly covered


#define COVERAGE_IGNORE /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void DistributedQuadraticMesh<DIM>::HelperMethod2(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
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

    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode0);
    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode1);
    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode2);
}
#undef COVERAGE_IGNORE /// \todo These helper methods aren't properly covered


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class DistributedQuadraticMesh<1>;
template class DistributedQuadraticMesh<2>;
template class DistributedQuadraticMesh<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DistributedQuadraticMesh)
