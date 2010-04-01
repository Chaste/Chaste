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

#include "MutableVertexMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::MutableVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                               std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements,
                                               double cellRearrangementThreshold,
                                               double edgeDivisionThreshold,
                                               double t2Threshold)
    : mCellRearrangementThreshold(cellRearrangementThreshold),
      mCellRearrangementRatio(1.5),
      mEdgeDivisionThreshold(edgeDivisionThreshold),
      mT2Threshold(t2Threshold)
{
    assert(cellRearrangementThreshold > 0.0);
    assert(edgeDivisionThreshold > 0.0);
    assert(t2Threshold > 0.0);

    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements
    for (unsigned node_index=0; node_index<nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index=0; elem_index<vertexElements.size(); elem_index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        this->mElements.push_back(p_temp_vertex_element);
    }

    // Register elements with nodes
    for (unsigned index=0; index<this->mElements.size(); index++)
    {
        VertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp_vertex_element = this->mElements[index];
        for (unsigned node_index=0; node_index<p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = true;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::MutableVertexMesh()
    : mCellRearrangementThreshold(0.01), // Overwritten as soon as archiving is complete
      mCellRearrangementRatio(1.5), // Overwritten as soon as archiving is complete
      mEdgeDivisionThreshold(DBL_MAX), // Overwritten as soon as archiving is complete
      mT2Threshold(0.001) // Overwritten as soon as archiving is complete
{
    this->mMeshChangesDuringSimulation = true;
    Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::~MutableVertexMesh()
{
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementThreshold() const
{
    return mCellRearrangementThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetEdgeDivisionThreshold() const
{
    return mEdgeDivisionThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetT2Threshold() const
{
    return mT2Threshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementThreshold(double cellRearrangementThreshold)
{
    mCellRearrangementThreshold = cellRearrangementThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetEdgeDivisionThreshold(double edgeDivisionThreshold)
{
    mEdgeDivisionThreshold = edgeDivisionThreshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetT2Threshold(double t2Threshold)
{
    mT2Threshold = t2Threshold;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    mDeletedNodeIndices.clear();
    mDeletedElementIndices.clear();

    VertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - mDeletedNodeIndices.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - mDeletedElementIndices.size();
}

//\todo deal with boundary elements in vertex meshes.
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
//{
//    return mBoundaryElements.size();
//}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(this->mNodes.size());
        this->mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete this->mNodes[index];
        this->mNodes[index] = pNewNode;
    }
    return pNewNode->GetIndex();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pNewElement)
{
    unsigned new_element_index = pNewElement->GetIndex();

    if (new_element_index == this->mElements.size())
    {
        this->mElements.push_back(pNewElement);
    }
    else
    {
        this->mElements[new_element_index] = pNewElement;
    }
    pNewElement->RegisterWithNodes();
    return pNewElement->GetIndex();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    this->mNodes[nodeIndex]->SetPoint(point);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DeleteElementPriorToReMesh(unsigned index)
{
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
    #undef COVERAGE_IGNORE

    // Mark any nodes that are contained only in this element as deleted
    for (unsigned i=0; i<this->mElements[index]->GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* p_node = this->mElements[index]->GetNode(i);

        if (p_node->rGetContainingElementIndices().size()==1)
        {
            p_node->MarkAsDeleted();
            mDeletedNodeIndices.push_back(p_node->GetIndex());
        }

        // Mark all the nodes contained in the removed element as boundary nodes
        p_node->SetAsBoundaryNode(true);

    }

    // Mark this element as deleted
    this->mElements[index]->MarkAsDeleted();
    mDeletedElementIndices.push_back(index);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(shared_elements.size() > 0);

    // Create a new node (position is not important as it will be changed)
    Node<SPACE_DIM>* p_new_node = new Node<SPACE_DIM>(GetNumNodes(), false, 0.0, 0.0);

    // Update the node location
    c_vector<double, SPACE_DIM> new_node_position = pNodeA->rGetLocation() + 0.5*GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation());
    ChastePoint<SPACE_DIM> point(new_node_position);
    p_new_node->SetPoint(new_node_position);

    // Add node to mesh
    this->mNodes.push_back(p_new_node);

    // Iterate over common elements
    for (std::set<unsigned>::iterator iter = shared_elements.begin();
         iter != shared_elements.end();
         ++iter)
    {
        // Find which node has the lower local index in this element
        unsigned local_indexA = this->GetElement(*iter)->GetNodeLocalIndex(pNodeA->GetIndex());
        unsigned local_indexB = this->GetElement(*iter)->GetNodeLocalIndex(pNodeB->GetIndex());

        unsigned index = local_indexB;
        if ( local_indexB > local_indexA )
        {
            index = local_indexA;
        }
        if ( (local_indexA == 0) && (local_indexB == this->GetElement(*iter)->GetNumNodes()-1))
        {
            index = local_indexB;
        }
        if ( (local_indexB == 0) && (local_indexA == this->GetElement(*iter)->GetNumNodes()-1))
        {
            index = local_indexA;
        }
        // Add new node to this element
        this->GetElement(*iter)->AddNode(index, p_new_node);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodesAndElements(VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2 || SPACE_DIM==3);
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Make sure the map is big enough
    rElementMap.Resize(this->GetNumAllElements());

    // Remove deleted elements
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> live_elements;
    for (unsigned i=0; i< this->mElements.size(); i++)
    {
        if (this->mElements[i]->IsDeleted())
        {
            delete this->mElements[i];
        }
        else
        {
            live_elements.push_back(this->mElements[i]);
            rElementMap.SetNewIndex(i, (unsigned)(live_elements.size()-1));
        }
    }

    assert(mDeletedElementIndices.size() == this->mElements.size() - live_elements.size());
    mDeletedElementIndices.clear();
    this->mElements = live_elements;

    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        this->mElements[i]->ResetIndex(i);
    }

    // Remove deleted nodes
    RemoveDeletedNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodes()
{
    std::vector<Node<SPACE_DIM>*> live_nodes;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            delete this->mNodes[i];
        }
        else
        {
            live_nodes.push_back(this->mNodes[i]);
        }
    }

    assert(mDeletedNodeIndices.size() == this->mNodes.size() - live_nodes.size());
    this->mNodes = live_nodes;
    mDeletedNodeIndices.clear();

    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->SetIndex(i);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM==2 || SPACE_DIM==3);
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    if (SPACE_DIM==2)
    {
        // Remove deleted nodes and elements
        RemoveDeletedNodesAndElements(rElementMap);

        /*
         * We do not need to call Clear() and remove all current data, since
         * cell birth, rearrangement and death result only in local remeshing
         * of a vertex-based mesh.
         */

        // Restart check after each T1/T2Swap as elements are changed
        bool recheck_mesh = true;
        while (recheck_mesh == true)
        {
            recheck_mesh = false;

            // Loop over elements to check for T2Swaps
            // Separate loops as need to check for T2Swaps first
           for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = this->GetElementIteratorBegin();
                 elem_iter != this->GetElementIteratorEnd();
                 ++elem_iter)
            {
                if (!recheck_mesh)
                {
                    if (elem_iter->GetNumNodes() == 3u)
                    {
                        /*
                         * Perform T2 swaps where necesary
                         * Check there are only 3 nodes and the element is small enough
                         */
                        if (GetAreaOfElement(elem_iter->GetIndex()) < GetT2Threshold())
                        {
                            PerformT2Swap(*elem_iter);
                            // Now remove the deleted nodes (if we don't do this then the search for T1Swap causes errors)
                            RemoveDeletedNodesAndElements(rElementMap);
                            recheck_mesh = true;
                            break;
                        }
                    }
                }
            }
            // Loop over elements to check for T1Swaps
            for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = this->GetElementIteratorBegin();
                 elem_iter != this->GetElementIteratorEnd();
                 ++elem_iter)
            {
                if (!recheck_mesh)
                {
                    unsigned num_nodes = elem_iter->GetNumNodes();
                    assert(num_nodes > 0); // if not element should be deleted

                    unsigned new_num_nodes = num_nodes;

                    /*
                     *  Perform T1 swaps and divide edges where necesary
                     *  Check there are > 3 nodes in both elements that contain the pair of nodes
                     *  and the edges are small enough
                     */

                    // Loop over element vertices
                    for (unsigned local_index=0; local_index<num_nodes; local_index++)
                    {
                        // Find locations of current node and anticlockwise node
                        Node<SPACE_DIM>* p_current_node = elem_iter->GetNode(local_index);
                        unsigned local_index_plus_one = (local_index+1)%new_num_nodes; /// \todo Should use iterators to tidy this up
                        Node<SPACE_DIM>* p_anticlockwise_node = elem_iter->GetNode(local_index_plus_one);

                        // Find distance between nodes
                        double distance_between_nodes = this->GetDistanceBetweenNodes(p_current_node->GetIndex(), p_anticlockwise_node->GetIndex());

                        std::set<unsigned> elements_of_node_a = p_current_node->rGetContainingElementIndices();
                        std::set<unsigned> elements_of_node_b = p_anticlockwise_node->rGetContainingElementIndices();

                        std::set<unsigned> all_elements;
                        std::set_union(elements_of_node_a.begin(), elements_of_node_a.end(),
                                       elements_of_node_b.begin(), elements_of_node_b.end(),
                                       std::inserter(all_elements, all_elements.begin()));

                        // Track if either node is in a trinagular element.
                        bool triangular_element = false;

                        for (std::set<unsigned>::const_iterator it = all_elements.begin();
                             it != all_elements.end();
                             ++it)
                        {
                            ///\todo this magic number needs to be investigated.
                            if (this->GetElement(*it)->GetNumNodes() <= 3)//&&(this->GetAreaOfElement(*it) < 2.0*GetT2Threshold()))
                            {
                                triangular_element = true;
                            }
                        }

                        // If the nodes are too close together and we don't have any triangular elements connected to the nodes, perform a swap
                        if ((!triangular_element) && (distance_between_nodes < mCellRearrangementThreshold))
                        {
                            // Identify the type of node swap/merge needed then call method to perform swap/merge
                            IdentifySwapType(p_current_node, p_anticlockwise_node, rElementMap);
                            recheck_mesh = true;
                            break;
                        }

                        if (distance_between_nodes > mEdgeDivisionThreshold)
                        {
                            // If the nodes are too far apart, divide the edge
                            DivideEdge(p_current_node, p_anticlockwise_node);
                            new_num_nodes++;
                        }
                    }
                }
                else
                {
                    break;
                }
            }
        }

        // ... end of element rearrangement code

        // areas and perimeters of elements are sorted in PerformT1Swap() method


        // Restart check after each T3Swap as elements are changed
        // \todo dont need to do this at present as # of elements crrently stays the same but will need this when tracking voids
        recheck_mesh = true;
        while (recheck_mesh == true)
        {
            recheck_mesh = false;

            // Check that no nodes have overlapped elements
            /// \todo Only need to check this next bit if the element/node is on the boundary (see #933 and #943)
            for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = this->GetElementIteratorBegin();
                 elem_iter != this->GetElementIteratorEnd();
                 ++elem_iter)
            {
                unsigned num_nodes = elem_iter->GetNumNodes();

                if (!recheck_mesh)
                {
                    // Loop over element vertices
                    for (unsigned local_index=0; local_index<num_nodes; local_index++)
                    {
                        // Find locations of current node and anticlockwise node
                        Node<SPACE_DIM>* p_current_node = elem_iter->GetNode(local_index);

                        if (p_current_node->IsBoundaryNode())
                        {
                            for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator other_iter = this->GetElementIteratorBegin();
                                 other_iter != this->GetElementIteratorEnd();
                                 ++other_iter)
                            {
                                if (other_iter != elem_iter)
                                {
                                    if (ElementIncludesPoint(p_current_node->rGetLocation(), other_iter->GetIndex()))
                                    {
                                        PerformT3Swap(p_current_node, other_iter->GetIndex());
                                        recheck_mesh = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else // 3D
    {
        #define COVERAGE_IGNORE
        EXCEPTION("Remeshing has not been implemented in 3D (see #827 and #860)\n");
        #undef COVERAGE_IGNORE
        /// \todo put code for remeshing in 3D here - see #866 and the paper doi:10.1016/j.jtbi.2003.10.001
    }

}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    VertexElementMap map(GetNumElements());
    ReMesh(map);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // this method only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();

    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    // Form the set union
    std::set<unsigned> all_indices, temp_set;
    std::set_union(nodeA_elem_indices.begin(), nodeA_elem_indices.end(),
                   nodeB_elem_indices.begin(), nodeB_elem_indices.end(),
                   std::inserter(temp_set, temp_set.begin()));
    all_indices.swap(temp_set); // temp_set will be deleted

    if ((nodeA_elem_indices.size()>3) || (nodeB_elem_indices.size()>3))
    {
        /*
         * Looks like
         *
         *  \
         *   \ A   B
         * ---o---o---
         *   /
         *  /
         *
         */
        EXCEPTION("A node is contained in more than three elements"); // the code can't handle this case
    }
    else // each node is contained in at most three elements
    {
        switch (all_indices.size())
        {
            case 1:
            {
                /*
                 * In this case, each node is contained in a single element, so the nodes
                 * lie on the boundary of the mesh:
                 *
                 *    A   B
                 * ---o---o---
                 *
                 * We merge the nodes.
                 */
                PerformNodeMerge(pNodeA, pNodeB);

                // Remove the deleted node and re-index
                RemoveDeletedNodes();
                break;
            }
            case 2:
            {
                if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                {
                    if (pNodeA->IsBoundaryNode() && pNodeB->IsBoundaryNode())
                    {
                        /*
                         * In this case, the node configuration looks like:
                         *
                         *   \   /
                         *    \ / Node A
                         * (1) |   (2)      (element number in brackets)
                         *    / \ Node B
                         *   /   \
                         *
                         * We perform a Type 1 swap and seperate the elements in this case.
                         */
                         PerformT1Swap(pNodeA, pNodeB,all_indices);
                    }
                    else
                    {
                        /*
                         * In this case, each node is contained in two elements, so the nodes
                         * lie on an internal edge:
                         *
                         *    A   B
                         * ---o---o---
                         *
                         * We merge the nodes in this case.
                         */
                        PerformNodeMerge(pNodeA, pNodeB);

                        // Remove the deleted node and re-index
                        RemoveDeletedNodes();
                    }
                }
                else
                {
                    /*
                     * In this case, the node configuration looks like:
                     *
                     * Outside
                     *         /
                     *   --o--o (2)
                     *     (1) \
                     *
                     * We merge the nodes in this case.
                     */
                     PerformNodeMerge(pNodeA, pNodeB);

                     // Remove the deleted node and re-index
                     RemoveDeletedNodes();
                }
                break;
            }
            case 3:
            {
                if (nodeA_elem_indices.size()==1 || nodeB_elem_indices.size()==1)
                {
                    /*
                     * In this case, one node is contained in one element and
                     * the other node is contained in three elements:
                     *
                     *    A   B
                     *
                     *  empty   /
                     *         / (3)
                     * ---o---o-----   (element number in brackets)
                     *  (1)    \ (2)
                     *          \
                     *
                     * We perform a node merge in this case.
                     */
                    PerformNodeMerge(pNodeA, pNodeB);

                    // Remove the deleted node and re-index
                    RemoveDeletedNodes();
                }
                //\todo use the fact that both nodes are boundary nodes here
                else if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                {
                    // element in nodeA_elem_indices which is not in nodeB_elem_indices contains a shared node with the element in nodeA_elem_indices which is not in nodeB_elem_indices.

                    std::set<unsigned> element_A_not_B, temp_set;
                    std::set_difference(all_indices.begin(), all_indices.end(),
                                          nodeB_elem_indices.begin(), nodeB_elem_indices.end(),
                                       std::inserter(temp_set, temp_set.begin()));
                    element_A_not_B.swap(temp_set);

                    //Should only be one such element
                    assert(element_A_not_B.size()==1);

                    std::set<unsigned> element_B_not_A;
                    std::set_difference(all_indices.begin(), all_indices.end(),
                                          nodeA_elem_indices.begin(), nodeA_elem_indices.end(),
                                       std::inserter(temp_set, temp_set.begin()));
                    element_B_not_A.swap(temp_set);

                    //Should only be one such element
                    assert(element_B_not_A.size()==1);


                    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_A_not_B = this->mElements[*element_A_not_B.begin()];
                    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_B_not_A = this->mElements[*element_B_not_A.begin()];

                    unsigned local_index_1 = p_element_A_not_B->GetNodeLocalIndex(pNodeA->GetIndex());
                    unsigned next_node_1 = p_element_A_not_B->GetNodeGlobalIndex((local_index_1 + 1)%(p_element_A_not_B->GetNumNodes()));
                    unsigned previous_node_1 = p_element_A_not_B->GetNodeGlobalIndex((local_index_1 + p_element_A_not_B->GetNumNodes() - 1)%(p_element_A_not_B->GetNumNodes()));
                    unsigned local_index_2 = p_element_B_not_A->GetNodeLocalIndex(pNodeB->GetIndex());
                    unsigned next_node_2 = p_element_B_not_A->GetNodeGlobalIndex((local_index_2 + 1)%(p_element_B_not_A->GetNumNodes()));
                    unsigned previous_node_2 = p_element_B_not_A->GetNodeGlobalIndex((local_index_2 + p_element_B_not_A->GetNumNodes() - 1)%(p_element_B_not_A->GetNumNodes()));

                    if (next_node_1 == previous_node_2 || next_node_2 == previous_node_1)
                     {
                        /*
                         * In this case, the node configuration looks like:
                         *
                         *     A  B                  A  B
                         *      /\                 \      /
                         *     /v \                 \(1) /
                         * (3) o--o (1)  or      (2) o--o (3)    (element number in brackets, v is a void)
                         *    / (2)\                 \v /
                         *   /      \                 \/
                         *
                         * We perform a T3 way merge, removing the void.
                         */
                        unsigned nodeC_index;
                        if (next_node_1 == previous_node_2 && next_node_2 != previous_node_1)
                        {
                            nodeC_index = next_node_1;
                        }
                        else if (next_node_2 == previous_node_1 && next_node_1 != previous_node_2)
                        {
                            nodeC_index = next_node_2;
                        }
                        else
                        {
                             assert(next_node_1 == previous_node_2 && next_node_2 == previous_node_1);
                             EXCEPTION("Triangular element next to triangular void, not implemented yet.");
                        }

                        Node<SPACE_DIM>* p_nodeC = this->mNodes[nodeC_index];

                        // \todo this should be a helper method "PerformVoidRemoval(pNodeA, pNodeB, p_nodeC)";

                        unsigned nodeA_index = pNodeA->GetIndex();
                        unsigned nodeB_index = pNodeB->GetIndex();

                        c_vector<double, SPACE_DIM> nodes_midpoint = pNodeA->rGetLocation() + this->GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation())/3.0
                                                                                            + this->GetVectorFromAtoB(pNodeA->rGetLocation(), p_nodeC->rGetLocation())/3.0;

                        Node<SPACE_DIM>*  p_low_node_A_B = (nodeA_index < nodeB_index) ? pNodeA : pNodeB; // Node with the lowest index out of A and B
                        Node<SPACE_DIM>*  p_low_node = (p_low_node_A_B->GetIndex() < nodeC_index) ? p_low_node_A_B : p_nodeC; // Node with the lowest index out of A, B and C.
                        PerformNodeMerge(pNodeA,pNodeB);
                        PerformNodeMerge(p_low_node_A_B,p_nodeC);

                        c_vector<double, SPACE_DIM>& r_low_node_location = p_low_node->rGetModifiableLocation();

                         r_low_node_location = nodes_midpoint;

                           // Sort out boundary nodes
                        p_low_node->SetAsBoundaryNode(false);

                        // Remove the deleted nodes and re-index
                        RemoveDeletedNodes();
                    }
                    else
                    {
                        /*
                         * In this case, the node configuration looks like:
                         *
                         *     A  B                  A  B
                         *   \ empty/              \      /
                         *    \    /                \(1) /
                         * (3) o--o (1)  or      (2) o--o (3)    (element number in brackets)
                         *    / (2)\                /    \
                         *   /      \              /empty \
                         *
                         * We perform a Type 1 swap in this case.
                         */
                        PerformT1Swap(pNodeA, pNodeB, all_indices);
                    }
                }
                else
                {
                    /*
                     * In this case, one of the nodes is contained in two elements
                     * and the other node is contained in three elements.
                     */
                    assert (   (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==3)
                            || (nodeA_elem_indices.size()==3 && nodeB_elem_indices.size()==2) );

                    /*
                     * Let node alpha be the node contained in two elements and node beta
                     * the node contained in three elements.
                     */
                    Node<SPACE_DIM>* p_node_alpha = (nodeA_elem_indices.size()==2) ? pNodeA : pNodeB;
                    Node<SPACE_DIM>* p_node_beta = (nodeA_elem_indices.size()==2) ? pNodeB : pNodeA;

                    // Get the set of elements containing node alpha and assert there are two such elements
                    std::set<unsigned> node_alpha_elem_indices = p_node_alpha->rGetContainingElementIndices();
                    assert(node_alpha_elem_indices.size() == 2u);

                    // Get a pointer to the first of these elements and assert it contains node alpha
                    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->mElements[*node_alpha_elem_indices.begin()];

                    unsigned node_alpha_local_index = p_element->GetNodeLocalIndex(p_node_alpha->GetIndex());
                    assert(node_alpha_local_index < UINT_MAX);

                    /*
                     * Get the local indices of the nodes neighbouring node alpha in this element.
                     * We add an extra temp in the second call to the % operator term to ensure we
                     * are operating on a positive number (otherwise the % operator may break).
                     */
                    unsigned temp = p_element->GetNumNodes();
                    unsigned node_alpha_local_index_before = (node_alpha_local_index+1)%temp;
                    unsigned node_alpha_local_index_after = (node_alpha_local_index+temp-1)%temp;

                    // Get pointers to these nodes and assert one of them is p_node_beta
                    Node<SPACE_DIM>* p_node1 = p_element->GetNode(node_alpha_local_index_before);
                    Node<SPACE_DIM>* p_node2 = p_element->GetNode(node_alpha_local_index_after);
                    assert(p_node1 == p_node_beta || p_node2 == p_node_beta);

                    // Get whichever of the nodes is NOT p_node_beta, call it node gamma
                    Node<SPACE_DIM>* p_node_gamma = (p_node1 == p_node_beta) ? p_node2 : p_node1;

                    // Get the set of elements containing nodes beta and gamma
                    std::set<unsigned> node_beta_elem_indices = p_node_beta->rGetContainingElementIndices();
                    std::set<unsigned> node_gamma_elem_indices = p_node_gamma->rGetContainingElementIndices();

                    // Form the set intersection
                    std::set<unsigned> intersection_indices, temp_set2;
                    std::set_intersection(node_beta_elem_indices.begin(), node_beta_elem_indices.end(),
                                   node_gamma_elem_indices.begin(), node_gamma_elem_indices.end(),
                                   std::inserter(temp_set2, temp_set2.begin()));
                    intersection_indices.swap(temp_set2); // temp_set2 will be deleted

                    // The number of intersecting indices must be 1 or 2
                    switch (intersection_indices.size())
                    {
                        case 1:
                        {
                            /*
                             * In this case, the node configuration looks like:
                             *
                             *     A  B                      A  B
                             *   \      /                  \      /
                             *    \ (1)/                    \(1) /
                             * (3) o--o (empty)  or  (empty) o--o (3)    (element number in brackets)
                             *    / (2)\                    /(2) \
                             *   /      \                  /      \
                             *
                             * We perform a Type 1 swap in this case.
                             */
                            PerformT1Swap(pNodeA, pNodeB, all_indices);
                            break;
                        }
                        case 2:
                        {
                            /*
                             * In this case, the node configuration looks like:
                             *
                             *     A  B             A  B
                             *   \                       /
                             *    \  (1)           (1)  /
                             * (3) o--o---   or  ---o--o (3)    (element number in brackets)
                             *    /  (2)           (2)  \
                             *   /                       \
                             *
                             * We perform a node merge in this case.
                             */
                            PerformNodeMerge(pNodeA, pNodeB);

                            // Remove the deleted node and re-index
                            RemoveDeletedNodes();
                            break;
                        }
                        default:
                            NEVER_REACHED;
                    }
                }
                break;
            }
            case 4:
            {
                /*
                 * In this case, the node configuration looks like:
                 *
                 *   \(1)/
                 *    \ / Node A
                 * (2) |   (4)      (element number in brackets)
                 *    / \ Node B
                 *   /(3)\
                 *
                 * We perform a Type 1 swap in this case.
                 */
                PerformT1Swap(pNodeA, pNodeB, all_indices);
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformNodeMerge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    // Sort the nodes by index
    unsigned nodeA_index = pNodeA->GetIndex();
    unsigned nodeB_index = pNodeB->GetIndex();

    unsigned lo_node_index = (nodeA_index < nodeB_index) ? nodeA_index : nodeB_index; // low index
    unsigned hi_node_index = (nodeA_index < nodeB_index) ? nodeB_index : nodeA_index; // high index

    // Get pointers to the nodes, sorted by index
    Node<SPACE_DIM>* p_lo_node = this->GetNode(lo_node_index);
    Node<SPACE_DIM>* p_hi_node = this->GetNode(hi_node_index);

    // Find the sets of elements containing each of the nodes, sorted by index
    std::set<unsigned> lo_node_elem_indices = p_lo_node->rGetContainingElementIndices();
    std::set<unsigned> hi_node_elem_indices = p_hi_node->rGetContainingElementIndices();

    // Move the low-index node to the mid-point
    c_vector<double, SPACE_DIM> node_midpoint = p_lo_node->rGetLocation() + 0.5*this->GetVectorFromAtoB(p_lo_node->rGetLocation(), p_hi_node->rGetLocation());
    c_vector<double, SPACE_DIM>& r_lo_node_location = p_lo_node->rGetModifiableLocation();
    r_lo_node_location = node_midpoint;

    // Update the elements previously containing the high-index node to contain the low-index node
    for (std::set<unsigned>::const_iterator it = hi_node_elem_indices.begin();
         it != hi_node_elem_indices.end();
         ++it)
    {
        // Find the local index of the high-index node in this element
        unsigned hi_node_local_index = this->mElements[*it]->GetNodeLocalIndex(hi_node_index);
        assert(hi_node_local_index < UINT_MAX); // this element should contain the high-index node

        /*
         * If this element already contains the low-index node, then just remove the high-index node.
         * Otherwise replace it with the low-index node in the element and remove it from mNodes.
         */
        if (lo_node_elem_indices.count(*it) > 0)
        {
            this->mElements[*it]->DeleteNode(hi_node_local_index); // think this method removes the high-index node from mNodes
        }
        else
        {
            // Replace the high-index node with the low-index node in this element
            this->mElements[*it]->UpdateNode(hi_node_local_index, p_lo_node);
        }
    }
    assert(!(this->mNodes[hi_node_index]->IsDeleted()));
    this->mNodes[hi_node_index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(hi_node_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT1Swap(Node<SPACE_DIM>* pNodeA,
                                                       Node<SPACE_DIM>* pNodeB,
                                                       std::set<unsigned>& rElementsContainingNodes)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    /*
     * Restructure elements - remember to update nodes and elements.
     *
     * We need to implement the following changes:
     *
     * The element whose index was in nodeA_elem_indices but not nodeB_elem_indices,
     * and the element whose index was in nodeB_elem_indices but not nodeA_elem_indices,
     * should now both contain nodes A and B.
     *
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which
     * node C lies inside, should now only contain node A.
     *
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which
     * node D lies inside, should now only contain node B.
     *
     * Iterate over all elements involved and identify which element they are
     * in the diagram then update the nodes as necessary.
     *
     *   \(1)/
     *    \ / Node A
     * (2) |   (4)     elements in brackets
     *    / \ Node B
     *   /(3)\
     *
     */

    /*
     * Compute the locations of two new nodes C, D, placed on either side of the
     * edge E_old formed by nodes current_node and anticlockwise_node, such
     * that the edge E_new formed by the new nodes is the perpendicular bisector
     * of E_old, with |E_new| 'just larger' (mCellRearrangementRatio) than mThresholdDistance.
     */

    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    double distance_between_nodes_CD = mCellRearrangementRatio*mCellRearrangementThreshold;

    c_vector<double, SPACE_DIM> nodeA_location = pNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> nodeB_location = pNodeB->rGetLocation();

    c_vector<double, SPACE_DIM> a_to_b = this->GetVectorFromAtoB(nodeA_location, nodeB_location);
    c_vector<double, SPACE_DIM> perpendicular_vector;
    perpendicular_vector(0) = -a_to_b(1);
    perpendicular_vector(1) = a_to_b(0);

    c_vector<double, SPACE_DIM> c_to_d = distance_between_nodes_CD / norm_2(a_to_b) * perpendicular_vector;
    c_vector<double, SPACE_DIM> nodeC_location = nodeA_location + 0.5*a_to_b - 0.5*c_to_d;
    c_vector<double, SPACE_DIM> nodeD_location = nodeC_location + c_to_d;

    /*
     * Move node A to C and node B to D
     */

    c_vector<double, SPACE_DIM>& r_nodeA_location = pNodeA->rGetModifiableLocation();
    r_nodeA_location = nodeC_location;

    c_vector<double, SPACE_DIM>& r_nodeB_location = pNodeB->rGetModifiableLocation();
    r_nodeB_location = nodeD_location;

    for (std::set<unsigned>::const_iterator it = rElementsContainingNodes.begin();
         it != rElementsContainingNodes.end();
         ++it)
    {
        if (nodeA_elem_indices.find(*it) == nodeA_elem_indices.end()) // not in nodeA_elem_indices so element 3
        {
            /*
             * In this case the element index was not in
             * nodeA_elem_indices, so this element
             * does not contain node A. Therefore we must add node A
             * (which has been moved to node C) to this element.
             *
             * Locate local index of node B in element then add node A after
             * in anticlockwise direction.
             */

            unsigned nodeB_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX); // this element should contain node B

            this->mElements[*it]->AddNode(nodeB_local_index, pNodeA);
        }
        else if (nodeB_elem_indices.find(*it) == nodeB_elem_indices.end()) // not in nodeB_elem_indices so element 1
        {
            /*
             * In this case the element index was not in
             * nodeB_elem_indices, so this element
             * does not contain node B. Therefore we must add node B
             * (which has been moved to node D) to this element.
             *
             * Locate local index of node A in element then add node B after
             * in anticlockwise direction.
             */
            unsigned nodeA_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX); // this element should contain node A
            this->mElements[*it]->AddNode(nodeA_local_index, pNodeB);
        }
        else
        {
            /*
             * In this case the element index was in both nodeB_elem_indices and nodeB_elem_indices
             * so is element 2 or 4
             */

            /*
             * Locate local index of nodeA and nodeB and use the oredering to
             * identify the element, if nodeB_index > nodeA_index then element 4
             * and if nodeA_index > nodeB_index then element 2
             */
            unsigned nodeA_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX); // this element should contain node A

            unsigned nodeB_local_index = this->mElements[*it]->GetNodeLocalIndex(pNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX); // this element should contain node B

            unsigned nodeB_local_index_plus_one = (nodeB_local_index + 1)%(this->mElements[*it]->GetNumNodes());

            if (nodeA_local_index == nodeB_local_index_plus_one)
            {
                /*
                 * In this case the local index of nodeA is the local index of
                 * nodeB plus one so we are in element 2 so we remove nodeB
                 */
                this->mElements[*it]->DeleteNode(nodeB_local_index);
            }
            else
            {
                assert(nodeB_local_index == (nodeA_local_index + 1)%(this->mElements[*it]->GetNumNodes())); // as A and B are next to each other
                /*
                 * In this case the local index of nodeA is the local index of
                 * nodeB minus one so we are in element 4 so we remove nodeA
                 */
                this->mElements[*it]->DeleteNode(nodeA_local_index);
            }
        }
    }

    // Sort out boundary nodes
    if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
    {
        if (pNodeA->GetNumContainingElements()==3)
        {
            pNodeA->SetAsBoundaryNode(false);
        }
        else
        {
            pNodeA->SetAsBoundaryNode(true);
        }
        if (pNodeB->GetNumContainingElements()==3)
        {
            pNodeB->SetAsBoundaryNode(false);
        }
        else
        {
            pNodeB->SetAsBoundaryNode(true);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT2Swap(VertexElement<ELEMENT_DIM,SPACE_DIM>& rElement)
{
   // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    /* For removing a small triangle
     *
     *  \
     *   \
     *   |\
     *   | \
     *   |  \_ _ _ _ _
     *   |  /
     *   | /
     *   |/
     *   /
     *  /
     *
     * To

     *  \
     *   \
     *    \
     *     \
     *      \_ _ _ _ _
     *      /
     *     /
     *    /
     *   /
     *  /
     *
     */

     // Assert that the triangle element has only three nodes (!)

     assert(rElement.GetNumNodes() == 3u);

     c_vector<double, SPACE_DIM>& new_node_location = rElement.GetNode(0)->rGetModifiableLocation();
     new_node_location = GetCentroidOfElement(rElement.GetIndex());

     c_vector<unsigned, 3> neighbouring_elem_nums;

     for (unsigned i=0; i<3; i++)
     {
         std::set<unsigned> elements_of_node_a = rElement.GetNode((i+1)%3)->rGetContainingElementIndices();
         std::set<unsigned> elements_of_node_b = rElement.GetNode((i+2)%3)->rGetContainingElementIndices();

         std::set<unsigned> common_elements;
         std::set_intersection(elements_of_node_a.begin(), elements_of_node_a.end(),
                                elements_of_node_b.begin(), elements_of_node_b.end(),
                                std::inserter(common_elements, common_elements.begin()));

         assert(common_elements.size() == 2u);
         common_elements.erase(rElement.GetIndex());
         assert(common_elements.size() == 1u);

         neighbouring_elem_nums(i) = *(common_elements.begin());
     }

     // Extract the neighbouring elements
     VertexElement<ELEMENT_DIM,SPACE_DIM>* p_neighbouring_element_0 = this->GetElement(neighbouring_elem_nums(0));
     VertexElement<ELEMENT_DIM,SPACE_DIM>* p_neighbouring_element_1 = this->GetElement(neighbouring_elem_nums(1));
     VertexElement<ELEMENT_DIM,SPACE_DIM>* p_neighbouring_element_2 = this->GetElement(neighbouring_elem_nums(2));

     // Need to check that none of the neighbouring elements are triangles
     if (    (p_neighbouring_element_0->GetNumNodes() > 3u)
         && (p_neighbouring_element_1->GetNumNodes() > 3u)
         && (p_neighbouring_element_2->GetNumNodes() > 3u) )
     {
         // Neighbour 0 - replace node 1 with node 0, delete node 2
         p_neighbouring_element_0->ReplaceNode(rElement.GetNode(1), rElement.GetNode(0));
         p_neighbouring_element_0->DeleteNode(p_neighbouring_element_0->GetNodeLocalIndex(rElement.GetNodeGlobalIndex(2)));

         // Neighbour 1 - delete node 2
         p_neighbouring_element_1->DeleteNode(p_neighbouring_element_1->GetNodeLocalIndex(rElement.GetNodeGlobalIndex(2)));

         // Neighbour 2 - delete node 1
         p_neighbouring_element_2->DeleteNode(p_neighbouring_element_2->GetNodeLocalIndex(rElement.GetNodeGlobalIndex(1)));

         // Also have to mark pElement, pElement->GetNode(1), pElement->GetNode(2) as deleted.
         mDeletedNodeIndices.push_back(rElement.GetNodeGlobalIndex(1));
         mDeletedNodeIndices.push_back(rElement.GetNodeGlobalIndex(2));
         rElement.GetNode(1)->MarkAsDeleted();
         rElement.GetNode(2)->MarkAsDeleted();

         mDeletedElementIndices.push_back(rElement.GetIndex());
         rElement.MarkAsDeleted();
     }
     else
     {
        EXCEPTION("One of the neighbours of a small triangular element is also a triangle - dealing with this has not been implemented yet");
     }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongGivenAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, c_vector<double, SPACE_DIM> axisOfDivision,
                                                                                bool placeOriginalElementBelow)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Get the centroid of the element
    c_vector<double, SPACE_DIM> centroid = GetCentroidOfElement(pElement->GetIndex());

    // Create a vector perpendicular to the axis of division
    c_vector<double, SPACE_DIM> perp_axis;
    perp_axis(0) = -axisOfDivision(1);
    perp_axis(1) = axisOfDivision(0);

    /*
     * Find which edges the axis of division crosses by finding any node
     * that lies on the opposite side of the axis of division to its next
     * neighbour.
     */
    unsigned num_nodes = pElement->GetNumNodes();
    std::vector<unsigned> intersecting_nodes;
    for (unsigned i=0; i<num_nodes; i++)
    {
        bool is_current_node_on_left = (inner_prod(GetVectorFromAtoB(pElement->GetNodeLocation(i), centroid), perp_axis) >= 0);
        bool is_next_node_on_left = (inner_prod(GetVectorFromAtoB(pElement->GetNodeLocation((i+1)%num_nodes), centroid), perp_axis) >= 0);

        if (is_current_node_on_left != is_next_node_on_left)
        {
            intersecting_nodes.push_back(i);
        }
    }

    // If the axis of division does not cross two edges then we cannot proceed
    if (intersecting_nodes.size() != 2)
    {
        EXCEPTION("Cannot proceed with element division: the given axis of division does not cross two edges of the element");
    }

    std::vector<unsigned> division_node_global_indices;
    unsigned nodes_added = 0;

    // Find the intersections between the axis of division and the element edges
    for (unsigned i=0; i<intersecting_nodes.size(); i++)
    {
        /*
         * Get pointers to the nodes forming the edge into which one new node will be inserted.
         *
         * Note that when we use the first entry of intersecting_nodes to add a node,
         * we change the local index of the second entry of intersecting_nodes in
         * pElement, so must account for this by moving one entry further on.
         */
        Node<SPACE_DIM>* p_node_A = pElement->GetNode((intersecting_nodes[i]+nodes_added)%pElement->GetNumNodes());
        Node<SPACE_DIM>* p_node_B = pElement->GetNode((intersecting_nodes[i]+nodes_added+1)%pElement->GetNumNodes());

        // Find the indices of the elements owned by each node on the edge into which one new node will be inserted
        std::set<unsigned> elems_containing_node_A = p_node_A->rGetContainingElementIndices();
        std::set<unsigned> elems_containing_node_B = p_node_B->rGetContainingElementIndices();


        c_vector<double, SPACE_DIM> position_a = p_node_A->rGetLocation();
        c_vector<double, SPACE_DIM> position_b = p_node_B->rGetLocation();

        c_vector<double, SPACE_DIM> a_to_b = GetVectorFromAtoB(position_a, position_b);

        // Find the location of the intersection
        double determinant = a_to_b[0]*axisOfDivision[1] - a_to_b[1]*axisOfDivision[0];

        c_vector<double, SPACE_DIM> moved_centroid = position_a + GetVectorFromAtoB(position_a, centroid); // allow for periodicity

        double alpha = (moved_centroid[0]*a_to_b[1] - position_a[0]*a_to_b[1]
                        -moved_centroid[1]*a_to_b[0] + position_a[1]*a_to_b[0])/determinant;

        c_vector<double, SPACE_DIM> intersection = moved_centroid + alpha*axisOfDivision;

        /*
         * If then new node is too close to one of the edge nodes, then reposition it
         * a distance mCellRearrangementRatio*mCellRearrangementThreshold further along the edge.
         */
        c_vector<double, SPACE_DIM> a_to_intersection = this->GetVectorFromAtoB(position_a, intersection);
        if (norm_2(a_to_intersection) < mCellRearrangementThreshold)
        {
            intersection = position_a + mCellRearrangementRatio*mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
        }

        c_vector<double, SPACE_DIM> b_to_intersection = this->GetVectorFromAtoB(position_b, intersection);
        if (norm_2(b_to_intersection) < mCellRearrangementThreshold)
        {
            intersection = position_b - mCellRearrangementRatio*mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
        }


        /*
         * The new node is boundary node if the 2 nodes are bounary nodes and the elements dont look like
         *   ___A___
         *  |   |   |
         *  |___|___|
         *      B
         */
        bool is_boundary = false;
        if (p_node_A->IsBoundaryNode() && p_node_B->IsBoundaryNode())
        {
            if (elems_containing_node_A.size() !=2 ||
                elems_containing_node_B.size() !=2 ||
                elems_containing_node_A != elems_containing_node_B)
            {
                is_boundary = true;
            }
        }

        // Add a new node to the mesh at the location of the intersection
        unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, is_boundary, intersection[0], intersection[1]));
        nodes_added++;

        // Now make sure node is added to neighbouring elements

        // Find common elements
        std::set<unsigned> shared_elements;
        std::set_intersection(elems_containing_node_A.begin(),
                              elems_containing_node_A.end(),
                              elems_containing_node_B.begin(),
                              elems_containing_node_B.end(),
                              std::inserter(shared_elements, shared_elements.begin()));

        // Iterate over common elements
        for (std::set<unsigned>::iterator iter = shared_elements.begin();
             iter != shared_elements.end();
             ++iter)
        {
            // Find which node has the lower local index in this element
            unsigned local_indexA = this->GetElement(*iter)->GetNodeLocalIndex(p_node_A->GetIndex());
            unsigned local_indexB = this->GetElement(*iter)->GetNodeLocalIndex(p_node_B->GetIndex());

            unsigned index = local_indexB;
            if (local_indexB > local_indexA)
            {
                index = local_indexA;
            }
            if ((local_indexA == 0) && (local_indexB == this->GetElement(*iter)->GetNumNodes()-1))
            {
                index = local_indexB;
            }
            if ((local_indexB == 0) && (local_indexA == this->GetElement(*iter)->GetNumNodes()-1))
            {
                index = local_indexA;
            }
            // Add new node to this element
            this->GetElement(*iter)->AddNode(index, this->GetNode(new_node_global_index));
        }
        // Store index of new node
        division_node_global_indices.push_back(new_node_global_index);

    }

    // Now call DivideElement() to divide the element using the new nodes
    unsigned new_element_index = DivideElement(pElement,
                                               pElement->GetNodeLocalIndex(division_node_global_indices[0]),
                                               pElement->GetNodeLocalIndex(division_node_global_indices[1]),
                                               placeOriginalElementBelow);
    return new_element_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongShortAxis(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                                                                bool placeOriginalElementBelow)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Find the short axis of the element
    c_vector<double, SPACE_DIM> short_axis = GetShortAxisOfElement(pElement->GetIndex());

    unsigned new_element_index = DivideElementAlongGivenAxis(pElement, short_axis, placeOriginalElementBelow);
    return new_element_index;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                                                  unsigned nodeAIndex,
                                                                  unsigned nodeBIndex,
                                                                  bool placeOriginalElementBelow)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Sort nodeA and nodeB such that nodeBIndex > nodeAindex
    assert(nodeBIndex != nodeAIndex);

    unsigned node1_index = (nodeAIndex < nodeBIndex) ? nodeAIndex : nodeBIndex; // low index
    unsigned node2_index = (nodeAIndex < nodeBIndex) ? nodeBIndex : nodeAIndex; // high index

    // Store the number of nodes in the element (this changes when nodes are deleted from the element)
    unsigned num_nodes = pElement->GetNumNodes();

    // Copy the nodes in this element
    std::vector<Node<SPACE_DIM>*> nodes_elem;
    for (unsigned i=0; i<num_nodes; i++)
    {
        nodes_elem.push_back(pElement->GetNode(i));
    }

    // Get the index of the new element
    unsigned new_element_index;
    if (mDeletedElementIndices.empty())
    {
        new_element_index = this->mElements.size();
    }
    else
    {
        new_element_index = mDeletedElementIndices.back();
        mDeletedElementIndices.pop_back();
        delete this->mElements[new_element_index];
    }

    // Add the new element to the mesh
    AddElement(new VertexElement<ELEMENT_DIM,SPACE_DIM>(new_element_index, nodes_elem));

    /**
     * Remove the correct nodes from each element. If placeOriginalElementBelow is true,
     * place the original element below (in the y direction) the new element; otherwise,
     * place it above.
     */

    /// Find lowest element \todo this could be more efficient
    double height_midpoint_1 = 0.0;
    double height_midpoint_2 = 0.0;
    unsigned counter_1 = 0;
    unsigned counter_2 = 0;

    for (unsigned i=0; i<num_nodes; i++)
    {
        if (i>=node1_index && i<=node2_index)
        {
            height_midpoint_1 += pElement->GetNode(i)->rGetLocation()[1];
            counter_1++;
        }
        if (i<=node1_index || i>=node2_index)
        {
            height_midpoint_2 += pElement->GetNode(i)->rGetLocation()[1];
            counter_2++;
        }
    }
    height_midpoint_1 /= (double)counter_1;
    height_midpoint_2 /= (double)counter_2;

    for (unsigned i=num_nodes; i>0; i--)
    {
        if (i-1 < node1_index || i-1 > node2_index)
        {
            if (height_midpoint_1 < height_midpoint_2)
            {
                if (placeOriginalElementBelow)
                {
                    pElement->DeleteNode(i-1);
                }
                else
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
            }
            else
            {
                if (placeOriginalElementBelow)
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
                else
                {
                    pElement->DeleteNode(i-1);
                }
            }
        }
        else if (i-1 > node1_index && i-1 < node2_index)
        {
            if (height_midpoint_1 < height_midpoint_2)
            {
                if (placeOriginalElementBelow)
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
                else
                {
                    pElement->DeleteNode(i-1);
                }
            }
            else
            {
                if (placeOriginalElementBelow)
                {
                    pElement->DeleteNode(i-1);
                }
                else
                {
                    this->mElements[new_element_index]->DeleteNode(i-1);
                }
            }
        }
    }

    return new_element_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT3Swap(Node<SPACE_DIM>* pNode, unsigned elementIndex)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    #define COVERAGE_IGNORE
    assert(SPACE_DIM == 2); // only works in 2D at present
    assert(ELEMENT_DIM == SPACE_DIM);
    #undef COVERAGE_IGNORE

    // Check pNode is a boundary node
    assert(pNode->IsBoundaryNode());

    // Store the index of the elements containing the intersecting node
    std::set<unsigned> elements_containing_intersecting_node = pNode->rGetContainingElementIndices();

    // Get the local index of the node in the intersected element after which the new node is to be added
    unsigned node_A_local_index = GetLocalIndexForElementEdgeClosestToPoint(pNode->rGetLocation(), elementIndex);

    // Get current node location
    c_vector<double, SPACE_DIM> node_location = pNode->rGetModifiableLocation();

    // Get element
    VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    /*
     * Get the nodes at either end of the edge to be divided
     */
    unsigned vertexA_index = p_element->GetNodeGlobalIndex(node_A_local_index);
    unsigned vertexB_index = p_element->GetNodeGlobalIndex((node_A_local_index+1)%num_nodes);


    // Check these nodes are also boundary nodes
    assert(this->mNodes[vertexA_index]->IsBoundaryNode());
    assert(this->mNodes[vertexB_index]->IsBoundaryNode());

    /*
     * Get the nodes at either end of the edge to be divided and calculate intersection
     */
    c_vector<double, SPACE_DIM> vertexA = p_element->GetNodeLocation(node_A_local_index);
    c_vector<double, SPACE_DIM> vertexB = p_element->GetNodeLocation((node_A_local_index+1)%num_nodes);
    c_vector<double, SPACE_DIM> vector_a_to_point = this->GetVectorFromAtoB(vertexA, node_location);

    c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);

    if (norm_2(vector_a_to_b) < 2.0*mCellRearrangementRatio*mCellRearrangementRatio*mCellRearrangementThreshold)
    {
        EXCEPTION("Trying to merge a node onto an edge which is too small.");
    }

    c_vector<double, SPACE_DIM> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);
    c_vector<double, SPACE_DIM> intersection = vertexA + edge_ab_unit_vector*inner_prod(vector_a_to_point, edge_ab_unit_vector);

    /*
     * If the intersection is within mCellRearrangementRatio^2*mCellRearangementThreshold of vertexA or vertexB move it
     * mCellRearrangementRatio^2*mCellRearrangementThreshold away.
     *
     * Note: this distance so that there is always enough room for new nodes (if necessary).
     */
    if (norm_2(intersection - vertexA) < mCellRearrangementRatio*mCellRearrangementRatio*mCellRearrangementThreshold)
    {
        intersection = vertexA + mCellRearrangementRatio*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
    }
    if (norm_2(intersection - vertexB) < mCellRearrangementRatio*mCellRearrangementRatio*mCellRearrangementThreshold)
    {
        intersection = vertexB - mCellRearrangementRatio*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
    }

    if (pNode->GetNumContainingElements() == 1)
    {
        // Get the index of the element containing the intersecting node
        unsigned intersecting_element_index = *elements_containing_intersecting_node.begin();

        // Get element
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_intersecting_element = this->GetElement(intersecting_element_index);

        unsigned local_index = p_intersecting_element->GetNodeLocalIndex(pNode->GetIndex());
        unsigned next_node = p_intersecting_element->GetNodeGlobalIndex((local_index + 1)%(p_intersecting_element->GetNumNodes()));
        unsigned previous_node = p_intersecting_element->GetNodeGlobalIndex((local_index + p_intersecting_element->GetNumNodes() - 1)%(p_intersecting_element->GetNumNodes()));

        // Check to see if the nodes adjacent to the intersecting node are contained in the intersected element VertexA and VertexB
        if (next_node == vertexA_index || previous_node == vertexA_index || next_node == vertexB_index || previous_node == vertexB_index)
        {
            /*
             *  From          To
             *   _             _
             *    |  <---       |
             *    |  /\         |\
             *    | /  \        | \
             *   _|/____\      _|__\
             *
             */

            /*
             *  The edge goes from vertexA--vertexB to vertexA--pNode--vertexB
             */
            // Move original node
            pNode->rGetModifiableLocation() = intersection;

            // Add the moved nodes to the element (this also updates the node)
            this->GetElement(elementIndex)->AddNode(node_A_local_index, pNode);

            // Check the nodes are updated correctly
            assert(pNode->GetNumContainingElements() == 2);
        }
        else
        {
            /*
             *  From          To
             *   ____        _______
             *                 / \
             *    /\   ^      /   \
             *   /  \  |
             *
             */

            /*
             *  The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
             */

            // Move original node
            pNode->rGetModifiableLocation() = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

            c_vector<double, SPACE_DIM> new_node_location = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

            // Add new node whic will always be a boundary node
            unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

            // Add the moved and new nodes to the element (this also updates the node)
            this->GetElement(elementIndex)->AddNode(node_A_local_index, pNode);
            this->GetElement(elementIndex)->AddNode(node_A_local_index, this->mNodes[new_node_global_index]);

            // Add the new node to the original element contiang pNode (this also updates the node)
            this->GetElement(intersecting_element_index)->AddNode(this->GetElement(intersecting_element_index)->GetNodeLocalIndex(pNode->GetIndex()), this->mNodes[new_node_global_index]);

            // Check the nodes are updated correctly
            assert(pNode->GetNumContainingElements() == 2);
            assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);
        }
    }
    else if (pNode->GetNumContainingElements() == 2)
    {
        // Find the nodes contained in elements containing the intersecting node
        std::set<unsigned>::const_iterator it = elements_containing_intersecting_node.begin();
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_1 = this->GetElement(*it);
        it++;
        VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_2 = this->GetElement(*it);

        unsigned local_index_1 = p_element_1->GetNodeLocalIndex(pNode->GetIndex());
        unsigned next_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + 1)%(p_element_1->GetNumNodes()));
        unsigned previous_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + p_element_1->GetNumNodes() - 1)%(p_element_1->GetNumNodes()));
        unsigned local_index_2 = p_element_2->GetNodeLocalIndex(pNode->GetIndex());
        unsigned next_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + 1)%(p_element_2->GetNumNodes()));
        unsigned previous_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + p_element_2->GetNumNodes() - 1)%(p_element_2->GetNumNodes()));

        // Check to see if the nodes adjacent to the intersecting node are contained in the intersected element VertexA and VertexB
        if (next_node_1 == vertexA_index || previous_node_1 == vertexA_index || next_node_2 == vertexA_index || previous_node_2 == vertexA_index)
        {
            /*
             *  From          To
             *   _ B              _ B
             *    |   <---         |
             *    |   /|\          |\
             *    |  / | \         | \
             *    | /  |  \        |\ \
             *   _|/___|___\      _|_\_\
             *   A                  A
             */

            /*
             *  The edge goes from vertexA--vertexB to vertexA--pNode--new_node--vertexB
             */

            // Move original node and change to non-boundary node
            pNode->rGetModifiableLocation() = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
            pNode->SetAsBoundaryNode(false);

             c_vector<double, SPACE_DIM> new_node_location = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

            // Add new node which will always be a boundary node
            unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

            // Add the moved nodes to the element (this also updates the node)
            this->GetElement(elementIndex)->AddNode(node_A_local_index, this->mNodes[new_node_global_index]);
            this->GetElement(elementIndex)->AddNode(node_A_local_index, pNode);

            // Add the new nodes to the original elements contiang pNode (this also updates the node)
            if (next_node_1 == previous_node_2)
            {
                p_element_1->AddNode((local_index_1 + p_element_1->GetNumNodes() - 1)%(p_element_1->GetNumNodes()), this->mNodes[new_node_global_index]);
            }
            else
            {
                assert(next_node_2 == previous_node_1);

                p_element_2->AddNode((local_index_2 + p_element_2->GetNumNodes() - 1)%(p_element_2->GetNumNodes()), this->mNodes[new_node_global_index]);
            }


            // Check the nodes are updated correctly
            assert(pNode->GetNumContainingElements() == 3);
            assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);
        }
        else if (next_node_1 == vertexB_index || previous_node_1 == vertexB_index || next_node_2 == vertexB_index || previous_node_2 == vertexB_index)
        {
            /*
             *  From          To
             *   _B_________      _B____
             *    |\   |   /       | / /
             *    | \  |  /        |/ /
             *    |  \ | /         | /
             *    |   \|/          |/
             *   _|   <---        _|
             *    A
             */

            /*
             *  The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
             */

            // Move original node and change to non-boundary node
            pNode->rGetModifiableLocation() = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
            pNode->SetAsBoundaryNode(false);

             c_vector<double, SPACE_DIM> new_node_location = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

            // Add new node which will always be a boundary node
            unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));

            // Add the moved nodes to the element (this also updates the node)
            this->GetElement(elementIndex)->AddNode(node_A_local_index, pNode);
            this->GetElement(elementIndex)->AddNode(node_A_local_index, this->mNodes[new_node_global_index]);

            // Add the new nodes to the original elements contiang pNode (this also updates the node)
            if (next_node_1 == previous_node_2)
            {
                p_element_2->AddNode(local_index_2, this->mNodes[new_node_global_index]);
            }
            else
            {
                assert(next_node_2 == previous_node_1);

                p_element_1->AddNode(local_index_1, this->mNodes[new_node_global_index]);
            }

            // Check the nodes are updated correctly
            assert(pNode->GetNumContainingElements() == 3);
            assert(this->mNodes[new_node_global_index]->GetNumContainingElements() == 2);
        }
        else
        {
            /*
             *  From          To
             *   _____         _______
             *                  / | \
             *    /|\   ^      /  |  \
             *   / | \  |
             *
             */

            /*
             *  The edge goes from vertexA--vertexB to vertexA--new_node_1--pNode--new_node_2--vertexB
             */

            // Move original node and change to non-boundary node
            pNode->rGetModifiableLocation() = intersection;
            pNode->SetAsBoundaryNode(false);

            c_vector<double, SPACE_DIM> new_node_1_location = intersection - mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;
            c_vector<double, SPACE_DIM> new_node_2_location = intersection + mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

            // Add new nodes which will always be boundary nodes
            unsigned new_node_1_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_1_location[0], new_node_1_location[1]));
            unsigned new_node_2_global_index = this->AddNode(new Node<SPACE_DIM>(0, true, new_node_2_location[0], new_node_2_location[1]));

            // Add the moved and new nodes to the element (this also updates the node)
            this->GetElement(elementIndex)->AddNode(node_A_local_index, this->mNodes[new_node_2_global_index]);
            this->GetElement(elementIndex)->AddNode(node_A_local_index, pNode);
            this->GetElement(elementIndex)->AddNode(node_A_local_index, this->mNodes[new_node_1_global_index]);

            // Add the new nodes to the original elements contiang pNode (this also updates the node)
            if (next_node_1 == previous_node_2)
            {
                p_element_1->AddNode((local_index_1 + p_element_1->GetNumNodes() - 1)%(p_element_1->GetNumNodes()), this->mNodes[new_node_2_global_index]);
                p_element_2->AddNode(local_index_2, this->mNodes[new_node_1_global_index]);
            }
            else
            {
                assert(next_node_2 == previous_node_1);

                p_element_1->AddNode(local_index_1, this->mNodes[new_node_1_global_index]);
                p_element_2->AddNode((local_index_2 + p_element_2->GetNumNodes() - 1)%(p_element_2->GetNumNodes()), this->mNodes[new_node_2_global_index]);
            }

            // Check the nodes are updated correctly
            assert(pNode->GetNumContainingElements() == 3);
            assert(this->mNodes[new_node_1_global_index]->GetNumContainingElements() == 2);
            assert(this->mNodes[new_node_2_global_index]->GetNumContainingElements() == 2);
        }
    }
    else
    {
        EXCEPTION("Trying to merge a node, contained in more than 2 elements, into another element, this is not possible with the vertex mesh.");
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class MutableVertexMesh<1,1>;
template class MutableVertexMesh<1,2>;
template class MutableVertexMesh<1,3>;
template class MutableVertexMesh<2,2>;
template class MutableVertexMesh<2,3>;
template class MutableVertexMesh<3,3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableVertexMesh);
