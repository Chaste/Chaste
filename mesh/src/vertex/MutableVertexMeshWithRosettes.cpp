/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "MutableVertexMeshWithRosettes.hpp"
#include "UblasCustomFunctions.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::MutableVertexMeshWithRosettes(std::vector<Node<SPACE_DIM>*> nodes,
                                                                                     std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements,
                                                                                     double cellRearrangementThreshold,
                                                                                     double t2Threshold,
                                                                                     double cellRearrangementRatio,
                                                                                     double protorosetteFormationProbability,
                                                                                     double protorosetteResolutionProbabilityPerTimestep,
                                                                                     double rosetteResolutionProbabilityPerTimestep)
        : MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>(nodes, vertexElements, cellRearrangementThreshold, t2Threshold, cellRearrangementRatio),
          mProtorosetteFormationProbability(protorosetteFormationProbability),
          mProtorosetteResolutionProbabilityPerTimestep(protorosetteResolutionProbabilityPerTimestep),
          mRosetteResolutionProbabilityPerTimestep(rosetteResolutionProbabilityPerTimestep)
{
    // Make a globally accessible random number generator, as many methods use it
    mpGen = RandomNumberGenerator::Instance();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::MutableVertexMeshWithRosettes()
        : MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::~MutableVertexMeshWithRosettes()
{
    this->Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::HandleHighOrderJunctions(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    unsigned node_a_rank = pNodeA->rGetContainingElementIndices().size();
    unsigned node_b_rank = pNodeB->rGetContainingElementIndices().size();

    if ((node_a_rank > 3) && (node_b_rank > 3))
    {
        // The code can't handle this case
        EXCEPTION("Both nodes involved in a swap event are contained in more than three elements");
    }
    else if ((node_a_rank <= 3) && (node_b_rank <= 3))
    {
        // This method shouldn't have been called
        EXCEPTION("Neither node is high order");
    }
    else // the rosette degree should increase in this case
    {
        assert(node_a_rank > 3 || node_b_rank > 3);
        this->PerformRosetteRankIncrease(pNodeA, pNodeB);
        this->RemoveDeletedNodes();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::HandleAdditionalRemodellingBehaviour(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, std::set<unsigned> elem_indices, unsigned case_num)
{
    /*
     * The case number refers to MutableVertexMesh::IdentifySwapType() switch case
     */
    if( case_num == 4 )
    {
        /*
         * We choose either to perform protorosette formation (node merge), or else a T1 swap, based on a random number
         */
        if (mProtorosetteFormationProbability >= mpGen->ranf())
        {
            this->PerformNodeMerge(pNodeA, pNodeB);
            this->RemoveDeletedNodes();
        }
        else
        {
            this->PerformT1Swap(pNodeA, pNodeB, elem_indices);
        }
    }
    else
    {
        EXCEPTION("No functionality for this case yet");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::HandleAdditionalReMeshingBehaviour()
{
    /**
     * Check for protorosettes and rosettes, to allow resolution if necessary.  Unlike the two methods above,
     * we do not recheck the mesh.  Instead, we look through the mesh to check for all high-order nodes and
     * allow for possible resolution all at once.  We do not allow a Rosette to undergo multiple resolution
     * events in a single timestep.
     */
    this->CheckForRosettes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::PerformRosetteRankIncrease(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    /*
     * One of the nodes will have 3 containing element indices, the other
     * will have at least four. We first identify which node is which.
     */

    unsigned node_a_index = pNodeA->GetIndex();
    unsigned node_b_index = pNodeB->GetIndex();

    unsigned node_a_rank = pNodeA->rGetContainingElementIndices().size();
    unsigned node_b_rank = pNodeB->rGetContainingElementIndices().size();

    unsigned lo_rank_index = (node_a_rank < node_b_rank) ? node_a_index : node_b_index;
    unsigned hi_rank_index = (node_a_rank < node_b_rank) ? node_b_index : node_a_index;

    // Get pointers to the nodes, sorted by index
    Node<SPACE_DIM>* p_lo_rank_node = this->GetNode(lo_rank_index);
    Node<SPACE_DIM>* p_hi_rank_node = this->GetNode(hi_rank_index);

    // Find the sets of elements containing each of the nodes, sorted by index
    std::set<unsigned> lo_rank_elem_indices = p_lo_rank_node->rGetContainingElementIndices();
    std::set<unsigned> hi_rank_elem_indices = p_hi_rank_node->rGetContainingElementIndices();

    /**
     * The picture below shows the situation we are in.  The central node (marked with an 'X')
     * is contained in (at least) four elements already (A, B, C and D), and this is the
     * rosette node which we have designated as hi_rank_node.
     *
     * The node shared by elements C, D and E has come within the cell rearrangement threshold
     * of the rosette node in order for this method to have been called.  We have designated
     * this node lo_rank_node.
     *
     * We now 'merge' hi_rank_node and lo_rank_node, but keep hi_rank_node where it is
     * (which is why we don't call PerformNodeMerge(), as that would move both nodes to
     * their average location).
     *
     * To accomplish this merge, we need do nothing to elements A and B.  We remove lo_rank_node
     * from elements C and D, and we replace lo_rank_node by hi_rank_node in element E.
     *
     *
     *      \  A  /
     *       \   /
     *        \ /
     *     B   X   C
     *        / \______
     *       /   \
     *      /  D  \  E
     */

    for (std::set<unsigned>::const_iterator it = lo_rank_elem_indices.begin();
         it != lo_rank_elem_indices.end();
         ++it)
    {
        // Find the local index of lo_rank_node in this element
        unsigned lo_rank_local_index = this->mElements[*it]->GetNodeLocalIndex(lo_rank_index);
        assert(lo_rank_local_index < UINT_MAX); // double check this element contains lo_rank_node

        /*
         * If this element already contains the hi_rank_node, we are in the situation of elements
         * C and D above, so we just remove lo_rank_node.
         *
         * Otherwise, we are in element E, so we must replace lo_rank_node with high-rank node,
         * and remove it from mNodes.
         *
         * We can check whether hi_rank_node is in this element using the set::count() method.
         */

        if (hi_rank_elem_indices.count(*it) > 0)
        {
            // Delete lo_rank_node from current element
            this->mElements[*it]->DeleteNode(lo_rank_local_index);
        }
        else
        {
            // Update lo_rank_node with all information (including index and location) of hi_rank_node
            this->mElements[*it]->UpdateNode(lo_rank_local_index, p_hi_rank_node);
        }
    }

    // Tidy up the mesh by ensuring the global instance of lo_rank_node is deleted
    assert(!(this->mNodes[lo_rank_index]->IsDeleted()));
    this->mNodes[lo_rank_index]->MarkAsDeleted();
    this->mDeletedNodeIndices.push_back(lo_rank_index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::PerformProtorosetteResolution(Node<SPACE_DIM>* pProtorosetteNode)
{
    // Double check we are dealing with a protorosette
    assert(pProtorosetteNode->rGetContainingElementIndices().size() == 4);

    // Get random number (0, 1, 2 or 3), as the resolution axis is assumed to be random
    unsigned random_elem_increment = mpGen->randMod(4);

    // Find global indices of elements around the protorosette node
    std::set<unsigned> protorosette_node_containing_elem_indices = pProtorosetteNode->rGetContainingElementIndices();

    // Select random element by advancing iterator a random number times
    std::set<unsigned>::const_iterator elem_index_iter(protorosette_node_containing_elem_indices.begin());
    advance(elem_index_iter, random_elem_increment);

    /**
     * Ordering elements as follows:
     *
     *      \  A  /
     *       \   /
     *        \ /
     *     B   X   D
     *        / \
     *       /   \
     *      /  C  \
     *
     * Element A is the randomly chosen element.  Element C, which is directly
     * opposite A, will end up separated from A, while the two elements B and
     * D which start adjacent to A will end up sharing a common edge:
     *
     *      \  A  /
     *       \   /
     *        \ /
     *         |
     *    B    |    D
     *         |
     *        / \
     *       /   \
     *      /  C  \
     *
     */

    /*
     * We need to find the global indices of elements B, C and D.  We do this with set intersections.
     */

    unsigned elem_a_idx = *elem_index_iter;
    unsigned elem_b_idx = UINT_MAX;
    unsigned elem_c_idx = UINT_MAX;
    unsigned elem_d_idx = UINT_MAX;

    // Get pointer to element we've chosen at random (element A)
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_a = this->GetElement(elem_a_idx);

    // Get all necessary info about element A and the protorosette node
    unsigned num_nodes_elem_a = p_elem_a->GetNumNodes();
    unsigned protorosette_node_global_idx = pProtorosetteNode->GetIndex();
    unsigned protorosette_node_local_idx = p_elem_a->GetNodeLocalIndex(protorosette_node_global_idx);

    // Find global indices of previous (cw) and next (ccw) nodes, locally, from the protorosette node, in element A
    unsigned prev_node_global_idx = p_elem_a->GetNodeGlobalIndex((protorosette_node_local_idx + num_nodes_elem_a - 1) % num_nodes_elem_a);
    unsigned next_node_global_idx = p_elem_a->GetNodeGlobalIndex((protorosette_node_local_idx + 1) % num_nodes_elem_a);

    // Get the set of elements the previous and next nodes are contained in
    Node<SPACE_DIM>* p_prev_node = this->GetNode(prev_node_global_idx);
    Node<SPACE_DIM>* p_next_node = this->GetNode(next_node_global_idx);
    std::set<unsigned> prev_node_elem_indices = p_prev_node->rGetContainingElementIndices();
    std::set<unsigned> next_node_elem_indices = p_next_node->rGetContainingElementIndices();

    // Perform set intersections with the set of element indices which the protorosette node is contained in
    std::set<unsigned> intersection_with_prev;
    std::set<unsigned> intersection_with_next;

    // This intersection should contain just global indices for elements A and B
    std::set_intersection(protorosette_node_containing_elem_indices.begin(),
                          protorosette_node_containing_elem_indices.end(),
                          prev_node_elem_indices.begin(),
                          prev_node_elem_indices.end(),
                          std::inserter(intersection_with_prev, intersection_with_prev.begin()));

    // This intersection should contain just global indices for elements A and D
    std::set_intersection(protorosette_node_containing_elem_indices.begin(),
                          protorosette_node_containing_elem_indices.end(),
                          next_node_elem_indices.begin(),
                          next_node_elem_indices.end(),
                          std::inserter(intersection_with_next, intersection_with_next.begin()));

    assert(intersection_with_prev.size() == 2);
    assert(intersection_with_next.size() == 2);

    // Get global index of element B
    if (*intersection_with_prev.begin() != elem_a_idx)
    {
        elem_b_idx = *(intersection_with_prev.begin());
    }
    else
    {
        elem_b_idx = *(++(intersection_with_prev.begin()));
    }
    assert(elem_b_idx < UINT_MAX);

    // Get global index of element D
    if (*intersection_with_next.begin() != elem_a_idx)
    {
        elem_d_idx = *(intersection_with_next.begin());
    }
    else
    {
        elem_d_idx = *(++(intersection_with_next.begin()));
    }
    assert(elem_d_idx < UINT_MAX);

    // By elimination, the remaining unassigned index in the original set must be global index of element C
    for (elem_index_iter = protorosette_node_containing_elem_indices.begin();
         elem_index_iter != protorosette_node_containing_elem_indices.end();
         ++elem_index_iter)
    {
        if ( (*elem_index_iter != elem_a_idx) && (*elem_index_iter != elem_b_idx) && (*elem_index_iter != elem_d_idx) )
        {
            elem_c_idx = *elem_index_iter;
        }
    }
    assert(elem_c_idx < UINT_MAX);

    /**
     * Next, we compute where to place the two nodes which will replace the single protorosette node.
     *
     * We place each node along the line joining the protorosette node to the centroid of element which will contain it,
     * and the distance along this line is half of the swap distance.
     *
     * To do this, we will move the existing protorosette node in to element A, and create a new node in element C.  We
     * then need to tidy up the nodes by adding the new node to elements B, C and D, and removing the protorosette node
     * from element C.
     *
     * NOTE: as the protorosette node was not necessarily on the line joining the centroids of elements A and C, unlike
     * in a T1 swap, the new node locations will not necessarily be the full swap distance apart.
     */

    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_b = this->GetElement(elem_b_idx);
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_c = this->GetElement(elem_c_idx);
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_d = this->GetElement(elem_d_idx);

    double swap_distance = (this->mCellRearrangementRatio) * (this->mCellRearrangementThreshold);

    // Get normalized vectors to centre of elements A and B from protorosette node
    c_vector<double, SPACE_DIM> node_to_elem_a_centre = this->GetCentroidOfElement(elem_a_idx) - pProtorosetteNode->rGetLocation();
    node_to_elem_a_centre /= norm_2(node_to_elem_a_centre);

    c_vector<double, SPACE_DIM> node_to_elem_c_centre = this->GetCentroidOfElement(elem_c_idx) - pProtorosetteNode->rGetLocation();
    node_to_elem_c_centre /= norm_2(node_to_elem_c_centre);

    // Calculate new node locations
    c_vector<double, SPACE_DIM> new_location_of_protorosette_node = pProtorosetteNode->rGetLocation() + (0.5 * swap_distance) * node_to_elem_a_centre;
    c_vector<double, SPACE_DIM> location_of_new_node 			  = pProtorosetteNode->rGetLocation() + (0.5 * swap_distance) * node_to_elem_c_centre;

    // Move protorosette node to new location
    pProtorosetteNode->rGetModifiableLocation() = new_location_of_protorosette_node;

    // Create new node in correct location
    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(this->GetNumNodes(), location_of_new_node, false));
    Node<SPACE_DIM>* p_new_node = this->GetNode(new_node_global_index);

    /**
     * Here, we add the new node to elements B, C and D.
     *
     * The method AddNode() takes the local index of a node i, where the new node is to be inserted between nodes
     * i and i+1, so we need to find the local idx of the node directly before (clockwise from) where we wish to
     * insert it.
     *
     * For elements C and D, we just need the local index of the protorosette node, but for element B we need the one
     * before, modulo number of elements.
     */
    unsigned local_idx_elem_b = p_elem_b->GetNodeLocalIndex(protorosette_node_global_idx);
    local_idx_elem_b = (local_idx_elem_b + p_elem_b->GetNumNodes() - 1) % p_elem_b->GetNumNodes();
    unsigned local_idx_elem_c = p_elem_c->GetNodeLocalIndex(protorosette_node_global_idx);
    unsigned local_idx_elem_d = p_elem_d->GetNodeLocalIndex(protorosette_node_global_idx);

    p_elem_b->AddNode(p_new_node, local_idx_elem_b);
    p_elem_c->AddNode(p_new_node, local_idx_elem_c);
    p_elem_d->AddNode(p_new_node, local_idx_elem_d);

    // All that is left is to remove the original protorosette node from element C
    p_elem_c->DeleteNode(p_elem_c->GetNodeLocalIndex(protorosette_node_global_idx));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::PerformRosetteRankDecrease(Node<SPACE_DIM>* pRosetteNode)
{
    unsigned rosette_rank = pRosetteNode->rGetContainingElementIndices().size();

    // Double check we're dealing with a rosette
    assert(rosette_rank > 4);

    // Get random number in [0, 1, ..., n) where n is rank of rosette, as the resolution axis is assumed to be random
    unsigned random_elem_increment = mpGen->randMod(rosette_rank);

    // Find global indices of elements around the protorosette node
    std::set<unsigned> rosette_node_containing_elem_indices = pRosetteNode->rGetContainingElementIndices();

    // Select random element by advancing iterator a random number times
    std::set<unsigned>::const_iterator elem_index_iter(rosette_node_containing_elem_indices.begin());
    advance(elem_index_iter, random_elem_increment);

    /**
     * We have now picked a vertex element at random from the rosette.
     * This element will now be disconnected from the rosette in a manner
     * analogous to performing a T1 swap.
     *
     * The node at the centre of the rosette will not move, and a new
     * node will be created a suitable distance away, along a line joining
     * the rosette node and the centroid of the element we have randomly
     * selected to move.
     *
     * Ordering of elements is as follows:
     *
     *      \  S  /
     *       \   /
     *     N  \ /  P
     *    -----X-----
     *        / \
     *       /   \
     *      /     \
     *
     * where element S is the selected element,
     * N is the next (counterclockwise) element from S, and
     * P is the previous (clockwise) element from S.
     *
     * Elements N and P will end up sharing an edge:
     *
     *      \  S  /
     *       \   /
     *        \ /
     *      N  |  P
     *         |
     *    -----*-----
     *        / \
     *       /   \
     *      /     \
     *
     */

    /*
     * We need to find the global indices of elements N and P.  We do this with set intersections.
     */

    // Get the vertex element S (which we randomly selected)
    unsigned elem_s_idx = *elem_index_iter;
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_s = this->GetElement(elem_s_idx);

    unsigned elem_n_idx = UINT_MAX;
    unsigned elem_p_idx = UINT_MAX;

    // Get all necessary info about element S and the rosette node
    unsigned num_nodes_elem_s = p_elem_s->GetNumNodes();
    unsigned rosette_node_global_idx = pRosetteNode->GetIndex();
    unsigned rosette_node_local_idx = p_elem_s->GetNodeLocalIndex(rosette_node_global_idx);

    // Find global indices of previous (cw) and next (ccw) nodes, locally, from the rosette node, in element S
    unsigned prev_node_global_idx = p_elem_s->GetNodeGlobalIndex((rosette_node_local_idx + num_nodes_elem_s - 1) % num_nodes_elem_s);
    unsigned next_node_global_idx = p_elem_s->GetNodeGlobalIndex((rosette_node_local_idx + 1) % num_nodes_elem_s);

    // Get the set of elements that the previous and next nodes are contained in
    Node<SPACE_DIM>* p_prev_node = this->GetNode(prev_node_global_idx);
    Node<SPACE_DIM>* p_next_node = this->GetNode(next_node_global_idx);
    std::set<unsigned> prev_node_elem_indices = p_prev_node->rGetContainingElementIndices();
    std::set<unsigned> next_node_elem_indices = p_next_node->rGetContainingElementIndices();

    // Perform set intersections with the set of element indices that the rosette node is contained in
    std::set<unsigned> intersection_with_prev;
    std::set<unsigned> intersection_with_next;

    // This intersection should contain just global indices for elements S and N
    std::set_intersection(rosette_node_containing_elem_indices.begin(),
                          rosette_node_containing_elem_indices.end(),
                          prev_node_elem_indices.begin(),
                          prev_node_elem_indices.end(),
                          std::inserter(intersection_with_prev, intersection_with_prev.begin()));

    // This intersection should contain just global indices for elements S and P
    std::set_intersection(rosette_node_containing_elem_indices.begin(),
                          rosette_node_containing_elem_indices.end(),
                          next_node_elem_indices.begin(),
                          next_node_elem_indices.end(),
                          std::inserter(intersection_with_next, intersection_with_next.begin()));

    assert(intersection_with_prev.size() == 2);
    assert(intersection_with_next.size() == 2);

    // Get global index of element N
    if (*intersection_with_prev.begin() != elem_s_idx)
    {
        elem_n_idx = *intersection_with_prev.begin();
    }
    else
    {
        elem_n_idx = *(++(intersection_with_prev.begin()));
    }
    assert(elem_n_idx < UINT_MAX);

    // Get global index of element P
    if (*intersection_with_next.begin() != elem_s_idx)
    {
        elem_p_idx = *intersection_with_next.begin();
    }
    else
    {
        elem_p_idx = *(++(intersection_with_next.begin()));
    }
    assert(elem_p_idx < UINT_MAX);

    /**
     * Next, we compute where to place the new node which will separate the rosette node from element S.
     *
     * We place this node along the line joining the rosette node to the centroid of element S, at a distance of the
     * swap distance ( (rearrangement ratio) x (rearrangement threshold) )
     *
     * To do this, we create a new node in element S.  We then need to tidy up the nodes by adding the new node to
     * elements S, N and P, and by removing the rosette node from element S.
     */

    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_n = this->GetElement(elem_p_idx);
    VertexElement<ELEMENT_DIM,SPACE_DIM>* p_elem_p = this->GetElement(elem_n_idx);

    double swap_distance = (this->mCellRearrangementRatio) * (this->mCellRearrangementThreshold);

    // Calculate location of new node
    c_vector<double, 2> node_to_selected_elem = this->GetCentroidOfElement(elem_s_idx) - pRosetteNode->rGetLocation();
    node_to_selected_elem /= norm_2(node_to_selected_elem);
    c_vector<double, 2> new_node_location = pRosetteNode->rGetLocation() + (swap_distance * node_to_selected_elem);

    // Create new node in correct location
    unsigned new_node_global_index = this->AddNode(new Node<SPACE_DIM>(this->GetNumNodes(), new_node_location, false));
    Node<SPACE_DIM>* p_new_node = this->GetNode(new_node_global_index);

    /**
     * Here, we add the new node to elements S, N and P, and remove the rosette node from element S.
     *
     * The method AddNode() takes the local index of a node i, where the new node is to be inserted between nodes
     * i and i+1, so we need to find the local idx of the node directly before (clockwise from) where we wish to
     * insert it.
     *
     * For elements S and P, we just need the local index of the rosette node, but for element N we need the one
     * before, modulo number of elements.
     */

    // Add new node, and remove rosette node, from element S
    unsigned node_local_idx_in_elem_s = p_elem_s->GetNodeLocalIndex(rosette_node_global_idx);
    p_elem_s->AddNode(p_new_node, node_local_idx_in_elem_s);
    p_elem_s->DeleteNode(node_local_idx_in_elem_s);

    // Add new node to element N
    unsigned node_local_idx_in_elem_n = p_elem_n->GetNodeLocalIndex(rosette_node_global_idx);
    node_local_idx_in_elem_n = (node_local_idx_in_elem_n + p_elem_n->GetNumNodes() - 1) % p_elem_n->GetNumNodes();
    p_elem_n->AddNode(p_new_node, node_local_idx_in_elem_n);

    // Add new node to element P
    unsigned node_local_idx_in_elem_p = p_elem_p->GetNodeLocalIndex(rosette_node_global_idx);
    p_elem_p->AddNode(p_new_node, node_local_idx_in_elem_p);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::CheckForRosettes()
{
    // Rosettes are currently only allowed in two dimensions
    assert(ELEMENT_DIM==2 && SPACE_DIM==2);

    /**
     * First, we loop over each node and populate vectors of protorosette and rosette nodes which need to undergo
     * resolution.
     *
     * We do not perform the resolution events in this initial loop because the resolution events involve changing
     * nodes in the mesh.
     */

    // Vectors to store the nodes that need resolution events
    std::vector<Node<SPACE_DIM>* > protorosette_nodes;
    std::vector<Node<SPACE_DIM>* > rosette_nodes;

    // First loop in which we populate these vectors
    unsigned num_nodes = this->GetNumAllNodes();
    for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
    {
        Node<SPACE_DIM>* current_node = this->GetNode(node_idx);
        unsigned node_rank = current_node->rGetContainingElementIndices().size();

        if (node_rank < 4)
        {
            // Nothing to do if the node is not high-rank
            continue;
        }
        else if (node_rank == 4)
        {
            // For protorosette nodes, we check against a random number to decide if resolution is necessary
            if (mProtorosetteResolutionProbabilityPerTimestep >= mpGen->ranf())
            {
                protorosette_nodes.push_back(current_node);
            }
        }
        else // if (node_rank > 4)
        {
            // For rosette nodes, we check against a random number to decide if resolution is necessary
            if (mRosetteResolutionProbabilityPerTimestep >= mpGen->ranf())
            {
                rosette_nodes.push_back(current_node);
            }
        }
    }

    /**
     * Finally, we loop over the contents of each node vector and perform the necessary resolution events.
     *
     * Because each resolution event changes nodes, we include several assertions to catch possible unconsidered
     * behaviour.
     */

    // First, resolve any protorosettes
    for (unsigned node_idx = 0 ; node_idx < protorosette_nodes.size() ; node_idx++)
    {
        Node<SPACE_DIM>* current_node = protorosette_nodes[node_idx];

        // Verify that node has not been marked for deletion, and that it is still contained in four elements
        assert( !(current_node->IsDeleted()) );
        assert( current_node->rGetContainingElementIndices().size() == 4 );

        // Perform protorosette resolution
        this->PerformProtorosetteResolution(current_node);
    }

    // Finally, resolve any rosettes
    for (unsigned node_idx = 0 ; node_idx < rosette_nodes.size() ; node_idx++)
    {
        Node<SPACE_DIM>* current_node = rosette_nodes[node_idx];

        // Verify that node has not been marked for deletion, and that it is still contained in at least four elements
        assert( !(current_node->IsDeleted()) );
        assert( current_node->rGetContainingElementIndices().size() > 4 );

        // Perform protorosette resolution
        this->PerformRosetteRankDecrease(current_node);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteFormationProbability(double protorosetteFormationProbability)
{
    // Check that the new value is in [0,1]
    if (protorosetteFormationProbability < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteFormationProbability > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteFormationProbability = protorosetteFormationProbability;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteResolutionProbabilityPerTimestep(double protorosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (protorosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteResolutionProbabilityPerTimestep = protorosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::SetRosetteResolutionProbabilityPerTimestep(double rosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (rosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (rosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mRosetteResolutionProbabilityPerTimestep = rosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteFormationProbability()
{
    return this->mProtorosetteFormationProbability;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteResolutionProbabilityPerTimestep()
{
    return this->mProtorosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableVertexMeshWithRosettes<ELEMENT_DIM, SPACE_DIM>::GetRosetteResolutionProbabilityPerTimestep()
{
    return this->mRosetteResolutionProbabilityPerTimestep;
}

// Explicit instantiation
template class MutableVertexMeshWithRosettes<1,1>;
template class MutableVertexMeshWithRosettes<1,2>;
template class MutableVertexMeshWithRosettes<1,3>;
template class MutableVertexMeshWithRosettes<2,2>;
template class MutableVertexMeshWithRosettes<2,3>;
template class MutableVertexMeshWithRosettes<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableVertexMeshWithRosettes)