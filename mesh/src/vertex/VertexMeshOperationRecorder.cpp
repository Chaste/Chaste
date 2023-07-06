/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "VertexMeshOperationRecorder.hpp"

#include "EdgeRemapInfo.hpp"

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::VertexMeshOperationRecorder()
    : mpEdgeHelper(nullptr)
{
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::~VertexMeshOperationRecorder()
{
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::SetEdgeHelper(EdgeHelper<SPACE_DIM>* pEdgeHelper)
{
    mpEdgeHelper = pEdgeHelper;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordT1Swap(T1SwapInfo<SPACE_DIM>& rSwapInfo)
{
    mT1Swaps.push_back(rSwapInfo);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
std::vector<T1SwapInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetT1SwapsInfo() const
{
    return mT1Swaps;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearT1SwapsInfo()
{
    mT1Swaps.clear();
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordT2Swap(T2SwapInfo<SPACE_DIM>& rSwapInfo)
{
    mT2Swaps.push_back(rSwapInfo);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
std::vector<T2SwapInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetT2SwapsInfo() const
{
    return mT2Swaps;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearT2SwapsInfo()
{
    mT2Swaps.clear();
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordT3Swap(T3SwapInfo<SPACE_DIM>& rSwapInfo)
{
    mT3Swaps.push_back(rSwapInfo);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
std::vector<T3SwapInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetT3SwapsInfo() const
{
    return mT3Swaps;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearT3SwapsInfo()
{
    mT3Swaps.clear();
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordCellDivisionInfo(CellDivisionInfo<SPACE_DIM>& rDivisionInfo)
{
    mCellDivisions.push_back(rDivisionInfo);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
std::vector<CellDivisionInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetCellDivisionInfo() const
{
    return mCellDivisions;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearCellDivisionInfo()
{
    mCellDivisions.clear();
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const std::vector<EdgeOperation>& VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetEdgeOperations()
{
    return mEdgeOperations;
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearEdgeOperations()
{
    mEdgeOperations.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordNodeMergeOperation(const std::vector<unsigned> oldIds,
                                                                                   VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                   const std::pair<unsigned, unsigned> merged_nodes_pair,
                                                                                   const bool elementIndexIsRemapped)
{
    const unsigned element_index = pElement->GetIndex();
    const unsigned element_num_edges = pElement->GetNumEdges();
    std::vector<long> edge_mapping(element_num_edges, -1);
    std::vector<unsigned> edge_status(element_num_edges, 0);

    // Marking unaffected edges
    for (unsigned i = 0; i < oldIds.size(); ++i)
    {
        long index = pElement->GetLocalEdgeIndex((*mpEdgeHelper)[oldIds[i]]);
        if (index >= 0)
        {
            edge_mapping[index] = i;
            edge_status[index] = 0;
        }
    }
    // Checking whether the deleted node is upper or lower node
    const unsigned node_A_index = merged_nodes_pair.first;
    const unsigned node_B_index = merged_nodes_pair.second;

    // Node B is also considered upper node if the last two nodes are merged
    bool is_B_upper = node_B_index > node_A_index;
    if (node_A_index == 0)
    {
        is_B_upper = node_B_index == 1;
    }
    if (node_B_index == 0)
    {
        // This line is excluded from coverage, as it is very difficult to test this case
        // This case becomes relevant in long time simulations of proliferating tissue
        // and difficult to reproduce
        is_B_upper = node_A_index != 1; // LCOV_EXCL_LINE
    }

    unsigned lower_node = node_A_index;
    unsigned upper_node = node_B_index;
    if (!is_B_upper)
    {
        lower_node = node_B_index;
        upper_node = node_A_index;
    }
    // Previous edge denotes the edge below the lower node index
    // and next_edge denotes the edge above the upper node index
    unsigned prev_edge = 0;
    if (is_B_upper)
    {
        if (upper_node == 0)
        {
            // This line is excluded from coverage, as it is very difficult to test this case
            // This case becomes relevant in long time simulations of proliferating tissue
            // and difficult to reproduce
            prev_edge = element_num_edges - 2; // LCOV_EXCL_LINE
        }
        else if (upper_node == 1)
        {
            prev_edge = element_num_edges - 1;
        }
        else
        {
            prev_edge = upper_node - 2;
        }
    }
    else
    {
        if (upper_node == 0 || upper_node == 1)
        {
            prev_edge = element_num_edges - 1;
        }
        else
        {
            prev_edge = upper_node - 2;
        }
    }
    const unsigned next_edge = (prev_edge + 1) % element_num_edges;

    // Marking edges below and above the deleted edge
    edge_status[prev_edge] = 3;
    edge_status[next_edge] = 3;

    // The edge below node A is in the old element if node B is upper node, and marked with status 0 in the loop above
    // Because node deletion during node merging also removes the pair of edges associated to it,
    // we need to make sure the edges associated with nodes about to be merges correctly map to the old element edges
    if (is_B_upper)
    {
        edge_mapping[next_edge] = node_B_index;
    }
    else
    {
        if (lower_node == 0)
        {
            // This line is excluded from coverage, as it is very difficult to test this case
            // This case becomes relevant in long time simulations of proliferating tissue
            // and difficult to reproduce
            edge_mapping[prev_edge] = element_num_edges; // LCOV_EXCL_LINE
        }
        else
        {
            edge_mapping[prev_edge] = lower_node - 1;
        }
    }
    // Sanity check
    for (unsigned i = 0; i < edge_mapping.size(); ++i)
    {
        assert(edge_mapping[i] >= 0);
    }

    const EdgeRemapInfo remap_info(edge_mapping, edge_status);
    mEdgeOperations.emplace_back(EDGE_OPERATION_NODE_MERGE, element_index, remap_info, elementIndexIsRemapped);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordEdgeSplitOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                   const unsigned edgeIndex,
                                                                                   const double insertedNodeRelPosition,
                                                                                   const bool elementIndexIsRemapped)
{
    const unsigned element_index = pElement->GetIndex();
    const unsigned element_num_edges = pElement->GetNumEdges();
    std::vector<double> thetas(element_num_edges);
    std::vector<long> edge_mapping(element_num_edges);
    std::vector<unsigned> edge_status(element_num_edges, 0);

    // Daughter edge indices
    const unsigned split_1 = edgeIndex;
    const unsigned split_2 = edgeIndex + 1;
    edge_status[split_1] = 1;
    edge_status[split_2] = 1;
    thetas[split_1] = insertedNodeRelPosition;
    thetas[split_2] = 1.0 - insertedNodeRelPosition;
    unsigned count = 0;
    for (unsigned i = 0; i < element_num_edges; ++i)
    {
        edge_mapping[i] = i - count;
        if (edge_status[i] == 1)
        {
            count = 1;
        }
    }

    EdgeRemapInfo remap_info(edge_mapping, edge_status);
    remap_info.SetSplitProportions(thetas);
    mEdgeOperations.emplace_back(EDGE_OPERATION_SPLIT, element_index, remap_info, elementIndexIsRemapped);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordCellDivideOperation(const std::vector<unsigned>& rOldIds,
                                                                                    VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement1,
                                                                                    VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement2)
{
    const unsigned num_edges_1 = pElement1->GetNumEdges();
    const unsigned num_edges_2 = pElement2->GetNumEdges();
    std::vector<long> edge_mapping_1(num_edges_1, -2);
    std::vector<long> edge_mapping_2(num_edges_2, -2);
    std::vector<unsigned> edge_status_1(num_edges_1);
    std::vector<unsigned> edge_status_2(num_edges_2);

    std::vector<unsigned> old_split_edges(rOldIds.size());
    // Keeps track of parent edges that are NOT retained in daughter cells
    for (unsigned i = 0; i < rOldIds.size(); ++i)
    {
        old_split_edges[i] = i;
    }
    unsigned counter_1 = 0;
    unsigned counter_2 = 0;

    // First find parent edges that correspond directly to daughter cells' edges
    // At the end of the loop, old_split_edges contains parent edge indices that are split
    for (unsigned i = 0; i < rOldIds.size(); ++i)
    {
        // Index of parent edge corresponding to daughter cell's edge.
        //-1 if not found.
        long index_1 = pElement1->GetLocalEdgeIndex((*mpEdgeHelper)[rOldIds[i]]);
        long index_2 = pElement2->GetLocalEdgeIndex((*mpEdgeHelper)[rOldIds[i]]);
        auto position = std::find(old_split_edges.begin(), old_split_edges.end(), i);

        // Modify edge map and status
        if (index_1 >= 0)
        {
            edge_mapping_1[index_1] = i;
            edge_status_1[index_1] = 0;
            old_split_edges.erase(position);
            counter_1++;
        }
        if (index_2 >= 0)
        {
            edge_mapping_2[index_2] = i;
            edge_status_2[index_2] = 0;
            old_split_edges.erase(position);
            counter_2++;
        }
    }
    // Two parent edges are split
    assert(old_split_edges.size() == 2);

    // Three edges in daughter cells are unmapped
    assert(counter_1 == num_edges_1 - 3);
    assert(counter_2 == num_edges_2 - 3);

    // Edge split proportions.
    std::vector<double> thetas_1(num_edges_1);
    std::vector<double> thetas_2(num_edges_2);

    // Go through unmapped edges of daughter cell to find a mapping between parent split edge and
    // daughter edge
    std::vector<unsigned> old_split_edges_1(old_split_edges);
    for (unsigned i = 0; i < num_edges_1; ++i)
    {
        if (edge_mapping_1[i] == -2)
        {
            auto p_node_1 = pElement1->GetEdge(i)->GetNode(0);
            auto p_node_2 = pElement1->GetEdge(i)->GetNode(1);
            bool split_edge_found = false;
            for (unsigned j = 0; j < old_split_edges_1.size(); ++j)
            {
                auto old_edge = (*mpEdgeHelper)[rOldIds[old_split_edges_1[j]]];
                if (old_edge->ContainsNode(p_node_1) || old_edge->ContainsNode(p_node_2))
                {
                    edge_mapping_1[i] = old_split_edges_1[j];
                    edge_status_1[i] = 1;
                    counter_1++;
                    split_edge_found = true;
                    thetas_1[i] = pElement1->GetEdge(i)->rGetLength() / old_edge->rGetLength();
                    old_split_edges_1.erase(old_split_edges_1.begin() + j);
                    break;
                }
            }
            if (!split_edge_found)
            {
                edge_mapping_1[i] = -1;
                edge_status_1[i] = 2;
                counter_1++;
            }
        }
    }

    for (unsigned i = 0; i < num_edges_2; ++i)
    {
        if (edge_mapping_2[i] == -2)
        {
            auto p_node_1 = pElement2->GetEdge(i)->GetNode(0);
            auto p_node_2 = pElement2->GetEdge(i)->GetNode(1);
            bool split_edge_found = false;
            for (unsigned j = 0; j < old_split_edges.size(); ++j)
            {
                auto old_edge = (*mpEdgeHelper)[rOldIds[old_split_edges[j]]];
                if (old_edge->ContainsNode(p_node_1) || old_edge->ContainsNode(p_node_2))
                {
                    edge_mapping_2[i] = old_split_edges[j];
                    edge_status_2[i] = 1;
                    counter_2++;
                    split_edge_found = true;
                    thetas_2[i] = pElement2->GetEdge(i)->rGetLength() / old_edge->rGetLength();
                    old_split_edges.erase(old_split_edges.begin() + j);
                    break;
                }
            }
            if (!split_edge_found)
            {
                edge_mapping_2[i] = -1;
                edge_status_2[i] = 2;
                counter_2++;
            }
        }
    }
    assert(old_split_edges_1.empty());
    assert(old_split_edges.empty());

    // Checking if all edges of daughter cells have been mapped.
    assert(counter_1 == num_edges_1);
    assert(counter_2 == num_edges_2);

    EdgeRemapInfo remap_info_1(edge_mapping_1, edge_status_1);
    EdgeRemapInfo remap_info_2(edge_mapping_2, edge_status_2);
    remap_info_1.SetSplitProportions(thetas_1);
    remap_info_2.SetSplitProportions(thetas_2);
    mEdgeOperations.emplace_back(pElement1->GetIndex(), pElement2->GetIndex(), remap_info_1, remap_info_2);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordNewEdgeOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                 const unsigned edgeIndex)
{
    const unsigned element_index = pElement->GetIndex();
    const unsigned num_edges = pElement->GetNumEdges();
    std::vector<long> edge_mapping(num_edges, 0);
    std::vector<unsigned> edge_status(num_edges);
    for (unsigned i = 0; i < edgeIndex; ++i)
    {
        edge_mapping[i] = i;
        edge_status[i] = 0;
    }
    edge_mapping[edgeIndex] = -1;
    edge_status[edgeIndex] = 2;
    for (unsigned i = edgeIndex + 1; i < num_edges; ++i)
    {
        edge_mapping[i] = i - 1;
        edge_status[i] = 0;
    }

    const EdgeRemapInfo remap_info(edge_mapping, edge_status);
    mEdgeOperations.emplace_back(EDGE_OPERATION_ADD, element_index, remap_info);
}

template<unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordEdgeMergeOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                   const unsigned nodeIndex)
{
    const unsigned element_index = pElement->GetIndex();
    const unsigned num_edges = pElement->GetNumEdges();
    std::vector<long> edge_mapping(num_edges, 0);
    std::vector<unsigned> edge_status(num_edges, 0);

    // Here we find the edge with the lower index.
    // High index edge is merged into low edge index
    unsigned low_edge = (nodeIndex + num_edges) % (num_edges + 1);
    unsigned high_edge = nodeIndex;

    // If the first edge was merged into the last edge
    if (low_edge > high_edge)
    {
        edge_status[low_edge - 1] = 4;
        for (unsigned i = 0; i < num_edges; ++i)
        {
            edge_mapping[i] = i + 1;
        }
        edge_mapping[low_edge - 1] = low_edge;
    }
    else
    {
        edge_status[low_edge] = 4;
        for (unsigned i = 0; i < high_edge; ++i)
        {
            edge_mapping[i] = i;
        }
        for (unsigned i = high_edge; i < num_edges; ++i)
        {
            edge_mapping[i] = i + 1;
        }
    }

    const EdgeRemapInfo remap_info(edge_mapping, edge_status);
    mEdgeOperations.emplace_back(EDGE_OPERATION_MERGE, element_index, remap_info);
}

template class VertexMeshOperationRecorder<1, 1>;
template class VertexMeshOperationRecorder<1, 2>;
template class VertexMeshOperationRecorder<1, 3>;
template class VertexMeshOperationRecorder<2, 2>;
template class VertexMeshOperationRecorder<2, 3>;
template class VertexMeshOperationRecorder<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMeshOperationRecorder)