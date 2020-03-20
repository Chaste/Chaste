#include "VertexMeshOperationRecorder.hpp"
#include "EdgeRemapInfo.hpp"
template<unsigned ELEMENT_DIM, unsigned int SPACE_DIM>
VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::VertexMeshOperationRecorder()
:mpEdgeHelper(nullptr)
{}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::~VertexMeshOperationRecorder()
{
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::SetEdgeHelper(EdgeHelper<SPACE_DIM> *pEdgeHelper)
{
    mpEdgeHelper = pEdgeHelper;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordT1Swap(T1SwapInfo<SPACE_DIM>& rSwap_info)
{
    mT1Swaps.push_back(rSwap_info);
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
std::vector<T1SwapInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetT1SwapsInfo() const
{
    return mT1Swaps;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearT1SwapsInfo()
{
    mT1Swaps.clear();
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordT2Swap(T2SwapInfo<SPACE_DIM>& rSwap_info)
{
    mT2Swaps.push_back(rSwap_info);
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
std::vector<T2SwapInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetT2SwapsInfo() const
{
    return mT2Swaps;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearT2SwapsInfo()
{
    mT2Swaps.clear();
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordT3Swap(T3SwapInfo<SPACE_DIM>& rSwap_info)
{
    mT3Swaps.push_back(rSwap_info);
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
std::vector<T3SwapInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetT3SwapsInfo() const
{
    return mT3Swaps;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearT3SwapsInfo()
{
    mT3Swaps.clear();
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordCellDivisionInfo(CellDivisionInfo<SPACE_DIM>& rDivision_info)
{
    mCellDivisions.push_back(rDivision_info);
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
std::vector<CellDivisionInfo<SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetCellDivisionInfo() const
{
    return mCellDivisions;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearCellDivisionInfo()
{
    mCellDivisions.clear();
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
const std::vector<EdgeOperation*> & VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetEdgeOperations()
{
    return mEdgeOperations;
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearEdgeOperations()
{
    for (auto operation : mEdgeOperations)
    {
        delete operation;
    }
    mEdgeOperations.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordNodeMergeOperation(const std::vector<unsigned int> oldIds,
                                                                                   VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                                                                   const std::pair<unsigned int, unsigned int> merged_nodes_pair,
                                                                                   const bool elementIndexIsRemapped)
{
    const unsigned int element_index = pElement->GetIndex();
    const unsigned int elementNumEdges = pElement->GetNumEdges();
    std::vector<long int> edge_mapping(elementNumEdges,-1);
    std::vector<unsigned int> edge_status(elementNumEdges,0);

    //Marking unaffected edges
    for (unsigned int i=0; i<oldIds.size(); ++i)
    {
        long index = pElement->GetLocalEdgeIndex((*mpEdgeHelper)[oldIds[i]]);
        if (index>=0)
        {
            edge_mapping[index] = i;
            edge_status[index] = 0;
        }
    }
    //Checking whether the deleted node is upper or lower node
    const unsigned int node_A_index = merged_nodes_pair.first;
    const unsigned int node_B_index = merged_nodes_pair.second;
    //Node B is also considered upper node if the last two nodes are merged
    const bool is_B_upper = node_B_index>node_A_index||(node_B_index==0&&node_A_index==oldIds.size()-1);
    unsigned int lower_node = node_A_index;
    unsigned int upper_node = node_B_index;
    if (!is_B_upper)
    {
        lower_node = node_B_index;
        upper_node = node_A_index;
    }
    //Previous edge denotes the edge below the lower node index
    //and nextEdge denotes the edge above the upper node index
    unsigned int prevEdge = 0;
    if (upper_node == 0)
    {
        prevEdge = elementNumEdges - 2;
    }
    else if (upper_node == 1)
    {
        prevEdge = elementNumEdges - 1;
    }
    else
    {
        prevEdge = upper_node - 2;
    }
    const unsigned int nextEdge = (prevEdge+1)%elementNumEdges;

    //Marking edges below and above the deleted edge
    edge_status[prevEdge] = 3;
    edge_status[nextEdge] = 3;
    //The edge below node A is in the old element if node B is upper node, and marked with status 0 in the loop above
    //Because node deletion during node merging also removes the pair of edges associated to it,
    //we need to make sure the edges associated with nodes about to be merges correctly map to the old element edges
    if (is_B_upper)
    {
        edge_mapping[nextEdge] = node_B_index;
    }
    else
    {
        if (lower_node==0)
        {
            edge_mapping[prevEdge] = elementNumEdges;
        }
        else
        {
            edge_mapping[prevEdge] = lower_node-1;
        }
    }
    //Sanity check
    for (unsigned int i=0; i<edge_mapping.size(); ++i)
    {
        assert(edge_mapping[i]>=0);
    }

    EdgeRemapInfo* remap_info = new EdgeRemapInfo(edge_mapping, edge_status);
    mEdgeOperations.push_back(new EdgeOperation(EDGE_OPERATION_NODE_MERGE, element_index, remap_info,elementIndexIsRemapped));
}

template <unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordEdgeSplitOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                   const unsigned int edge_index,
                                                                                   const double inserted_node_rel_position,
                                                                                   const bool elementIndexIsRemapped)
{
    const unsigned int element_index = pElement->GetIndex();
    const unsigned int elementNumEdges = pElement->GetNumEdges();
    std::vector<double> thetas(elementNumEdges);
    std::vector<long int> edge_mapping(elementNumEdges);
    std::vector<unsigned int> edge_status(elementNumEdges,0);
    //Daughter edge indices
    const unsigned int split_1 = edge_index;
    const unsigned int split_2 = edge_index+1;
    edge_status[split_1] = 1;
    edge_status[split_2] = 1;
    thetas[split_1] = inserted_node_rel_position;
    thetas[split_2] = 1.0-inserted_node_rel_position;
    unsigned int count = 0;
    for (unsigned int i=0; i < elementNumEdges; ++i)
    {
        edge_mapping[i] = i - count;
        if (edge_status[i]==1)
            count = 1;
    }
    EdgeRemapInfo* newEdges = new EdgeRemapInfo(edge_mapping, edge_status);
    newEdges->SetSplitProportions(thetas);
    mEdgeOperations.push_back(new EdgeOperation(EDGE_OPERATION_SPLIT, element_index, newEdges, elementIndexIsRemapped));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordCellDivideOperation(const std::vector<unsigned int>& oldIds,
                                                                                    VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement1,
                                                                                    VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement2)
{
    const unsigned int n_edges_1 = pElement1->GetNumEdges();
    const unsigned int n_edges_2 = pElement2->GetNumEdges();
    std::vector<long> edge_mapping_1(n_edges_1,-2);
    std::vector<long> edge_mapping_2(n_edges_2,-2);
    std::vector<unsigned int> edge_status_1(n_edges_1);
    std::vector<unsigned int> edge_status_2(n_edges_2);

    std::vector<unsigned int> old_split_edges(oldIds.size());
    //Keeps track of parent edges that are NOT retained in daughter cells
    for (unsigned int i=0; i<oldIds.size(); ++i)
        old_split_edges[i] = i;
    unsigned int counter_1 = 0;
    unsigned int counter_2 = 0;
    //First find parent edges that correspond directly to daughter cells' edges
    //At the end of the loop, old_split_edges contains parent edge indices that are split
    for (unsigned int i=0; i<oldIds.size(); ++i)
    {
        //Index of parent edge corresponding to daughter cell's edge.
        //-1 if not found.
        long index_1 = pElement1->GetLocalEdgeIndex((*mpEdgeHelper)[oldIds[i]]);
        long index_2 = pElement2->GetLocalEdgeIndex((*mpEdgeHelper)[oldIds[i]]);
        auto position = std::find(old_split_edges.begin(), old_split_edges.end(),i);
        //Modify edge map and status
        if (index_1>=0)
        {
            edge_mapping_1[index_1] = i;
            edge_status_1[index_1] = 0;
            old_split_edges.erase(position);
            counter_1++;
        }
        if (index_2>=0)
        {
            edge_mapping_2[index_2] = i;
            edge_status_2[index_2] = 0;
            old_split_edges.erase(position);
            counter_2++;
        }
    }
    if (old_split_edges.size()!=2)
        EXCEPTION("edge split size is wrong");
    //Two parent edges are split
    assert(old_split_edges.size()==2);
    //Three edges in daughter cells are unmapped
    assert(counter_1==n_edges_1-3);
    assert(counter_2==n_edges_2-3);
    //Edge split proportions.
    std::vector<double> thetas_1(n_edges_1);
    std::vector<double> thetas_2(n_edges_2);
    //Go through unmapped edges of daughter cell to find a mapping between parent split edge and
    //daughter edge
    std::vector<unsigned int> old_split_edges_1(old_split_edges);
    for (unsigned int i=0; i<n_edges_1; ++i)
    {
        if (edge_mapping_1[i]==-2)
        {
            auto node_1 = pElement1->GetEdge(i)->GetNode(0);
            auto node_2 = pElement1->GetEdge(i)->GetNode(1);
            bool split_edge_found = false;
            for (unsigned int j=0; j<old_split_edges_1.size(); ++j)
            {
                auto old_edge = (*mpEdgeHelper)[oldIds[old_split_edges_1[j]]];
                if (old_edge->ContainsNode(node_1)||old_edge->ContainsNode(node_2))
                {
                    edge_mapping_1[i] = old_split_edges_1[j];
                    edge_status_1[i] = 1;
                    counter_1++;
                    split_edge_found = true;
                    thetas_1[i] = pElement1->GetEdge(i)->rGetLength()/old_edge->rGetLength();
                    old_split_edges_1.erase(old_split_edges_1.begin()+j);
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

    for (unsigned int i=0; i<n_edges_2; ++i)
    {
        if (edge_mapping_2[i]==-2)
        {
            auto node_1 = pElement2->GetEdge(i)->GetNode(0);
            auto node_2 = pElement2->GetEdge(i)->GetNode(1);
            bool split_edge_found = false;
            for (unsigned int j=0; j<old_split_edges.size(); ++j)
            {
                auto old_edge = (*mpEdgeHelper)[oldIds[old_split_edges[j]]];
                if (old_edge->ContainsNode(node_1)||old_edge->ContainsNode(node_2))
                {
                    edge_mapping_2[i] = old_split_edges[j];
                    edge_status_2[i] = 1;
                    counter_2++;
                    split_edge_found = true;
                    thetas_2[i] = pElement2->GetEdge(i)->rGetLength()/old_edge->rGetLength();
                    old_split_edges.erase(old_split_edges.begin()+j);
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
    //Checking if all edges of daughter cells have been mapped.
    assert(counter_1==n_edges_1);
    assert(counter_2==n_edges_2);

    EdgeRemapInfo* remap_info_1 = new EdgeRemapInfo(edge_mapping_1, edge_status_1);
    EdgeRemapInfo* remap_info_2 = new EdgeRemapInfo(edge_mapping_2, edge_status_2);
    remap_info_1->SetSplitProportions(thetas_1);
    remap_info_2->SetSplitProportions(thetas_2);
    mEdgeOperations.push_back(new EdgeOperation(pElement1->GetIndex(), pElement2->GetIndex(), remap_info_1, remap_info_2));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordNewEdgeOperation(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                                                                 const unsigned int edge_index)
{
    const unsigned int element_index = pElement->GetIndex();
    const unsigned int n_edges = pElement->GetNumEdges();
    std::vector<long> edge_mapping(n_edges,0);
    std::vector<unsigned int> edge_status(n_edges);
    for (unsigned int i=0; i<edge_index; ++i)
    {
        edge_mapping[i] = i;
        edge_status[i] = 0;
    }
    edge_mapping[edge_index] = -1;
    edge_status[edge_index] = 2;
    for (unsigned int i=edge_index+1; i<n_edges; ++i)
    {
        edge_mapping[i] = i-1;
        edge_status[i] = 0;
    }
    EdgeRemapInfo* remap_info = new EdgeRemapInfo(edge_mapping, edge_status);
    mEdgeOperations.push_back(new EdgeOperation(EDGE_OPERATION_ADD, element_index, remap_info));
}

template <unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::RecordEdgeMergeOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                         const unsigned int node_index)
{
    const unsigned int element_index = pElement->GetIndex();
    const unsigned int n_edges = pElement->GetNumEdges();
    std::vector<long> edge_mapping(n_edges,0);
    std::vector<unsigned int> edge_status(n_edges,0);

    //Here we find the edge with the lower index.
    //High index edge is merged into low edge index

    unsigned int low_edge = (node_index+n_edges)%(n_edges+1);
    unsigned int high_edge = node_index;

    //If the first edge was merged into the last edge
    if (low_edge>high_edge)
    {
        edge_status[low_edge-1] = 4;
        for (unsigned int i=0; i<n_edges; ++i)
        {
            edge_mapping[i] = i+1;
        }
        edge_mapping[low_edge-1] = low_edge;
    }
    else
    {
        edge_status[low_edge] = 4;
        for (unsigned int i=0; i<high_edge; ++i)
        {
            edge_mapping[i]=i;
        }
        for (unsigned int i=high_edge; i<n_edges; ++i)
        {
            edge_mapping[i] = i+1;
        }
    }

    EdgeRemapInfo* remap_info = new EdgeRemapInfo(edge_mapping, edge_status);
    mEdgeOperations.push_back(new EdgeOperation(EDGE_OPERATION_MERGE, element_index, remap_info));
}

template class VertexMeshOperationRecorder<1,1>;
template class VertexMeshOperationRecorder<1,2>;
template class VertexMeshOperationRecorder<1,3>;
template class VertexMeshOperationRecorder<2,2>;
template class VertexMeshOperationRecorder<2,3>;
template class VertexMeshOperationRecorder<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMeshOperationRecorder)

