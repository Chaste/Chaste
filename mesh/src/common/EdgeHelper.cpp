//
// Created by twin on 11/02/19.
//

#include "EdgeHelper.hpp"


template<unsigned int SPACE_DIM>
EdgeHelper<SPACE_DIM>::EdgeHelper(): holdEdgeOperations(false) {

}

template<unsigned int SPACE_DIM>
EdgeHelper<SPACE_DIM>::~EdgeHelper() {

}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::Clear() {
    // Iterate over edges and free the memory
    for(auto edge: mEdges)
    {
        delete edge;
    }
    mEdges.clear();
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::GetEdgeFromNodes(Node<SPACE_DIM> *node0, Node<SPACE_DIM> *node1) {
    assert(node0->GetIndex() != node1->GetIndex());

    //Swap node so we always have the lower index as node 0
    if(node0->GetIndex() > node1->GetIndex())
    {
        auto swapNode = node0;
        node0 = node1;
        node1 = swapNode;
    }

    auto edgeMapIndices = UIndexPair(node0->GetIndex(), node1->GetIndex());

    //Check that an edge hasn't been created already
    Edge<SPACE_DIM>* edge = nullptr;

    auto edgeItt = mEdgesMap.find(edgeMapIndices);
    if(edgeItt == mEdgesMap.end() || edgeItt->second->IsDeleted())
    {
        edge = new Edge<SPACE_DIM>(mEdges.size(), node0, node1);
        mEdgesMap[edgeMapIndices] = edge;
        mEdges.push_back(edge);
    }
    else
    {
        edge = edgeItt->second;
    }


    return edge;
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *
EdgeHelper<SPACE_DIM>::GetEdgeFromNodes(unsigned elementIndex, Node<SPACE_DIM> *node0, Node<SPACE_DIM> *node1) {
    auto edge = GetEdgeFromNodes(node0, node1);
    edge->AddElement(elementIndex);
    return edge;
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::GetEdge(unsigned index) {
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::GetEdge(unsigned index) const {
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::operator[](unsigned index) {
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
Edge<SPACE_DIM> *EdgeHelper<SPACE_DIM>::operator[](unsigned index) const {
    return mEdges[index];
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::RemoveDeletedEdges() {
    // Remove any nodes that have been marked for deletion and store all other nodes in a temporary structure
    std::vector<Edge<SPACE_DIM>*> live_edges;
    for (unsigned i=0; i<this->mEdges.size(); i++)
    {
        if (this->mEdges[i]->GetNumElements() == 0)
        {
            delete this->mEdges[i];
        }
        else
        {
            live_edges.push_back(this->mEdges[i]);
        }
    }

    // Repopulate the nodes vector and reset the list of deleted node indices
    this->mEdges = live_edges;

    // Finally, reset the node indices to run from zero
    for (unsigned i=0; i<this->mEdges.size(); i++)
    {
        this->mEdges[i]->SetIndex(i);
    }
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::UpdateEdgesMapKey() {
    mEdgesMap.clear();

    for(auto edge: mEdges)
    {
        mEdgesMap[edge->GetMapIndex()] = edge;
    }
}

template<unsigned int SPACE_DIM>
unsigned EdgeHelper<SPACE_DIM>::GetNumEdges() const {
    return mEdges.size();
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::InsertAddEdgeOperation(unsigned elementIndex, unsigned localEdgeIndex) {
    if(!holdEdgeOperations)
        mEdgeOperations.push_back(EdgeOperation(EDGE_OPERATION_ADD, elementIndex, localEdgeIndex));
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::InsertDeleteEdgeOperation(unsigned elementIndex,
                                                      unsigned localEdgeIndex) {
    if(!holdEdgeOperations)
        mEdgeOperations.push_back(EdgeOperation(EDGE_OPERATION_DELETE, elementIndex, localEdgeIndex));
}

template<unsigned int SPACE_DIM>
void EdgeHelper<SPACE_DIM>::InsertCellDivideOperation(unsigned elementIndex,
                                                      unsigned elementIndex2,
                                                      std::vector<long int> newEdges, std::vector<long int> newEdges2) {
    if(!holdEdgeOperations)
        mEdgeOperations.push_back(EdgeOperation(elementIndex, elementIndex2, newEdges, newEdges2));
}

template class EdgeHelper<1>;
template class EdgeHelper<2>;
template class EdgeHelper<3>;