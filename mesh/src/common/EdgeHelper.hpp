//
// Created by twin on 11/02/19.
//

#ifndef EDGEHELPER_HPP_
#define EDGEHELPER_HPP_

#include <vector>
#include <map>
#include "Node.hpp"
#include "Edge.hpp"

/**
 * Class for facilitating the creation and management of unique edges in a vertex mesh
 */
 template <unsigned SPACE_DIM>
class EdgeHelper {

private:

    std::vector<Edge<SPACE_DIM>*> mEdges;
    std::map< UIndexPair, Edge<SPACE_DIM>*> mEdgesMap;

public:


    EdgeHelper()
    {

    }

    ~EdgeHelper()
    {

    }

    void Clear()
    {
        // Iterate over edges and free the memory
        for(auto edge: mEdges)
        {
            delete edge;
        }
        mEdges.clear();
    }



    Edge<SPACE_DIM>* GetEdgeFromNodes(Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1)
    {
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
        if(edgeItt == mEdgesMap.end())
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


    Edge<SPACE_DIM>* GetEdgeFromNodes(unsigned elementIndex, Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1)
    {
        auto edge = GetEdgeFromNodes(node0, node1);
        edge->AddElement(elementIndex);
        return edge;
    }


    Edge<SPACE_DIM>* GetEdge(unsigned index)
    {
        return mEdges[index];
    }

    Edge<SPACE_DIM>* GetEdge(unsigned index) const
    {
        return mEdges[index];
    }


    Edge<SPACE_DIM>* operator[](unsigned index)
    {
        return mEdges[index];
    }

    Edge<SPACE_DIM>* operator[](unsigned index) const
    {
        return mEdges[index];
    }

    /**
     * Rebuilds node-node to edge map
     */
    void UpdateEdgesMapKey()
    {
        mEdgesMap.clear();

        for(auto edge: mEdges)
        {
            mEdgesMap[edge->GetMapIndex()] = edge;
        }
    }



    unsigned GetNumEdges() const
    {
        return mEdges.size();
    }

    typename std::vector<Edge<SPACE_DIM>*>::iterator begin()
    {
        return mEdges.begin();
    }

    typename std::vector<Edge<SPACE_DIM>*>::iterator end()
    {
        return mEdges.end();
    }



};


#endif //CHASTE_EDGEHELPER_HPP
