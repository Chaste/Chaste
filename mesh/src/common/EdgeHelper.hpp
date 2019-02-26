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


    EdgeHelper();

    ~EdgeHelper();

    void Clear();

    Edge<SPACE_DIM>* GetEdgeFromNodes(Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1);
    Edge<SPACE_DIM>* GetEdgeFromNodes(unsigned elementIndex, Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1);

    Edge<SPACE_DIM>* GetEdge(unsigned index);
    Edge<SPACE_DIM>* GetEdge(unsigned index) const;


    Edge<SPACE_DIM>* operator[](unsigned index);
    Edge<SPACE_DIM>* operator[](unsigned index) const;

    void RemoveDeletedEdges();

    /**
     * Rebuilds node-node to edge map
     */
    void UpdateEdgesMapKey();

    unsigned GetNumEdges() const;

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
