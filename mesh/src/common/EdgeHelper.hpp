//
// Created by twin on 11/02/19.
//

#ifndef EDGEHELPER_HPP_
#define EDGEHELPER_HPP_

#include <vector>
#include "Node.hpp"
#include "Edge.hpp"

/**
 * Class for facilitating the creation and management of unique edges in a vertex mesh
 */
 template <unsigned SPACE_DIM>
class EdgeHelper {

private:

    std::vector<Edge<SPACE_DIM>*>* mEdges;

public:

    EdgeHelper(): mEdges(nullptr)
    {

    }

    void SetEdges(std::vector<Edge<SPACE_DIM>*>* edgesVector)
    {
        mEdges = edgesVector;
    }

    Edge<SPACE_DIM>* GetEdgeFromNodes(Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1)
    {


        return nullptr;
    }


    Edge<SPACE_DIM>* GetEdgeFromNodes(unsigned elementIndex, Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1)
    {


        return nullptr;
    }



};


#endif //CHASTE_EDGEHELPER_HPP
