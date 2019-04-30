//
// Created by twin on 25/03/19.
//

#ifndef EDGEOPERATION_HPP_
#define EDGEOPERATION_HPP_

#include <vector>
#include "EdgeRemapInfo.hpp"

enum EDGE_OPERATION {
    EDGE_OPERATION_ADD,
    EDGE_OPERATION_DELETE,
    EDGE_OPERATION_DIVIDE
};


class EdgeOperation {

    EDGE_OPERATION mOperation;


    unsigned mElementIndex;
    unsigned mElementIndex2;


    unsigned mLocalEdgeIndex;


    EdgeRemapInfo* mNewEdges;
    EdgeRemapInfo* mNewEdges2;

public:


    /**
     * Constructor for either ADD or DELETE operation
     * @param operation
     * @param elementIndex
     * @param localEdgeIndex
     */
    EdgeOperation(EDGE_OPERATION operation, unsigned elementIndex, unsigned localEdgeIndex)
    {
        assert(operation == EDGE_OPERATION_ADD || operation == EDGE_OPERATION_DELETE);

        this->mOperation = operation;
        this->mElementIndex = elementIndex;
        this->mElementIndex2 = 0;
        this->mLocalEdgeIndex = localEdgeIndex;
        this->mNewEdges = nullptr;
        this->mNewEdges2 = nullptr;
    }

     /**
      * Constructor for the DIVIDE operation
      * @param elementIndex
      * @param elementIndex2
      * @param newEdges
      * @param newEdges2
      */
    EdgeOperation(unsigned elementIndex,
                  unsigned elementIndex2,
                  EdgeRemapInfo* newEdges,
                  EdgeRemapInfo* newEdges2)
    {
        this->mOperation = EDGE_OPERATION_DIVIDE;
        this->mElementIndex = elementIndex;
        this->mElementIndex2 = elementIndex2;
        this->mNewEdges = newEdges;
        this->mNewEdges2 = newEdges2;
    }

    EDGE_OPERATION GetOperation()
    {
        return mOperation;
    }

    unsigned GetElementIndex()
    {
        return this->mElementIndex;
    }

    unsigned GetElementIndex2()
    {
        return this->mElementIndex2;
    }

    unsigned GetLocalEdgeIndex()
    {
        return mLocalEdgeIndex;
    }

    EdgeRemapInfo* GetNewEdges()
    {
        return mNewEdges;
    }

    EdgeRemapInfo* GetNewEdges2()
    {
        return mNewEdges;
    }



};

#endif //EDGEOPERATION_HPP_
