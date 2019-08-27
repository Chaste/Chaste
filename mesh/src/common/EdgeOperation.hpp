/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef EDGEOPERATION_HPP_
#define EDGEOPERATION_HPP_

#include <vector>
#include "EdgeRemapInfo.hpp"

enum EDGE_OPERATION {
    EDGE_OPERATION_ADD,
    EDGE_OPERATION_DELETE,
    EDGE_OPERATION_DIVIDE,
    EDGE_OPERATION_SPLIT
};

/**
 * Class for storing edge change during remeshing
 */
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

    /**
     * Constructor for split operation ()
     * @param operation
     * @param elementIndex
     * @param newEdges
     */
    EdgeOperation(EDGE_OPERATION operation,
                  const unsigned elementIndex,
                  EdgeRemapInfo* newEdges)
    {
        assert(operation == EDGE_OPERATION_SPLIT);
        this->mOperation = operation;
        this->mElementIndex = elementIndex;
        this->mNewEdges = newEdges;
    }

    ~EdgeOperation()
    {
        delete this->mNewEdges;
        delete this->mNewEdges2;
    }

    EDGE_OPERATION GetOperation() const
    {
        return mOperation;
    }

    unsigned GetElementIndex() const
    {
        return this->mElementIndex;
    }

    unsigned GetElementIndex2() const
    {
        return this->mElementIndex2;
    }

    unsigned GetLocalEdgeIndex() const
    {
        return mLocalEdgeIndex;
    }

    EdgeRemapInfo* GetNewEdges() const
    {
        return mNewEdges;
    }

    EdgeRemapInfo* GetNewEdges2() const
    {
        return mNewEdges2;
    }
};

#endif //EDGEOPERATION_HPP_
