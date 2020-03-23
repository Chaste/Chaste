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
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

enum EDGE_OPERATION {
    EDGE_OPERATION_ADD,
    EDGE_OPERATION_DIVIDE,
    EDGE_OPERATION_SPLIT,
    EDGE_OPERATION_NODE_MERGE,
    EDGE_OPERATION_MERGE
};

/**
 * Class for storing edge operation during remeshing
 */
class EdgeOperation {
private:
    EDGE_OPERATION mOperation;

    unsigned mElementIndex;
    unsigned mElementIndex2;

    EdgeRemapInfo* mpRemapInfo;
    EdgeRemapInfo* mpRemapInfo2;

    bool mIsElementIndexRemapped;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mOperation;
        archive & mElementIndex;
        archive & mElementIndex2;
        archive & mIsElementIndexRemapped;
        archive & mpRemapInfo;
        archive & mpRemapInfo2;
    }
public:
    /**
     * Default constructor. Here for boost serialization purposes
     */
    EdgeOperation();

    /**
     * Constructor for add, split, node and edge merge operations
     * @param operation - which operation was performed
     * @param elementIndex - index of the element
     * @param pRemapInfo - remapping info
     * @param isIndexRemapped - indicates whether the operation has been recorded before elements changed their indices
     * e.g. when T2 swap occurs, node merging operation is recorded with element index that will be modified in
     * RemoveDeletedNodesAndELements() function. See also Update() function in VertexBasedCellPopulation class
     */
    EdgeOperation(EDGE_OPERATION operation, unsigned elementIndex,
                  EdgeRemapInfo* pRemapInfo, const bool isIndexRemapped = false);

     /**
      * Constructor for the DIVIDE operation
      * @param elementIndex
      * @param elementIndex2
      * @param pRemapInfo
      * @param pRemapInfo2
      */
    EdgeOperation(unsigned elementIndex,
                  unsigned elementIndex2,
                  EdgeRemapInfo* pRemapInfo,
                  EdgeRemapInfo* pRemapInfo2);
    /**
     * Destructor
     */
    ~EdgeOperation();

    /**
     * @return edge operations
     */
    EDGE_OPERATION GetOperation() const;

    /**
     * @return Element index on which edge operation has been performed.
     * Also index (inherited from mother cell) of the first daughter cell.
     */
    unsigned GetElementIndex() const;

    /**
     * Modify element index on which edge operation has been performed
     */
    void SetElementIndex(const unsigned int index);

    /**
     * @return Element index of the second daughter cell
     */
    unsigned GetElementIndex2() const;

    /**
     * Modify element index of the second daughter cell
     */
    void SetElementIndex2(const unsigned int index);

    /**
     * @return Edge remap info
     */
    EdgeRemapInfo* GetRemapInfo() const;

    /**
     * @return Edge remap info after cell division
     */
    EdgeRemapInfo* GetRemapInfo2() const;

    /**
     * @return mIsElementIndexRemapped
     */
    bool IsElementIndexRemapped() const;
};

#endif //EDGEOPERATION_HPP_
