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

#ifndef EDGEOPERATION_HPP_
#define EDGEOPERATION_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include "ChasteSerialization.hpp"
#include "EdgeRemapInfo.hpp"

/**
 * Possible types of edge operation.
 */
enum EDGE_OPERATION {
    EDGE_OPERATION_ADD,
    EDGE_OPERATION_DIVIDE,
    EDGE_OPERATION_SPLIT,
    EDGE_OPERATION_NODE_MERGE,
    EDGE_OPERATION_MERGE
};

/**
 * Class for storing edge operation during remeshing.
 */
class EdgeOperation
{
private:

    /** Type of edge operation. */
    EDGE_OPERATION mOperation;

    /** Index of one element sharing the edge. */
    unsigned mElementIndex;

    /** Index of another element sharing the edge. */
    unsigned mElementIndex2;

    /** Edge remap info. */
    EdgeRemapInfo mRemapInfo;

    /** Edge remap info after cell division. */    
    EdgeRemapInfo mRemapInfo2;

    /**
     * If operation is recorded before element indices are changed. For example, 
     * if the operations recorded during T2 swap.
     */
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
        archive & mRemapInfo;
        archive & mRemapInfo2;
    }
  
public:

    /**
     * Default constructor.
     */
    EdgeOperation();

    /**
     * Constructor for add, split, node and edge merge operations.
     * 
     * @param operation which operation was performed
     * @param elementIndex index of the element
     * @param remapInfo remapping info
     * @param isIndexRemapped indicates whether the operation has been recorded 
     *                        before elements changed their indices (defaults to 
     *                        false); e.g. when T2 swap occurs, node merging 
     *                        operation is recorded with element index that will 
     *                        be modified in RemoveDeletedNodesAndELements(). 
     *                        See also VertexBasedCellPopulation::Update().
     */
    EdgeOperation(EDGE_OPERATION operation, 
                  unsigned elementIndex,
                  EdgeRemapInfo remapInfo, 
                  const bool isIndexRemapped = false);

     /**
      * Constructor for the DIVIDE operation.
      * 
      * @param elementIndex an element index
      * @param elementIndex2 another element index
      * @param pRemapInfo edge remap info
      * @param pRemapInfo2 edge remap info after cell division
      */
    EdgeOperation(unsigned elementIndex,
                  unsigned elementIndex2,
                  EdgeRemapInfo remapInfo,
                  EdgeRemapInfo remapInfo2);

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
     * Modify element index on which edge operation has been performed.
     * 
     * @paran index an element index 
     */
    void SetElementIndex(const unsigned index);

    /**
     * @return Element index of the second daughter cell
     */
    unsigned GetElementIndex2() const;

    /**
     * Modify element index of the second daughter cell.
     * 
     * @paran index an element index 
     */
    void SetElementIndex2(const unsigned index);

    /**
     * @return Edge remap info
     */
    const EdgeRemapInfo& rGetRemapInfo() const;

    /**
     * @return Edge remap info after cell division
     */
    const EdgeRemapInfo& rGetRemapInfo2() const;

    /**
     * @return mIsElementIndexRemapped
     */
    bool IsElementIndexRemapped() const;
};

#endif /* EDGEOPERATION_HPP_ */