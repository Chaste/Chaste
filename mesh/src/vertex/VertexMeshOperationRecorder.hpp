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

#ifndef VERTEXMESHOPERATIONRECORDER_HPP_
#define VERTEXMESHOPERATIONRECORDER_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "CellDivisionInfo.hpp"
#include "ChasteSerialization.hpp"
#include "EdgeHelper.hpp"
#include "EdgeOperation.hpp"
#include "T1SwapInfo.hpp"
#include "T2SwapInfo.hpp"
#include "T3SwapInfo.hpp"
#include "VertexElement.hpp"

/**
 * This class records operations performed on the mesh. In particular, this 
 * class records operations on edges and nodes during e.g. T1 transition or 
 * cell division.
 *
 * The sequence of operations as well as their nature are needed for remapping 
 * of old (prior to an operation) edge based quantities into new state. For 
 * example, when an edge is split into two or, shrinks, the edge quantities may 
 * change according to the kind of operation.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshOperationRecorder
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & mT1Swaps;
        archive & mT2Swaps;
        archive & mT3Swaps;
        archive & mCellDivisions;
        archive & mEdgeOperations;
    }

    /** Storage for T1 swap info. */
    std::vector<T1SwapInfo<SPACE_DIM> > mT1Swaps;

    /** Storage for T2 swap info.*/
    std::vector<T2SwapInfo<SPACE_DIM> > mT2Swaps;

    /** Storage for T3 swap info.*/
    std::vector<T3SwapInfo<SPACE_DIM> > mT3Swaps;

    /** Storage for cell division info.*/
    std::vector<CellDivisionInfo<SPACE_DIM> > mCellDivisions;

    /** Stores all mesh operations.*/
    std::vector<EdgeOperation> mEdgeOperations;

    /**
     * Pointer to edge handler.
     */
    EdgeHelper<SPACE_DIM>* mpEdgeHelper;

public:

    /**
     * Constructor.
     */
    VertexMeshOperationRecorder();

    /**
     * Destructor.
     */
    ~VertexMeshOperationRecorder();
    
    /**
     * Sets edge helper associated with the vertex mesh.
     * 
     * @param pEdgeHelper pointer to an edge helper
     */
    void SetEdgeHelper(EdgeHelper<SPACE_DIM>* pEdgeHelper);

    /**
     * Record T1 swap info.
     * 
     * @param rSwapInfo information about a T1 swap
     */
    void RecordT1Swap(T1SwapInfo<SPACE_DIM>& rSwapInfo);

    /**
     * @return Information about all T1 swaps
     */
    std::vector<T1SwapInfo<SPACE_DIM> > GetT1SwapsInfo() const;

    /**
     * Clear information about T1 swaps
     */
    void ClearT1SwapsInfo();

    /**
     * Record T2 swap info.
     * 
     * @param rSwapInfo information about a T2 swap
     */
    void RecordT2Swap(T2SwapInfo<SPACE_DIM>& rSwapInfo);

    /**
     * @return Information about T2 swaps
     */
    std::vector<T2SwapInfo<SPACE_DIM> > GetT2SwapsInfo() const;

    /**
     * Clear information about T1 swaps
     */
    void ClearT2SwapsInfo();

    /**
     * Record T3 swap info.
     * 
     * @param rSwapInfo information about a T3 swap
     */
    void RecordT3Swap(T3SwapInfo<SPACE_DIM>& rSwapInfo);

    /**
     * @return Information about T3 swaps
     */
    std::vector<T3SwapInfo<SPACE_DIM> > GetT3SwapsInfo() const;

    /**
     * Clear information about T3 swaps
     */
    void ClearT3SwapsInfo();

    /**
     * Record cell division event.
     * 
     * @param rDivisionInfo information about a cell division event.
     */
    void RecordCellDivisionInfo(CellDivisionInfo<SPACE_DIM>& rDivisionInfo);

    /**
     * @return Information about cell divisions
     */
    std::vector<CellDivisionInfo<SPACE_DIM> > GetCellDivisionInfo() const;

    /*
     * Clear information about T3 swaps.
     */
    void ClearCellDivisionInfo();

    /**
     * @return Edge operations since last clearing
     */
    const std::vector<EdgeOperation>& GetEdgeOperations();

    /**
     * Clear edge operations.
     */
    void ClearEdgeOperations();

    /**
     * Record node merging (or edge shrinkage) event.
     * 
     * @param oldIds Global Edge IDs prior to node merge
     * @param pElement Pointer to element associated with node merge
     * @param mergedNodesPair the index of the deleted node is stored in the 
     *                        second position
     * @param elementIndexIsRemapped indicates whether the operation has been 
     *                               recorded before element indices have been 
     *                               remapped, e.g. before elements are deleted 
     *                               (defaults to false)
     */
    void RecordNodeMergeOperation(const std::vector<unsigned> oldIds,
                                  VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                  const std::pair<unsigned, unsigned> mergedNodesPair,
                                  const bool elementIndexIsRemapped = false);

    /**
     * Record edge split operation in element pElement.
     * 
     * @param pElement Pointer to element associated with edge split
     * @param edgeIndex index of the edge being split
     * @param insertedNodeRelPosition position of the inserted node relative to 
     *                                the lower index node of the edge
     * @param elementIndexIsRemapped indicates whether the operation has been 
     *                               recorded before element indices have been 
     *                               remapped, e.g. before elements are deleted 
     *                               (defaults to false)
     */
    void RecordEdgeSplitOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                  const unsigned edgeIndex,
                                  const double insertedNodeRelPosition,
                                  const bool elementIndexIsRemapped = false);

    /**
     * Record cell divisions for VertexBasedPopulationSrn class to remap SRNs.
     * 
     * @param rOldIds Global Edge IDs of parent cell prior cell division
     * @param pElement1 Daughter element
     * @param pElement2 Daughter element
     */
    void RecordCellDivideOperation(const std::vector<unsigned>& rOldIds,
                                   VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement1,
                                   VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement2);

    /**
     * Records edge formation
     * @param pElement
     * @param edgeIndex - new edge index
     */
    void RecordNewEdgeOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                const unsigned edgeIndex);
    /**
     * Record merging of adjacent edges due to deletion of the shared node 
     * nodeIndex in element pElement
     * 
     * @param pElement Pointer to element associated with edge merge
     * @param nodeIndex Global index of node
     */
    void RecordEdgeMergeOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                  const unsigned nodeIndex);
};

#endif /* VERTEXMESHOPERATIONRECORDER_HPP_ */