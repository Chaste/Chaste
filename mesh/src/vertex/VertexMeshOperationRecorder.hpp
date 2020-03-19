/*
 * VertexMeshOperationRecorder.hpp
 *
 *  Created on: 10 Dec 2019
 *      Author: aydar
 */

#ifndef VERTEXMESHOPERATIONRECORDER_HPP_
#define VERTEXMESHOPERATIONRECORDER_HPP_

#include "VertexElement.hpp"
#include "EdgeHelper.hpp"
#include "EdgeOperation.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * Contains information about a T1 swap.
 * mLocation - central point of the shrinking edge
 * mPreSwapEdge - orientation of the shrinking edge
 * mPostSwapEdge - orientation of the post swap edge
 */
template<unsigned int SPACE_DIM>
struct T1SwapInfo
{
    c_vector<double, SPACE_DIM> mLocation;
    c_vector<double, SPACE_DIM> mPreSwapEdge;
    c_vector<double, SPACE_DIM> mPostSwapEdge;
};
/**
 * Contains information about a T2 swap.
 * mLocation - centroid of the element undergoing T2 swap
 */
template<unsigned int SPACE_DIM>
struct T2SwapInfo
{
    unsigned int mCellId;
    c_vector<double, SPACE_DIM> mLocation;
};

template<unsigned int SPACE_DIM>
struct T3SwapInfo
{
    c_vector<double, SPACE_DIM> mLocation;
};

/**
 * Contains information about cell division event
 * mLocation - centroid of the mother cell
 * mDaughtaurLocation1, mDaughtaurLocation1 - centroids of the two daughter cells
 * mDaughterOrientation1, mDaughterOrientation2 - orientation of the two daughter cells
 * mDivisionAxis - orientation of the division axis
 */
template<unsigned int SPACE_DIM>
struct CellDivisionInfo
{
    c_vector<double, SPACE_DIM> mLocation;
    c_vector<double, SPACE_DIM> mDaughterLocation1;
    c_vector<double, SPACE_DIM> mDaughterLongAxis1;
    c_vector<double, SPACE_DIM> mDaughterLocation2;
    c_vector<double, SPACE_DIM> mDaughterLongAxis2;
    c_vector<double, SPACE_DIM> mDivisionAxis;
};

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
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
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        /*archive & mT1Swaps;
        archive & mT2Swaps;
        archive & mT3Swaps;
        archive & mCellDivisions;*/
    }
    std::vector<T1SwapInfo<SPACE_DIM> > mT1Swaps;
    std::vector<T2SwapInfo<SPACE_DIM> > mT2Swaps;
    std::vector<T3SwapInfo<SPACE_DIM> > mT3Swaps;
    std::vector<CellDivisionInfo<SPACE_DIM> > mCellDivisions;

    std::vector<EdgeOperation*> mEdgeOperations;

    EdgeHelper<SPACE_DIM> *mpEdgeHelper;

public:
    VertexMeshOperationRecorder();
    ~VertexMeshOperationRecorder();

    /**
     * Sets edge helper associated with the vertex mesh
     * @param pEdgeHelper
     */
    void SetEdgeHelper(EdgeHelper<SPACE_DIM> *pEdgeHelper);

    /**
     * Record T1 swap info
     * @param rSwap_info - information about T1 swap
     */
    void RecordT1Swap(T1SwapInfo<SPACE_DIM>& rSwap_info);

    /**
     * @return Information about all T1 swaps
     */
    std::vector<T1SwapInfo<SPACE_DIM> > GetT1SwapsInfo() const;

    /**
     * Clear information about T1 swaps
     */
    void ClearT1SwapsInfo();

    /**
     * Record T2 swap info
     * @param rSwap_info
     */
    void RecordT2Swap(T2SwapInfo<SPACE_DIM>& rSwap_info);

    /**
     * @return Information about T2 swaps
     */
    std::vector<T2SwapInfo<SPACE_DIM> > GetT2SwapsInfo() const;

    /**
     * Clear information about T1 swaps
     */
    void ClearT2SwapsInfo();

    /**
     * Record T3 swap info
     * @param rSwap_info
     */
    void RecordT3Swap(T3SwapInfo<SPACE_DIM>& rSwap_info);

    /**
     * @return Information about T3 swaps
     */
    std::vector<T3SwapInfo<SPACE_DIM> > GetT3SwapsInfo() const;

    /**
     * Clear information about T3 swaps
     */
    void ClearT3SwapsInfo();

    /**
     * Record cell division event
     * @param rSwap_info
     */
    void RecordCellDivisionInfo(CellDivisionInfo<SPACE_DIM>& rDivision_info);

    /**
     * @return Information about cell divisions
     */
    std::vector<CellDivisionInfo<SPACE_DIM> > GetCellDivisionInfo() const;

    /**
     * Clear information about T3 swaps
     */
    void ClearCellDivisionInfo();

    /**
     * @return Edge operations since last clearing
     */
    const std::vector<EdgeOperation*> & GetEdgeOperations();

    /**
     * Clears edge operations
     */
    void ClearEdgeOperations();

    /**
     * Records node merging (or edge shrinkage) event
     * @param oldIds
     * @param pElement
     * @param merged_nodes_pair - the index of the deleted node is stored in the second position
     * @param elementIndexIsRemapped - indicates whether the operation has been recorded before
     *  element indices have been remapped (e.g. before elements are deleted)
     */
    void RecordNodeMergeOperation(const std::vector<unsigned int> oldIds,
                                  VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                  const std::pair<unsigned int, unsigned int> merged_nodes_pair,
                                  const bool elementIndexIsRemapped = false);

    /**
     * Records edge split operation in element pElement.
     * @param pElement
     * @param edge_index - index of the edge being split
     * @param inserted_node_rel_position - position of the inserted node relative to the lower index node of the edge
     * @param elementIndexIsRemapped - indicates whether the operation has been recorded before
     *  element indices have been remapped (e.g. before elements are deleted)
     */
    void RecordEdgeSplitOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                  const unsigned int edge_index,
                                  const double inserted_node_rel_position,
                                  const bool elementIndexIsRemapped = false);

    /**
     * Records cell divisions for VertexBasedPopulationSrn class to remap SRNs
     * @param oldIds Global Edge IDs of parent cell prior cell division
     * @param pElement1 Daughter element
     * @param pElement2 Daughter element
     */
    void RecordCellDivideOperation(const std::vector<unsigned int>& oldIds,
                                   VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement1,
                                   VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement2);

    /**
     * Records edge formation
     * @param pElement
     * @param edge_index - new edge index
     */
    void RecordNewEdgeOperation(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement,
                                const unsigned int edge_index);
    /**
     * Records merging of adjacent edges due to deletion of the shared node node_index in element pElement
     * @param pElement
     * @param edge
     */
    void RecordEdgeMergeOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                  const unsigned int node_index);
};

#endif /* VERTEXMESHOPERATIONRECORDER_HPP_ */
