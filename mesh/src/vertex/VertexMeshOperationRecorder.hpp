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



template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
class VertexMeshOperationRecorder
{
private:
    std::vector<c_vector<double, SPACE_DIM > > mLocationsOfT1Swaps;
    std::vector<c_vector<double, SPACE_DIM > > mLocationsOfT2Swaps;
    std::vector<c_vector<double, SPACE_DIM > > mLocationsOfT3Swaps;
    std::vector<EdgeOperation*> mEdgeOperations;

    EdgeHelper<SPACE_DIM> *mpEdgeHelper;
public:
    VertexMeshOperationRecorder();
    ~VertexMeshOperationRecorder();

    void SetEdgeHelper(EdgeHelper<SPACE_DIM> *pEdgeHelper);

    void InsertT1SwapLocation(const c_vector<double, SPACE_DIM> location);
    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT1Swaps() const;
    void ClearLocationsOfT1Swaps();

    void InsertT2SwapLocation(const c_vector<double, SPACE_DIM> location);
    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT2Swaps() const;
    void ClearLocationsOfT2Swaps();

    void InsertT3SwapLocation(const c_vector<double, SPACE_DIM> location);
    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT3Swaps() const;
    void ClearLocationsOfT3Swaps();

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
     */
    void RecordEdgeSplitOperation(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                  const unsigned int edge_index,
                                  const double inserted_node_rel_position);

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
