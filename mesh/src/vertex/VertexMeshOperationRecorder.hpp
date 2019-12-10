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

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
class VertexMeshOperationRecorder
{
private:
    std::vector<c_vector<double, SPACE_DIM > > mLocationsOfT1Swaps;
    std::vector<c_vector<double, SPACE_DIM > > mLocationsOfT2Swaps;
    std::vector<c_vector<double, SPACE_DIM > > mLocationsOfT3Swaps;

    EdgeHelper<SPACE_DIM> *mpEdgeHelper;
public:
    VertexMeshOperationRecorder();
    ~VertexMeshOperationRecorder();

    void SetEdgeHelper(EdgeHelper<SPACE_DIM> *pEdgeHelper);

    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT1Swaps() const;
    void ClearLocationsOfT1Swaps();
    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT2Swaps() const;
    void ClearLocationsOfT2Swaps();
    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT3Swaps() const;
    void ClearLocationsOfT3Swaps();

};

#endif /* VERTEXMESHOPERATIONRECORDER_HPP_ */
