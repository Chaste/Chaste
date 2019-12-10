
#include "VertexMeshOperationRecorder.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::VertexMeshOperationRecorder()
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::~VertexMeshOperationRecorder()
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::SetEdgeHelper(EdgeHelper<SPACE_DIM> *pEdgeHelper)
{
    mpEdgeHelper = pEdgeHelper;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetLocationsOfT1Swaps() const
{
    return mLocationsOfT1Swaps;
}

template<unsigned ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearLocationsOfT1Swaps()
{
    mLocationsOfT1Swaps.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetLocationsOfT2Swaps() const
{
    return mLocationsOfT2Swaps;
}

template<unsigned ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearLocationsOfT2Swaps()
{
    mLocationsOfT2Swaps.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::GetLocationsOfT3Swaps() const
{
    return mLocationsOfT3Swaps;
}

template<unsigned ELEMENT_DIM, unsigned int SPACE_DIM>
void VertexMeshOperationRecorder<ELEMENT_DIM, SPACE_DIM>::ClearLocationsOfT3Swaps()
{
    mLocationsOfT3Swaps.clear();
}

