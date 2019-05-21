#include <cell_based/src/population/VertexBasedCellPopulation.hpp>
#include "AbstractCellEdgeBasedSimulationModifier.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellEdges(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation) {

    VertexBasedCellPopulation<SPACE_DIM>& vertexCellPopulation = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>&>(rCellPopulation);

    for(auto operation: vertexCellPopulation.GetCellEdgeChangeOperations())
    {
        switch(operation->GetOperation())
        {
            case EDGE_OPERATION_ADD:
                this->EdgeAdded(rCellPopulation, operation->GetElementIndex(), operation->GetLocalEdgeIndex());
                break;

            case EDGE_OPERATION_DELETE:
                this->EdgeRemoved(rCellPopulation, operation->GetElementIndex(), operation->GetLocalEdgeIndex());
                break;

            case EDGE_OPERATION_DIVIDE:
                this->CellDivisionEdgeUpdate(rCellPopulation, operation->GetElementIndex(), operation->GetNewEdges(), operation->GetElementIndex2(),  operation->GetNewEdges2());
                break;

            default:
                EXCEPTION("Invalid edge operation");
        }
    }

    vertexCellPopulation.ClearCellEdgeOperations();

}

// Explicit instantiation
template class AbstractCellEdgeBasedSimulationModifier<1,1>;
template class AbstractCellEdgeBasedSimulationModifier<2,2>;
template class AbstractCellEdgeBasedSimulationModifier<3,3>;