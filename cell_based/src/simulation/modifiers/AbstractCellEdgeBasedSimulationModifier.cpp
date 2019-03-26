#include <cell_based/src/population/VertexBasedCellPopulation.hpp>
#include "AbstractCellEdgeBasedSimulationModifier.hpp"

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation) {

    auto vertexBasedCellPopulation = static_cast<VertexBasedCellPopulation<SPACE_DIM>&>(rCellPopulation);

    auto mesh = vertexBasedCellPopulation.rGetMesh();

    //TODO: Goes through the edge change list and broadcast to cells


}

