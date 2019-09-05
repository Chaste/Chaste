#include <cell_based/src/population/VertexBasedCellPopulation.hpp>
#include "AbstractCellEdgeBasedSimulationModifier.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellSrnLayout(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation)
{

}


// Explicit instantiation
template class AbstractCellEdgeBasedSimulationModifier<1,1>;
template class AbstractCellEdgeBasedSimulationModifier<2,2>;
template class AbstractCellEdgeBasedSimulationModifier<3,3>;
