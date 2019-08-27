#include <cell_based/src/population/VertexBasedCellPopulation.hpp>
#include "AbstractCellEdgeBasedSimulationModifier.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellSrnLayout(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation)
{

    VertexBasedCellPopulation<SPACE_DIM>& vertexCellPopulation
    = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>&>(rCellPopulation);

    for (auto operation : vertexCellPopulation.GetCellEdgeChangeOperations())
    {
        switch (operation->GetOperation())
        {
            case EDGE_OPERATION_ADD:
            {
                // Get the cell and the SRN model
                auto cell = rCellPopulation.GetCellUsingLocationIndex(operation->GetElementIndex());
                auto srn_cell = static_cast<SrnCellModel*>(cell->GetSrnModel());

                // Add a blank delta notch edge SRN model
                auto new_srn_edge = boost::shared_ptr<AbstractSrnModel>(this->CreateEmptySrnEdgeModel());
                new_srn_edge->SetCell(srn_cell->GetCell());
                new_srn_edge->Initialise();
                srn_cell->InsertEdgeSrn(operation->GetLocalEdgeIndex(), new_srn_edge);

                // Call the event handler
                this->EdgeAdded(rCellPopulation, operation->GetElementIndex(), operation->GetLocalEdgeIndex(), new_srn_edge);
                break;
            }
            case EDGE_OPERATION_DELETE:
            {
                // Get the cell and the SRN model
                auto cell = rCellPopulation.GetCellUsingLocationIndex(operation->GetElementIndex());
                auto srn_cell = static_cast<SrnCellModel*>(cell->GetSrnModel());

                // Remove the edge delta notch srn
                auto old_srn_edge = srn_cell->GetEdgeSrn(operation->GetLocalEdgeIndex());
                srn_cell->RemoveEdgeSrn(operation->GetLocalEdgeIndex());
                this->EdgeRemoved(rCellPopulation, operation->GetElementIndex(), operation->GetLocalEdgeIndex(), old_srn_edge);
                break;
            }
            case EDGE_OPERATION_DIVIDE:
                break;
            default:
            {
                EXCEPTION("Invalid edge operation");
            }
        }
    }

    vertexCellPopulation.ClearCellEdgeOperations();
}


// Explicit instantiation
template class AbstractCellEdgeBasedSimulationModifier<1,1>;
template class AbstractCellEdgeBasedSimulationModifier<2,2>;
template class AbstractCellEdgeBasedSimulationModifier<3,3>;
