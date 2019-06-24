#include "DeltaNotchEdgeSrnModel.hpp"
#include "CellEdgeDeltaNotchTrackingModifier.hpp"
#include <cell_based/src/population/VertexBasedCellPopulation.hpp>
#include "AbstractCellEdgeBasedSimulationModifier.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellEdges(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation)
{

    VertexBasedCellPopulation<SPACE_DIM>& vertexCellPopulation = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>&>(rCellPopulation);

    for(auto operation: vertexCellPopulation.GetCellEdgeChangeOperations())
    {
        switch(operation->GetOperation())
        {
            case EDGE_OPERATION_ADD:
            {
                //Gets the cell and the srn model
                auto cell = rCellPopulation.GetCellUsingLocationIndex(operation->GetElementIndex());
                auto srn_cell = static_cast<SrnCellModel*>(cell->GetSrnModel());
                //Adds a blank delta notch edge srn model
                auto new_srn_edge = boost::shared_ptr<AbstractSrnModel>(this->CreateEmptySrnEdgeModel());
                new_srn_edge->SetCell(srn_cell->GetCell());
                new_srn_edge->Initialise();
                srn_cell->InsertEdgeSrn(operation->GetLocalEdgeIndex(), new_srn_edge);
                //Calls the event handler
                this->EdgeAdded(rCellPopulation, operation->GetElementIndex(), operation->GetLocalEdgeIndex(), new_srn_edge);
            }
                break;

            case EDGE_OPERATION_DELETE:
            {
                //Gets the cell and the srn model
                auto cell = rCellPopulation.GetCellUsingLocationIndex(operation->GetElementIndex());
                auto srn_cell = static_cast<SrnCellModel*>(cell->GetSrnModel());
                //Removes the edge delta notch srn
                auto old_srn_edge = srn_cell->GetEdgeSrn(operation->GetLocalEdgeIndex());
                srn_cell->RemoveEdgeSrn(operation->GetLocalEdgeIndex());
                this->EdgeRemoved(rCellPopulation, operation->GetElementIndex(), operation->GetLocalEdgeIndex(), old_srn_edge);
            }
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

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::CellDivisionEdgeUpdate(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation, unsigned locationIndex,
        EdgeRemapInfo *edgeChange, unsigned locationIndex2, EdgeRemapInfo *edgeChange2) {

    //Gets the cell
    auto cell1 = rCellPopulation.GetCellUsingLocationIndex(locationIndex);
    auto cell2 = rCellPopulation.GetCellUsingLocationIndex(locationIndex2);

    //Gets the edge srn model and a copy of the srn models in the edges
    auto old_model = static_cast<SrnCellModel*>(cell1->GetSrnModel());
    std::vector<AbstractSrnModelPtr> old_srn_edges = old_model->GetEdges();

    auto srn_model1 = new SrnCellModel();
    auto srn_model2 = new SrnCellModel();

    cell1->SetSrnModel(srn_model1);
    cell2->SetSrnModel(srn_model2);

    PerformEdgeRemap(rCellPopulation, locationIndex, old_srn_edges, srn_model1, edgeChange);
    PerformEdgeRemap(rCellPopulation, locationIndex2, old_srn_edges, srn_model2, edgeChange2);
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::PerformEdgeRemap(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation,
        unsigned locationIndex,
        std::vector<AbstractSrnModelPtr>& old_edges, SrnCellModel *newModel,
        EdgeRemapInfo *edgeChange)
{
    //Goes through the SRN model
    for(unsigned i = 0 ; i < edgeChange->GetEdgesMapping().size(); i++)
    {
        //The remapIndex, if +ve refers to the srn index of the oldModel, if -ve then it's a new edge
        auto remapIndex = edgeChange->GetEdgesMapping()[i];

        //remapStatus can be the following:
        //0 - Direct remapping, the edge srn of the oldModel can be transferred directly to the new model
        //1 - The edge is a split point between the diving cells, in this example we divide all concentration in half
        //2 - This is a new edge i.e. the dividing line in the middle of the old and new cells
        auto remapStatus = edgeChange->GetEdgesStatus()[i];

        if((remapStatus == 0 || remapStatus == 1) && remapIndex < 0)
            EXCEPTION("Remap index cannot be negative when it's a direct remap or an edge split");

        switch(remapStatus)
        {
            /* Direct remap - move the edge to the new srn cell model */
            case 0:
            {
                auto current_edge_srn = old_edges[remapIndex];
                current_edge_srn->SetCell(newModel->GetCell()); //Move the edge to the current cell
                newModel->AddEdgeSrn(current_edge_srn);
            }
                break;

            /* Split or New edge, either way a new srn edge is created */
            case 1:
            case 2:
            {
                auto new_edge_srn = boost::shared_ptr<AbstractSrnModel>(this->CreateEmptySrnEdgeModel());
                new_edge_srn->SetCell(newModel->GetCell()); //New model does not have an existing cell so we must set this
                newModel->AddEdgeSrn(new_edge_srn);
                new_edge_srn->Initialise();
            }
                break;

        }
    }

    //Goes through the SRN model
    for(unsigned i = 0 ; i < edgeChange->GetEdgesMapping().size(); i++)
    {
        //The remapIndex, if +ve refers to the srn index of the oldModel, if -ve then it's a new edge
        auto remapIndex = edgeChange->GetEdgesMapping()[i];

        //remapStatus can be the following:
        //0 - Direct remapping, the edge srn of the oldModel can be transferred directly to the new model
        //1 - The edge is a split point between the diving cells, in this example we divide all concentration in half
        //2 - This is a new edge i.e. the dividing line in the middle of the old and new cells
        auto remapStatus = edgeChange->GetEdgesStatus()[i];

        if((remapStatus == 0 || remapStatus == 1) && remapIndex < 0)
            EXCEPTION("Remap index cannot be negative when it's a direct remap or an edge split");

        switch(remapStatus)
        {
            /* Direct remap - move the edge to the new srn cell model */
            case 0:
            {
                /* Don't do anything */
            }
                break;

                /* Split - Divide the concentration in half */
            case 1:
            {
                this->EdgeDivide(old_edges[remapIndex], newModel->GetEdgeSrn(i));
            }
                break;

                /* New edge - Initialise new edge model */
            case 2:
            {
                this->EdgeAdded(rCellPopulation, locationIndex, i, newModel->GetEdgeSrn(i));
            }
                break;

        }
    }

}



// Explicit instantiation
template class AbstractCellEdgeBasedSimulationModifier<1,1>;
template class AbstractCellEdgeBasedSimulationModifier<2,2>;
template class AbstractCellEdgeBasedSimulationModifier<3,3>;