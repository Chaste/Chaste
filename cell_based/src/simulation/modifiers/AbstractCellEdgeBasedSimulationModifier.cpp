#include "DeltaNotchSrnEdgeModel.hpp"
#include "DeltaNotchCellEdgeTrackingModifier.hpp"
#include <cell_based/src/population/VertexBasedCellPopulation.hpp>
#include "AbstractCellEdgeBasedSimulationModifier.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellSrnLayout(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation)
{

    VertexBasedCellPopulation<SPACE_DIM>& vertexCellPopulation = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>&>(rCellPopulation);

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
            {
                this->CellDivisionEdgeUpdate(rCellPopulation, operation->GetElementIndex(), operation->GetNewEdges(), operation->GetElementIndex2(),  operation->GetNewEdges2());
                break;
            }
            default:
            {
                EXCEPTION("Invalid edge operation");
            }
        }
    }

    vertexCellPopulation.ClearCellEdgeOperations();
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::CellDivisionEdgeUpdate(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation, unsigned locationIndex,
        EdgeRemapInfo *pEdgeChange, unsigned locationIndex2, EdgeRemapInfo *pEdgeChange2)
{
    // Get the cell
    auto cell1 = rCellPopulation.GetCellUsingLocationIndex(locationIndex);
    auto cell2 = rCellPopulation.GetCellUsingLocationIndex(locationIndex2);

    // Ges the edge SRN model and a copy of the SRN models in the edges and interior
    auto old_model = static_cast<SrnCellModel*>(cell1->GetSrnModel());
    auto old_srn_interior = old_model->GetInteriorSrn();
    std::vector<AbstractSrnModelPtr> old_srn_edges = old_model->GetEdges();

    auto srn_model1 = new SrnCellModel();
    auto srn_model2 = new SrnCellModel();

    cell1->SetSrnModel(srn_model1);
    cell2->SetSrnModel(srn_model2);

    PerformEdgeRemap(rCellPopulation, locationIndex, old_srn_edges, srn_model1, pEdgeChange);
    PerformEdgeRemap(rCellPopulation, locationIndex2, old_srn_edges, srn_model2, pEdgeChange2);

    //If there's interior srn then also allow custom srn split
    if (old_srn_interior != nullptr)
    {
        srn_model1->SetInteriorSrnModel(boost::shared_ptr<AbstractSrnModel>(this->CreateEmptySrnInteriorModel()));
        srn_model1->GetInteriorSrn()->SetCell(cell1);
        srn_model1->GetInteriorSrn()->Initialise();
        srn_model2->SetInteriorSrnModel(boost::shared_ptr<AbstractSrnModel>(this->CreateEmptySrnInteriorModel()));
        srn_model2->GetInteriorSrn()->SetCell(cell2);
        srn_model2->GetInteriorSrn()->Initialise();
        this->InteriorDivide(old_srn_interior, srn_model1->GetInteriorSrn());
        this->InteriorDivide(old_srn_interior, srn_model2->GetInteriorSrn());
    }



}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void AbstractCellEdgeBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::PerformEdgeRemap(
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation,
        unsigned locationIndex,
        std::vector<AbstractSrnModelPtr>& old_edges,
        SrnCellModel *pNewModel,
        EdgeRemapInfo *pEdgeChange)
{
    // Goes through the SRN model
    for (unsigned i = 0; i < pEdgeChange->GetEdgesMapping().size(); i++)
    {
        //The remapIndex, if +ve refers to the SRN index of the oldModel, if -ve then it's a new edge
        auto remapIndex = pEdgeChange->GetEdgesMapping()[i];

        //remapStatus can be the following:
        //0 - Direct remapping, the edge SRN of the oldModel can be transferred directly to the new model
        //1 - The edge is a split point between the diving cells, in this example we divide all concentration in half
        //2 - This is a new edge i.e. the dividing line in the middle of the old and new cells
        auto remapStatus = pEdgeChange->GetEdgesStatus()[i];

        if ((remapStatus == 0 || remapStatus == 1) && remapIndex < 0)
        {
            EXCEPTION("Remap index cannot be negative when it's a direct remap or an edge split");
        }

        if (remapStatus == 0)
        {
            /* Direct remap - move the edge to the new SRN cell model */
            auto current_edge_srn = old_edges[remapIndex];
            current_edge_srn->SetCell(pNewModel->GetCell()); //Move the edge to the current cell
            pNewModel->AddEdgeSrnModel(current_edge_srn);
        }
        else
        {
            /* Split or new edge, either way a new SRN edge is created */
            auto new_edge_srn = boost::shared_ptr<AbstractSrnModel>(this->CreateEmptySrnEdgeModel());
            new_edge_srn->SetCell(pNewModel->GetCell()); //New model does not have an existing cell so we must set this
            pNewModel->AddEdgeSrnModel(new_edge_srn);
            new_edge_srn->Initialise();
        }
    }

    // Goes through the SRN model
    for (unsigned i = 0; i < pEdgeChange->GetEdgesMapping().size(); i++)
    {
        //The remapIndex, if +ve refers to the SRN index of the oldModel, if -ve then it's a new edge
        auto remapIndex = pEdgeChange->GetEdgesMapping()[i];

        //remapStatus can be the following:
        //0 - Direct remapping, the edge SRN of the oldModel can be transferred directly to the new model
        //1 - The edge is a split point between the diving cells, in this example we divide all concentration in half
        //2 - This is a new edge i.e. the dividing line in the middle of the old and new cells
        auto remapStatus = pEdgeChange->GetEdgesStatus()[i];

        if ((remapStatus == 0 || remapStatus == 1) && remapIndex < 0)
        {
            EXCEPTION("Remap index cannot be negative when it's a direct remap or an edge split");
        }

        switch (remapStatus)
        {
            case 0:
            {
                /* Direct remap - move the edge to the new SRN cell model */
                break;
            }
            case 1:
            {
                /* Split - e.g. Divide the concentration in half */
                this->EdgeDivide(old_edges[remapIndex], pNewModel->GetEdgeSrn(i));
                break;
            }
            case 2:
            {
                /* New edge - Initialise new edge model */
                this->EdgeAdded(rCellPopulation, locationIndex, i, pNewModel->GetEdgeSrn(i));
                break;
            }
        }
    }
}

// Explicit instantiation
template class AbstractCellEdgeBasedSimulationModifier<1,1>;
template class AbstractCellEdgeBasedSimulationModifier<2,2>;
template class AbstractCellEdgeBasedSimulationModifier<3,3>;
