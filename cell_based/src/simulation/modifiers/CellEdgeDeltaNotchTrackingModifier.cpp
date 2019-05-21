//
// Created by twin on 08/02/19.
//

#include "CellEdgeDeltaNotchTrackingModifier.hpp"
#include "CellEdgeSrnModel.hpp"
#include "DeltaNotchEdgeSrnModel.hpp"

template<unsigned DIM>
CellEdgeDeltaNotchTrackingModifier<DIM>::CellEdgeDeltaNotchTrackingModifier()
        : AbstractCellEdgeBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CellEdgeDeltaNotchTrackingModifier<DIM>::~CellEdgeDeltaNotchTrackingModifier()
{
}

template<unsigned DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    AbstractCellEdgeBasedSimulationModifier<DIM,DIM>::UpdateCellEdges(rCellPopulation);
//    this->UpdateCellEdges(rCellPopulation);

    //Updates the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Handles all edge changes
    this->UpdateCellEdges(rCellPopulation);

    // First recover each cell's Notch and Delta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<CellEdgeSrnModel*>(cell_iter->GetSrnModel());

        std::vector<double> notch_vec;
        std::vector<double> delta_vec;

        for(unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<DeltaNotchEdgeSrnModel> p_model = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(p_cell_edge_model->GetEdgeSrn(i));
            double this_delta = p_model->GetDelta();
            double this_notch = p_model->GetNotch();

            delta_vec.push_back(this_delta);
            notch_vec.push_back(this_notch);

        }

        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("notch", notch_vec);
        cell_iter->GetCellEdgeData()->SetItem("delta", delta_vec);

    }


    // Next iterate over the population to compute and store each cell's neighbouring Delta concentration in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<CellEdgeSrnModel*>(cell_iter->GetSrnModel());
        std::vector<double> delta_vec = cell_iter->GetCellEdgeData()->GetItem("delta");
        std::vector<double> mean_delta_vec(p_cell_edge_model->GetNumEdgeSrn());

        // Iterate through each Edge SRN
        for(unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            double mean_delta = 0.0;
            double delta_counter = 0;

            // Edge neighbour outside of cell
            auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(*cell_iter, i);
            for(auto neighbourIndex: elemNeighbours)
            {
                auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);

                std::vector<double> neighbour_delta_vec = neighbourCell->GetCellEdgeData()->GetItem("delta");

                mean_delta += neighbour_delta_vec[neighbourIndex.second];
                delta_counter ++;

            }

            // Edge neighbour inside of cell
            // In this case we're only looking at immediate neighbour to the current edge
            auto prevNeighbour = ((int)i -1) % p_cell_edge_model->GetNumEdgeSrn();
            auto nextNeighbour = (i +1) % p_cell_edge_model->GetNumEdgeSrn();
            delta_counter += 2;
            mean_delta += delta_vec[prevNeighbour];
            mean_delta += delta_vec[nextNeighbour];

            mean_delta = mean_delta/delta_counter;

            // Stores the delta
            mean_delta_vec[i] = mean_delta;

        }

        cell_iter->GetCellEdgeData()->SetItem("mean delta", mean_delta_vec);
    }
}

template<unsigned DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned int DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::EdgeAdded(AbstractCellPopulation<DIM,DIM> &rCellPopulation,
                                                        unsigned locationIndex, unsigned edgeLocalIndex) {

    //Gets the cell
    auto cell = rCellPopulation.GetCellUsingLocationIndex(locationIndex);

    //Gets the edge srn model
    auto p_cell_edge_model = static_cast<CellEdgeSrnModel*>(cell->GetSrnModel());

    //Adds a blank delta notch edge srn model
    boost::shared_ptr<AbstractOdeSrnModel> delta_notch_srn(new DeltaNotchEdgeSrnModel());
    p_cell_edge_model->InsertEdgeSrn(edgeLocalIndex, delta_notch_srn);




}

template<unsigned int DIM>
void
CellEdgeDeltaNotchTrackingModifier<DIM>::EdgeRemoved(AbstractCellPopulation<DIM,DIM> &rCellPopulation,
                                                     unsigned locationIndex, unsigned edgeLocalIndex) {

    //Gets the cell
    auto cell = rCellPopulation.GetCellUsingLocationIndex(locationIndex);

    //Gets the edge srn model
    auto p_cell_edge_model = static_cast<CellEdgeSrnModel*>(cell->GetSrnModel());

    //Removes the edge delta notch srn
    p_cell_edge_model->RemoveEdgeSrn(edgeLocalIndex);

}

template<unsigned int DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::CellDivisionEdgeUpdate(
        AbstractCellPopulation<DIM,DIM>& rCellPopulation,
        unsigned locationIndex, EdgeRemapInfo* edgeChange,
        unsigned locationIndex2, EdgeRemapInfo* edgeChange2) {

    //Gets the cell
    auto cell = rCellPopulation.GetCellUsingLocationIndex(locationIndex);
    auto cell2 = rCellPopulation.GetCellUsingLocationIndex(locationIndex2);



    //Gets the edge srn model
    auto old_model = static_cast<CellEdgeSrnModel*>(cell->GetSrnModel());
    auto model1 = new CellEdgeSrnModel();
    auto model2 = new CellEdgeSrnModel();

    PerformEdgeRemap(old_model, model1, edgeChange);
    PerformEdgeRemap(old_model, model2, edgeChange2);

    cell->SetSrnModel(model1);
    cell->SetSrnModel(model2);


}

template<unsigned int DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::PerformEdgeRemap(CellEdgeSrnModel *oldModel, CellEdgeSrnModel *newModel,
                                                               EdgeRemapInfo *edgeChange) {

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
            //Direct remap
            case 0:
            {
                auto current_edge_srn = oldModel->GetEdgeSrn(remapIndex);
                newModel->AddEdgeSrn(current_edge_srn);
            }

                break;

            //Split - Divide the concentration in half
            case 1:
            {
                auto current_edge_srn = boost::dynamic_pointer_cast<DeltaNotchEdgeSrnModel>(oldModel->GetEdgeSrn(remapIndex));
                boost::shared_ptr<DeltaNotchEdgeSrnModel> p_srn_model(new DeltaNotchEdgeSrnModel());

                p_srn_model->SetDelta(current_edge_srn->GetDelta()/2.0);
                p_srn_model->SetNotch(current_edge_srn->GetNotch()/2.0);

                newModel->AddEdgeSrn(p_srn_model);
            }
                break;

            //New edge - Initialise new model
            case 2:
                boost::shared_ptr<DeltaNotchEdgeSrnModel> p_srn_model(new DeltaNotchEdgeSrnModel());
                newModel->AddEdgeSrn(p_srn_model);
                break;

        }
    }

}


// Explicit instantiation
template class CellEdgeDeltaNotchTrackingModifier<1>;
template class CellEdgeDeltaNotchTrackingModifier<2>;
template class CellEdgeDeltaNotchTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellEdgeDeltaNotchTrackingModifier)
