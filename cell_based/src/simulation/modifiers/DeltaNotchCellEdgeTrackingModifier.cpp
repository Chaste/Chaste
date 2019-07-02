//
// Created by twin on 08/02/19.
//

#include "DeltaNotchCellEdgeTrackingModifier.hpp"
#include "SrnCellModel.hpp"
#include "DeltaNotchSrnEdgeModel.hpp"

template<unsigned DIM>
DeltaNotchCellEdgeTrackingModifier<DIM>::DeltaNotchCellEdgeTrackingModifier()
        : AbstractCellEdgeBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
DeltaNotchCellEdgeTrackingModifier<DIM>::~DeltaNotchCellEdgeTrackingModifier()
{
}

template<unsigned DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    AbstractCellEdgeBasedSimulationModifier<DIM,DIM>::UpdateCellEdges(rCellPopulation);

    //Updates the cell
    this->UpdateCellData(rCellPopulation);

    //Debug cell concentration
//    printf("Time: %f ", SimulationTime::Instance()->GetTime());
//    double totalConcentration = 0;
//    for(unsigned i = 0 ; i < rCellPopulation.GetNumAllCells(); i ++)
//    {
//        auto cell = rCellPopulation.GetCellUsingLocationIndex(i);
//        //auto edgesSrn = static_cast<SrnCellModel*>(cell->GetSrnModel());
//        auto notchConcentration = cell->GetCellEdgeData()->GetItem("notch");
//
//
//
//        printf("|| NC %i:", i);
//        for( auto notch: notchConcentration)
//        {
//            totalConcentration += notch;
//            printf(" %f ", notch);
//        }
//    }
//    printf(" || TC %f", totalConcentration);
//    printf("\n");

}

template<unsigned DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        auto p_cell_edge_model = static_cast<SrnCellModel*>(cell_iter->GetSrnModel());

        std::vector<double> notch_vec;
        std::vector<double> delta_vec;

        for(unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<DeltaNotchSrnEdgeModel> p_model = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_edge_model->GetEdgeSrn(i));
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
        auto p_cell_edge_model = static_cast<SrnCellModel*>(cell_iter->GetSrnModel());
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
            auto prevNeighbour = ((int)i -1);
            if(prevNeighbour < 0)
                prevNeighbour = p_cell_edge_model->GetNumEdgeSrn() -1;
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
void DeltaNotchCellEdgeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned int DIM>
AbstractSrnModel *DeltaNotchCellEdgeTrackingModifier<DIM>::CreateEmptySrnEdgeModel() {
    return new DeltaNotchSrnEdgeModel();
}

template<unsigned int DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::EdgeAdded(AbstractCellPopulation<DIM,DIM> &rCellPopulation,
                                                        unsigned locationIndex, unsigned edgeLocalIndex, AbstractSrnModelPtr addedEdge) {

    auto deltaNotchNewSrnEdge = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(addedEdge);
    //New edges have a concentration of 0
    deltaNotchNewSrnEdge->SetDelta(0);
    deltaNotchNewSrnEdge->SetNotch(0);

}

template<unsigned int DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::EdgeRemoved(AbstractCellPopulation<DIM,DIM> &rCellPopulation,
                                                     unsigned locationIndex, unsigned edgeLocalIndex, AbstractSrnModelPtr oldSrnEdge) {



}

template<unsigned int DIM>
void DeltaNotchCellEdgeTrackingModifier<DIM>::EdgeDivide(AbstractSrnModelPtr oldSrnEdge, AbstractSrnModelPtr newSrnEdge) {

    //Convert to delta notch srn type
    auto deltaNotchOldSrnEdge = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(oldSrnEdge);
    auto deltaNotchNewSrnEdge = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(newSrnEdge);

    //In this example we're just halving the concentrations
    deltaNotchNewSrnEdge->SetDelta(deltaNotchOldSrnEdge->GetDelta()/2.0);
    deltaNotchNewSrnEdge->SetNotch(deltaNotchOldSrnEdge->GetNotch()/2.0);
}




// Explicit instantiation
template class DeltaNotchCellEdgeTrackingModifier<1>;
template class DeltaNotchCellEdgeTrackingModifier<2>;
template class DeltaNotchCellEdgeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
#include "AbstractCellEdgeBasedSimulationModifier.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchCellEdgeTrackingModifier)
