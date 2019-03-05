//
// Created by twin on 08/02/19.
//

#include "CellEdgeDeltaNotchTrackingModifier.hpp"
#include "CellEdgeSrnModel.hpp"
#include "DeltaNotchSrnModel.hpp"

template<unsigned DIM>
CellEdgeDeltaNotchTrackingModifier<DIM>::CellEdgeDeltaNotchTrackingModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CellEdgeDeltaNotchTrackingModifier<DIM>::~CellEdgeDeltaNotchTrackingModifier()
{
}

template<unsigned DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
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

    // First recover each cell's Notch and Delta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_celledge_model = static_cast<CellEdgeSrnModel*>(cell_iter->GetSrnModel());

        std::vector<double> notches;
        std::vector<double> deltas;

        for(unsigned i = 0 ; i  < p_celledge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<DeltaNotchSrnModel> p_model = boost::static_pointer_cast<DeltaNotchSrnModel>(p_celledge_model->GetEdgeSrn(i));
            double this_delta = p_model->GetDelta();
            double this_notch = p_model->GetNotch();

            deltas.push_back(this_delta);
            notches.push_back(this_notch);

        }

        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("notch", notches);
        cell_iter->GetCellEdgeData()->SetItem("delta", deltas);

    }


    // Next iterate over the population to compute and store each cell's neighbouring Delta concentration in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

        // Compute this cell's average neighbouring Delta concentration and store in CellData
        if (!neighbour_indices.empty())
        {
            double mean_delta = 0.0;
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                double this_delta = p_cell->GetCellData()->GetItem("delta");
                mean_delta += this_delta/neighbour_indices.size();
            }
            cell_iter->GetCellData()->SetItem("mean delta", mean_delta);
        }
        else
        {
            // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
            cell_iter->GetCellData()->SetItem("mean delta", 0.0);
        }
    }
}

template<unsigned DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CellEdgeDeltaNotchTrackingModifier<1>;
template class CellEdgeDeltaNotchTrackingModifier<2>;
template class CellEdgeDeltaNotchTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellEdgeDeltaNotchTrackingModifier)
