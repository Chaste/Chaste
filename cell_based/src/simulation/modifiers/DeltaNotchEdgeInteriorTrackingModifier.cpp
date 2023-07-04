/*

Copyright (c) 2005-2023, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "DeltaNotchEdgeInteriorTrackingModifier.hpp"
#include "CellSrnModel.hpp"
#include "DeltaNotchInteriorSrnModel.hpp"
#include "DeltaNotchEdgeSrnModel.hpp"

template<unsigned DIM>
DeltaNotchEdgeInteriorTrackingModifier<DIM>::DeltaNotchEdgeInteriorTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
DeltaNotchEdgeInteriorTrackingModifier<DIM>::~DeltaNotchEdgeInteriorTrackingModifier()
{
}

template<unsigned DIM>
void DeltaNotchEdgeInteriorTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchEdgeInteriorTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchEdgeInteriorTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // First recover each cell's Notch and Delta concentrations from the ODEs and store
    // in CellData and CellEdgeData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());

        /* Cell Edges delta notch collection */
        std::vector<double> notch_vec;
        std::vector<double> delta_vec;


        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<DeltaNotchEdgeSrnModel> p_model
                    = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(p_cell_edge_model->GetEdgeSrn(i));
            double edge_delta = p_model->GetDelta();
            double edge_notch = p_model->GetNotch();

            delta_vec.push_back(edge_delta);
            notch_vec.push_back(edge_notch);
        }
        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge delta", delta_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge notch", notch_vec);

        /* Interior delta notch collection */
        //Filling interior delta/notch value, interior model must be specified for this edge-interior example
        assert(p_cell_edge_model->GetInteriorSrn() != nullptr);
        boost::shared_ptr<DeltaNotchInteriorSrnModel> p_interior_model
                = boost::static_pointer_cast<DeltaNotchInteriorSrnModel>(p_cell_edge_model->GetInteriorSrn());

        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellData()->SetItem("interior delta", p_interior_model->GetDelta());
        cell_iter->GetCellData()->SetItem("interior notch", p_interior_model->GetNotch());
    }

    // After the edge data is filled, fill the edge neighbour data
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        const unsigned num_cell_edges = p_cell_edge_model->GetNumEdgeSrn();
        std::vector<double> neigh_mean_delta(num_cell_edges);

        // Cell interior
        double total_edge_delta = 0;
        double total_edge_notch = 0;

        auto edges_delta = cell_iter->GetCellEdgeData()->GetItem("edge delta");
        auto edges_notch = cell_iter->GetCellEdgeData()->GetItem("edge notch");

        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            total_edge_delta += edges_delta[i];
            total_edge_notch += edges_notch[i];
        }

        cell_iter->GetCellData()->SetItem("total neighbour edge delta", total_edge_delta);
        cell_iter->GetCellData()->SetItem("total edge notch", total_edge_notch);

        // Cell edge
        for (unsigned i=0; i<num_cell_edges; ++i)
        {
            // Get neighbouring cell's values of delta on this
            auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(*cell_iter, i);
            double mean_delta = 0;
            for (auto neighbourIndex: elemNeighbours)
            {
                auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                std::vector<double> neighbour_delta_vec = neighbourCell->GetCellEdgeData()->GetItem("edge delta");
                mean_delta += neighbour_delta_vec[neighbourIndex.second];
            }
            if (elemNeighbours.size() > 0)
            {
                mean_delta = mean_delta/elemNeighbours.size();
            }
            neigh_mean_delta[i] = mean_delta;
        }
        cell_iter->GetCellEdgeData()->SetItem("neighbour delta", neigh_mean_delta);
    }
}

template<unsigned DIM>
void DeltaNotchEdgeInteriorTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}


// Explicit instantiation
template class DeltaNotchEdgeInteriorTrackingModifier<1>;
template class DeltaNotchEdgeInteriorTrackingModifier<2>;
template class DeltaNotchEdgeInteriorTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchEdgeInteriorTrackingModifier)
