/*

Copyright (c) 2005-2019, University of Oxford.
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

    // Update the cell
    this->UpdateCellData(rCellPopulation);
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

        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
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
        for (unsigned i = 0; i < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            double mean_delta = 0.0;
            double delta_counter = 0;

            // Edge neighbour outside of cell
            auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(*cell_iter, i);
            for (auto neighbourIndex: elemNeighbours)
            {
                auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);

                std::vector<double> neighbour_delta_vec = neighbourCell->GetCellEdgeData()->GetItem("delta");

                mean_delta += neighbour_delta_vec[neighbourIndex.second];
                delta_counter++;
            }

            // Edge neighbour inside of cell
            // In this case we're only looking at immediate neighbour to the current edge
            auto prevNeighbour = ((int)i -1);
            if (prevNeighbour < 0)
            {
                prevNeighbour = p_cell_edge_model->GetNumEdgeSrn() - 1;
            }
            auto nextNeighbour = (i +1) % p_cell_edge_model->GetNumEdgeSrn();
            delta_counter += 2;
            mean_delta += delta_vec[prevNeighbour];
            mean_delta += delta_vec[nextNeighbour];

            mean_delta = mean_delta/delta_counter;

            // Store the delta
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
AbstractSrnModel *CellEdgeDeltaNotchTrackingModifier<DIM>::CreateEmptySrnEdgeModel()
{
    return new DeltaNotchEdgeSrnModel();
}

template<unsigned int DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::EdgeAdded(AbstractCellPopulation<DIM,DIM> &rCellPopulation,
                                                        unsigned locationIndex, unsigned edgeLocalIndex, AbstractSrnModelPtr addedEdge)
{
    auto deltaNotchNewSrnEdge = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(addedEdge);

    // New edges have a concentration of 0
    deltaNotchNewSrnEdge->SetDelta(0);
    deltaNotchNewSrnEdge->SetNotch(0);
}

template<unsigned int DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::EdgeRemoved(AbstractCellPopulation<DIM,DIM> &rCellPopulation,
                                                     unsigned locationIndex, unsigned edgeLocalIndex, AbstractSrnModelPtr oldSrnEdge)
{
}

template<unsigned int DIM>
void CellEdgeDeltaNotchTrackingModifier<DIM>::EdgeDivide(AbstractSrnModelPtr oldSrnEdge, AbstractSrnModelPtr newSrnEdge)
{
    // Convert to delta notch SRN type
    auto deltaNotchOldSrnEdge = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(oldSrnEdge);
    auto deltaNotchNewSrnEdge = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(newSrnEdge);

    // In this example we're just halving the concentrations
    deltaNotchNewSrnEdge->SetDelta(0.5*deltaNotchOldSrnEdge->GetDelta());
    deltaNotchNewSrnEdge->SetNotch(0.5*deltaNotchOldSrnEdge->GetNotch());
}

// Explicit instantiation
template class DeltaNotchCellEdgeTrackingModifier<1>;
template class DeltaNotchCellEdgeTrackingModifier<2>;
template class DeltaNotchCellEdgeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
#include "AbstractCellEdgeBasedSimulationModifier.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchCellEdgeTrackingModifier)
