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

#include "DeltaNotchEdgeInteriorTrackingModifier.hpp"
#include "SrnCellModel.hpp"
#include "DeltaNotchSrnInteriorModel.hpp"
#include "DeltaNotchSrnEdgeModel.hpp"

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

    // Make sure the cell population is updated
    rCellPopulation.Update();

    // First recover each cell's Notch and Delta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<SrnCellModel*>(cell_iter->GetSrnModel());
        boost::shared_ptr<DeltaNotchSrnInteriorModel> p_interior_model
                    = boost::static_pointer_cast<DeltaNotchSrnInteriorModel>(p_cell_edge_model->GetInteriorSrn());

        const double this_notch = p_interior_model->GetNotch();
        const double this_delta = p_interior_model->GetDelta();
        double total_edge_delta = 0;
        double total_edge_notch = 0;
        const unsigned int n_cell_edges = p_cell_edge_model->GetNumEdgeSrn();
        std::vector<double> edge_delta(n_cell_edges);
        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<DeltaNotchSrnEdgeModel> p_model
            = boost::static_pointer_cast<DeltaNotchSrnEdgeModel>(p_cell_edge_model->GetEdgeSrn(i));
            total_edge_delta += p_model->GetNeighbouringDelta();
            total_edge_notch += p_model->GetNotch();
        }
        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellData()->SetItem("interior notch", this_notch);
        cell_iter->GetCellData()->SetItem("interior delta", this_delta);
        cell_iter->GetCellData()->SetItem("total neighbour edge delta", total_edge_delta);
        cell_iter->GetCellData()->SetItem("total edge notch", total_edge_notch);
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
