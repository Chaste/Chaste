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

#include "VertexBasedPopulationSrn.hpp"
#include "CellSrnModel.hpp"

template <unsigned DIM>
VertexBasedPopulationSrn<DIM>::VertexBasedPopulationSrn()
    : mpCellPopulation(nullptr)
{
}

template <unsigned DIM>
VertexBasedPopulationSrn<DIM>::~VertexBasedPopulationSrn()
{
}

template <unsigned DIM>
void VertexBasedPopulationSrn<DIM>::SetVertexCellPopulation(VertexBasedCellPopulation<DIM>* pCellPopulation)
{
    mpCellPopulation = pCellPopulation;
}

template <unsigned DIM>
void VertexBasedPopulationSrn<DIM>::UpdateSrnAfterBirthOrDeath(VertexElementMap& rElementMap)
{
    // Get recorded edge operations
    const std::vector<EdgeOperation>& r_operations = mpCellPopulation->rGetMesh().GetOperationRecorder()->GetEdgeOperations();
    for (auto operation : r_operations)
    {
        /*
         * An operation with deleted element may be recorded, e.g. following a 
         * T2 swap, neighbouring element may become triangular due to node 
         * merging and may subsequently be marked for deletion, despite being 
         * marked to perform edge operations (e.g. node merging).
         */
        switch (operation.GetOperation())
        {
            // All operations, except cell division, require the same information
            case EDGE_OPERATION_ADD:
            case EDGE_OPERATION_NODE_MERGE:
            case EDGE_OPERATION_SPLIT:
            case EDGE_OPERATION_MERGE:
            {
                const unsigned stored_index = operation.GetElementIndex();
                unsigned location_index = stored_index;

                /*
                 * If operation is recorded before element indices are changed. 
                 * For example, if the operations recorded during T2 swap.
                 */
                if (operation.IsElementIndexRemapped())
                {
                    // If the element is deleted, then ignore this operation
                    if (rElementMap.IsDeleted(location_index))
                    {
                        // LCOV_EXCL_START
                        break; // This line is difficult to test
                        // LCOV_EXCL_STOP
                    }
                    else
                    {
                        location_index = rElementMap.GetNewIndex(stored_index);
                    }
                }
                // Get the necessary information to perform Remap
                const EdgeRemapInfo& r_edge_change = operation.rGetRemapInfo();
                CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(location_index);
                auto p_old_model = static_cast<CellSrnModel*>(p_cell->GetSrnModel());
                std::vector<AbstractSrnModelPtr> old_srn_edges = p_old_model->GetEdges();
                RemapCellSrn(old_srn_edges, p_old_model, r_edge_change);
                old_srn_edges.clear();
                break;
            }
            case EDGE_OPERATION_DIVIDE:
            {
                const unsigned location_index_1 = rElementMap.GetNewIndex(operation.GetElementIndex());
                const unsigned location_index_2 = rElementMap.GetNewIndex(operation.GetElementIndex2());

                const EdgeRemapInfo& r_edge_change_1 = operation.rGetRemapInfo();
                const EdgeRemapInfo& r_edge_change_2 = operation.rGetRemapInfo2();

                CellPtr p_cell_1 = mpCellPopulation->GetCellUsingLocationIndex(location_index_1);
                CellPtr p_cell_2 = mpCellPopulation->GetCellUsingLocationIndex(location_index_2);

                auto p_old_model_1 = static_cast<CellSrnModel*>(p_cell_1->GetSrnModel());
                std::vector<AbstractSrnModelPtr> parent_srn_edges = p_old_model_1->GetEdges();
                auto p_old_model_2 = static_cast<CellSrnModel*>(p_cell_2->GetSrnModel());

                RemapCellSrn(parent_srn_edges, p_old_model_1, r_edge_change_1);
                RemapCellSrn(parent_srn_edges, p_old_model_2, r_edge_change_2);
                parent_srn_edges.clear();
                break;
            }
        }
    }
    mpCellPopulation->rGetMesh().GetOperationRecorder()->ClearEdgeOperations();
}

template <unsigned DIM>
void VertexBasedPopulationSrn<DIM>::RemapCellSrn(std::vector<AbstractSrnModelPtr> parentSrnEdges,
                                                 CellSrnModel* pCellSrn,
                                                 const EdgeRemapInfo& rEdgeChange)
{
    auto edge_mapping = rEdgeChange.GetEdgesMapping();
    std::vector<AbstractSrnModelPtr> new_edge_srn(edge_mapping.size());
    const std::vector<double> split_proportions = rEdgeChange.GetSplitProportions();
    const unsigned num_edges = edge_mapping.size();
    std::vector<unsigned> shrunk_edges;
  
    /*
     * Go through the edges, check its status and the index corresponding to the 
     * edge status before rearrangement. Go through the SRN model.
     */
    for (unsigned i = 0; i < num_edges; i++)
    {
        // The remap_index, if +ve refers to the SRN index of the oldModel, if -ve then it's a new edge
        const long remap_index = edge_mapping[i];

        /*
         * remap_status can be the following:
         * 0 - Direct remapping, the edge SRN of the oldModel can be transferred 
         *         directly to the new model
         * 1 - The edge is a split point between the dividing cells
         * 2 - This is a new edge i.e. the dividing line in the middle of the 
         *         old and new cells
         * 3 - Edge below or above the edge that was deleted due to node merging
         * 4 - Edge above was merged into this edge
         */
        const unsigned remap_status = rEdgeChange.GetEdgesStatus()[i];

        // Remap index cannot be negative when it's a direct remap or an edge split
        assert(!((remap_status == 0 || remap_status == 1) && remap_index < 0));
        switch(remap_status)
        {
            // Edge SRN remains the same
            case 0:
            {
                new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parentSrnEdges[remap_index]->CreateSrnModel());
                break;
            }
            case 1:
            {
                /* 
                 * Edge split depends on the relative splitting node position 
                 * because currently we assume that edge concentrations are 
                 * uniform on the edge.
                 */
                new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parentSrnEdges[remap_index]->CreateSrnModel());
                assert(split_proportions[i] >= 0);
                new_edge_srn[i]->SplitEdgeSrn(split_proportions[i]);
                break;
            }
            case 2:
            {
                // The edge is new
                new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parentSrnEdges[0]->CreateSrnModel());
                new_edge_srn[i]->SetCell(pCellSrn->GetCell());
                new_edge_srn[i]->InitialiseDaughterCell();
                break;
            }
            case 3:
            {
                // If the edge above or below this edge was deleted
                new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parentSrnEdges[remap_index]->CreateSrnModel());

                const bool is_prev_edge = rEdgeChange.GetEdgesStatus()[(i+1)%num_edges]==3;

                // Find the shrunk edge
                unsigned shrunk_edge = 0;

                // If this edge is below the shrunk edge
                if (is_prev_edge)
                {
                    shrunk_edge = (remap_index+1)%(num_edges+1);
                    shrunk_edges.push_back(shrunk_edge);
                }
                else
                {
                    shrunk_edge = remap_index==0 ? num_edges : remap_index-1;
                }
                new_edge_srn[i]->AddShrunkEdgeSrn(parentSrnEdges[shrunk_edge].get());
                break;
            }
            case 4:
            {
                // Add SRNs from the edge above to this edge
                new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parentSrnEdges[remap_index]->CreateSrnModel());
                const unsigned next_edge_index = (remap_index+1)%(num_edges+1);
                new_edge_srn[i]->AddMergedEdgeSrn(parentSrnEdges[next_edge_index].get());
                break;
            }
        }
        // Setting the new local edge index and the cell
        new_edge_srn[i]->SetEdgeLocalIndex(i);
        new_edge_srn[i]->SetCell(pCellSrn->GetCell());
    }
  
    // For the case when edge quantities are returned into interior when edge shrinks due to node merging
    auto p_interior_srn = boost::shared_ptr<AbstractSrnModel>(pCellSrn->GetInteriorSrn());
    if (p_interior_srn != nullptr)
    {
        for (unsigned shrunk_edge_index : shrunk_edges)
        {
            p_interior_srn->AddShrunkEdgeToInterior(parentSrnEdges[shrunk_edge_index].get());
        }
    }

    pCellSrn->AddEdgeSrn(new_edge_srn);
    assert(num_edges == pCellSrn->GetNumEdgeSrn());
}

template class VertexBasedPopulationSrn<1>;
template class VertexBasedPopulationSrn<2>;
template class VertexBasedPopulationSrn<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedPopulationSrn)