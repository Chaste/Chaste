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
#include "SrnCellModel.hpp"
#include "VertexBasedPopulationSrn.hpp"
template <unsigned DIM>
VertexBasedPopulationSrn<DIM>::VertexBasedPopulationSrn()
{}

template <unsigned DIM>
void VertexBasedPopulationSrn<DIM>::SetVertexCellPopulation(VertexBasedCellPopulation<DIM>* p_vertex_population)
{
    mpCellPopulation = p_vertex_population;
}

template <unsigned DIM>
void VertexBasedPopulationSrn<DIM>::UpdateSrnAfterBirthOrDeath()
{
    for (auto operation:mpCellPopulation->GetCellEdgeChangeOperations())
    {
        switch (operation->GetOperation())
        {
        case EDGE_OPERATION_ADD:
        case EDGE_OPERATION_NODE_MERGE:
        case EDGE_OPERATION_SPLIT:
        case EDGE_OPERATION_MERGE:
        case EDGE_OPERATION_NEW_NEIGHBOUR:
        {
            const unsigned int location_index = operation->GetElementIndex();
            EdgeRemapInfo *pEdgeChange = operation->GetNewEdges();
            CellPtr cell = mpCellPopulation->GetCellUsingLocationIndex(location_index);
            auto old_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            std::vector<AbstractSrnModelPtr> old_srn_edges = old_model->GetEdges();
            RemapCellSrn(old_srn_edges, old_model, pEdgeChange);
            old_srn_edges.clear();
            break;
        }
        case EDGE_OPERATION_DELETE:
            break;
        case EDGE_OPERATION_DIVIDE:
        {
            const unsigned int location_index_1 = operation->GetElementIndex();
            const unsigned int location_index_2 = operation->GetElementIndex2();
            EdgeRemapInfo *pEdgeChange_1 = operation->GetNewEdges();
            EdgeRemapInfo *pEdgeChange_2 = operation->GetNewEdges2();
            CellPtr cell_1 = mpCellPopulation->GetCellUsingLocationIndex(location_index_1);
            CellPtr cell_2 = mpCellPopulation->GetCellUsingLocationIndex(location_index_2);

            auto old_model_1 = static_cast<SrnCellModel*>(cell_1->GetSrnModel());
            std::vector<AbstractSrnModelPtr> parent_srn_edges = old_model_1->GetEdges();
            auto old_model_2 = static_cast<SrnCellModel*>(cell_2->GetSrnModel());

            RemapCellSrn(parent_srn_edges, old_model_1, pEdgeChange_1);
            RemapCellSrn(parent_srn_edges, old_model_2, pEdgeChange_2);
            parent_srn_edges.clear();
            break;
        }
        }
    }
    mpCellPopulation->ClearCellEdgeOperations();
}

template <unsigned DIM>
void VertexBasedPopulationSrn<DIM>::RemapCellSrn(std::vector<AbstractSrnModelPtr> parent_srn_edges,
                                                 SrnCellModel* pSrnCell,
                                                 EdgeRemapInfo* pEdgeChange)
{
    auto edge_mapping = pEdgeChange->GetEdgesMapping();
    std::vector<AbstractSrnModelPtr> new_edge_srn(edge_mapping.size());
    const std::vector<double> split_proportions = pEdgeChange->GetSplitProportions();
    const unsigned int n_edges = edge_mapping.size();
    std::vector<unsigned int> shrunk_edges;
    // Go through the SRN model
    for (unsigned i = 0; i < n_edges; i++)
    {
        //The remapIndex, if +ve refers to the SRN index of the oldModel, if -ve then it's a new edge
        const long int remapIndex = edge_mapping[i];

        //remapStatus can be the following:
        //0 - Direct remapping, the edge SRN of the oldModel can be transferred directly to the new model
        //1 - The edge is a split point between the diving cells, in this example we divide all concentration in half
        //2 - This is a new edge i.e. the dividing line in the middle of the old and new cells
        //3 - Edge below or above an edge that was deleted due to node merging
        //4 - Edge above was merged into this edge
        //5 - Edge acquired a new neighbour. For example after a T2 swap
        const unsigned int remapStatus = pEdgeChange->GetEdgesStatus()[i];
        if ((remapStatus == 0 || remapStatus == 1) && remapIndex < 0)
        {
            EXCEPTION("Remap index cannot be negative when it's a direct remap or an edge split");
        }
        switch(remapStatus)
        {
        case 0:
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[remapIndex]->CreateSrnModel());
            break;
        case 1:
            //Edge split depends on the relative splitting node position because
            //currently we assume that edge concentrations are uniform on the edge
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[remapIndex]->CreateSrnModel());
            assert(split_proportions[i]>=0);
            new_edge_srn[i]->SplitEdgeSrn(split_proportions[i]);
            break;
        case 2:
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[0]->CreateSrnModel());
            break;
        case 3:
        {
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[remapIndex]->CreateSrnModel());
            const bool isPrevEdge = pEdgeChange->GetEdgesStatus()[(i+1)%n_edges]==3;
            //Find the shrunk edge
            unsigned int shrunkEdge = 0;
            if (isPrevEdge)
            {
                shrunkEdge = (remapIndex+1)%(n_edges+1);
                shrunk_edges.push_back(shrunkEdge);
            }
            else
            {
                shrunkEdge = remapIndex-1;
                if (shrunkEdge<0)
                {
                    shrunkEdge = n_edges;
                }
            }
            new_edge_srn[i]->AddShrunkEdgeSrn(parent_srn_edges[shrunkEdge].get());
            break;
        }
        case 4:
        {
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[remapIndex]->CreateSrnModel());
            const unsigned int next_edge_index = (remapIndex+1)%(n_edges+1);
            new_edge_srn[i] -> AddMergedEdgeSrn(parent_srn_edges[next_edge_index].get());
            break;
        }
        case 5:
        {
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[remapIndex]->CreateSrnModel());
            new_edge_srn[i]->UpdateSrnAfterNewNeighbour();
            break;
        }
        }

        //Setting the new local edge index and the cell
        new_edge_srn[i]->SetEdgeLocalIndex(i);
        new_edge_srn[i]->SetCell(pSrnCell->GetCell());
        if (remapStatus == 2)
            new_edge_srn[i]->InitialiseDaughterCell();
    }

    //For the case when edge quantities are returned into interior when edge shrinks due to node merging
    boost::shared_ptr<AbstractSrnModel> interior_srn
    =boost::shared_ptr<AbstractSrnModel>(pSrnCell->GetInteriorSrn());
    if (interior_srn != nullptr)
    {
        for (unsigned int shrunkEdgeIndex:shrunk_edges)
        {
            interior_srn->AddShrunkEdgeToInterior(parent_srn_edges[shrunkEdgeIndex].get());
        }
        interior_srn->UpdateSrnAfterNewNeighbour();
    }

    pSrnCell->AddEdgeSrn(new_edge_srn);
}

template class VertexBasedPopulationSrn<1>;
template class VertexBasedPopulationSrn<2>;
template class VertexBasedPopulationSrn<3>;
