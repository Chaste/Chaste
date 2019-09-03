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
            break;
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
            //
            /*MutableVertexMesh<DIM, DIM>* pMesh
                    = static_cast<MutableVertexMesh<DIM, DIM>* >(&(mpCellPopulation->rGetMesh()));
            auto element = pMesh->GetElement(location_index_1);
            unsigned int el_edges = element->GetNumEdges();
            std::cout<<"Cell id divide: "<<cell_1->GetCellId()
                                     <<" elem edges: "<<el_edges<<std::endl;
            auto element_2 = pMesh->GetElement(location_index_2);
            el_edges = element_2->GetNumEdges();
            std::cout<<"Cell id divide: "<<cell_2->GetCellId()
                     <<" elem edges: "<<el_edges<<std::endl;
            unsigned int other_index;
            if (location_index_1==0)
                other_index = 1;
            else
                other_index = 0;
            auto element_3 = pMesh->GetElement(other_index);
            el_edges = element_3->GetNumEdges();
            std::cout<<"Cell id non-divided: "<<other_index
                                         <<" elem edges: "<<el_edges<<std::endl;*/
            //
            auto old_model_1 = static_cast<SrnCellModel*>(cell_1->GetSrnModel());
            std::vector<AbstractSrnModelPtr> parent_srn_edges = old_model_1->GetEdges();
            //checking if the new cell has empty edge srns
            auto old_model_2 = static_cast<SrnCellModel*>(cell_2->GetSrnModel());
            assert(old_model_2->GetEdges().empty());
            RemapCellSrn(parent_srn_edges, old_model_1, pEdgeChange_1);
            RemapCellSrn(parent_srn_edges, old_model_2, pEdgeChange_2);
            parent_srn_edges.clear();
            break;
        }
        case EDGE_OPERATION_SPLIT:
        {
            const unsigned int location_index = operation->GetElementIndex();
            EdgeRemapInfo *pEdgeChange = operation->GetNewEdges();
            CellPtr cell = mpCellPopulation->GetCellUsingLocationIndex(location_index);
            auto old_model = static_cast<SrnCellModel*>(cell->GetSrnModel());
            std::vector<AbstractSrnModelPtr> old_srn_edges = old_model->GetEdges();
            for (unsigned int i=0; i<pEdgeChange->GetEdgesMapping().size(); ++i)
            {
                auto remapIndex = pEdgeChange->GetEdgesMapping()[i];
                auto remapStatus = pEdgeChange->GetEdgesStatus()[i];
                if ((remapStatus == 0 || remapStatus == 1) && remapIndex < 0)
                {
                    EXCEPTION("Remap index cannot be negative when it's a direct remap or an edge split");
                }
                if (remapStatus == 2)
                {
                    EXCEPTION("A blank new edge has been created. SPLIT operation recorded wrongly");
                }
            }
            RemapCellSrn(old_srn_edges, old_model, pEdgeChange);

            old_srn_edges.clear();
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
    // Goes through the SRN model
    for (unsigned i = 0; i < edge_mapping.size(); i++)
    {
        //The remapIndex, if +ve refers to the SRN index of the oldModel, if -ve then it's a new edge
        const long int remapIndex = edge_mapping[i];

        //remapStatus can be the following:
        //0 - Direct remapping, the edge SRN of the oldModel can be transferred directly to the new model
        //1 - The edge is a split point between the diving cells, in this example we divide all concentration in half
        //2 - This is a new edge i.e. the dividing line in the middle of the old and new cells
        const unsigned int remapStatus = pEdgeChange->GetEdgesStatus()[i];
        if ((remapStatus == 0 || remapStatus == 1) && remapIndex < 0)
        {
            EXCEPTION("Remap index cannot be negative when it's a direct remap or an edge split");
        }

        if (remapStatus ==0)
        {
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[remapIndex]->CreateSrnModel());
        }
        if (remapStatus == 1)
        {
            //Edge split depends on the relative splitting node position because
            //currently we assume that edge concentrations are uniform on the edge
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[remapIndex]->CreateSrnModel());
            new_edge_srn[i]->ScaleSrnVariables(split_proportions[i]);
        }
        if (remapStatus == 2)
        {
            new_edge_srn[i] = boost::shared_ptr<AbstractSrnModel>(parent_srn_edges[0]->CreateSrnModel());
        }

        //Setting the new local edge index and the cell
        new_edge_srn[i]->SetEdgeLocalIndex(i);
        new_edge_srn[i]->SetCell(pSrnCell->GetCell());
        if (remapStatus == 2)
            new_edge_srn[i]->InitialiseDaughterCell();
    }
    pSrnCell->AddEdgeSrn(new_edge_srn);
}

template class VertexBasedPopulationSrn<1>;
template class VertexBasedPopulationSrn<2>;
template class VertexBasedPopulationSrn<3>;
