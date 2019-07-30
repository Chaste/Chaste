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

#include "NodeBasedCellPopulationWithParticles.hpp"
#include "VtkMeshWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"

template<unsigned DIM>
NodeBasedCellPopulationWithParticles<DIM>::NodeBasedCellPopulationWithParticles(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh)
    : NodeBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh, false)
{
    EXCEPT_IF_NOT(PetscTools::IsSequential());
    if (!locationIndices.empty())
    {
        // Create a set of node indices corresponding to particles
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices;
        std::set<unsigned> particle_indices;

        for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
             node_iter != rMesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            node_indices.insert(node_iter->GetIndex());
        }
        for (unsigned i=0; i<locationIndices.size(); i++)
        {
            location_indices.insert(locationIndices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices.begin(), location_indices.end(),
                            std::inserter(particle_indices, particle_indices.begin()));

        // This method finishes and then calls Validate()
        SetParticles(particle_indices);
    }
    else
    {
        for (typename NodesOnlyMesh<DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
             node_iter != rMesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            (*node_iter).SetIsParticle(false);
        }
        NodeBasedCellPopulationWithParticles::Validate();
    }
}

template<unsigned DIM>
NodeBasedCellPopulationWithParticles<DIM>::NodeBasedCellPopulationWithParticles(NodesOnlyMesh<DIM>& rMesh)
    : NodeBasedCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
bool NodeBasedCellPopulationWithParticles<DIM>::IsParticle(unsigned index)
{
    return this->GetNode(index)->IsParticle();
}

template<unsigned DIM>
std::set<unsigned> NodeBasedCellPopulationWithParticles<DIM>::GetParticleIndices()
{
    std::set<unsigned> particle_indices;

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        if (node_iter->IsParticle())
        {
            particle_indices.insert(node_iter->GetIndex());
        }
    }

    return particle_indices;
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::SetParticles(const std::set<unsigned>& rParticleIndices)
{
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        if (rParticleIndices.find(node_iter->GetIndex()) != rParticleIndices.end())
        {
            node_iter->SetIsParticle(true);
        }
    }
    NodeBasedCellPopulationWithParticles::Validate();
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::UpdateParticlesAfterReMesh(NodeMap& rMap)
{
}

template<unsigned DIM>
CellPtr NodeBasedCellPopulationWithParticles<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    assert(pNewCell);

    // Add new cell to population
    CellPtr p_created_cell = AbstractCentreBasedCellPopulation<DIM>::AddCell(pNewCell, pParentCell);
    assert(p_created_cell == pNewCell);

    // Then set the new cell radius in the NodesOnlyMesh
    unsigned node_index = this->GetLocationIndexUsingCell(p_created_cell);
    this->GetNode(node_index)->SetRadius(0.5);
    this->GetNode(node_index)->SetIsParticle(false);

    // Return pointer to new cell
    return p_created_cell;
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::Validate()
{
    std::map<unsigned, bool> validated_nodes;
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        validated_nodes[node_iter->GetIndex()] = node_iter->IsParticle();
    }

    // Look through all of the cells and record what node they are associated with.
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell((*cell_iter));

        // If the node attached to this cell is labelled as a particle, then throw an error
        if (this->GetNode(node_index)->IsParticle())
        {
            EXCEPTION("Node " << node_index << " is labelled as a particle and has a cell attached");
        }
        validated_nodes[node_index] = true;
    }

    for (std::map<unsigned, bool>::iterator map_iter = validated_nodes.begin();
         map_iter != validated_nodes.end();
         map_iter++)
    {
        if (!map_iter->second)
        {
            EXCEPTION("Node " << map_iter->first << " does not appear to be a particle or has a cell associated with it");
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::AcceptCellWritersAcrossPopulation()
{
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // If it isn't a particle then there might be cell writers attached
        if (! this->IsParticle(node_iter->GetIndex()))
        {
            for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
                 cell_writer_iter != this->mCellWriters.end();
                 ++cell_writer_iter)
            {
                CellPtr cell_from_node = this->GetCellUsingLocationIndex(node_iter->GetIndex());
                this->AcceptCellWriter(*cell_writer_iter, cell_from_node);
            }
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    // Store the present time as a string
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();

    // Make sure the nodes are ordered contiguously in memory
    NodeMap map(1 + this->mpNodesOnlyMesh->GetMaximumNodeIndex());
    this->mpNodesOnlyMesh->ReMesh(map);

    // Store the number of cells for which to output data to VTK
    unsigned num_nodes = this->GetNumNodes();
    std::vector<double> rank(num_nodes);
    std::vector<double> particles(num_nodes);

    unsigned num_cell_data_items = 0;
    std::vector<std::string> cell_data_names;

    // We assume that the first cell is representative of all cells
    if (num_nodes > 0)
    {
        num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
        cell_data_names = this->Begin()->GetCellData()->GetKeys();
    }

    std::vector<std::vector<double> > cell_data;
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cell_data_var(num_nodes);
        cell_data.push_back(cell_data_var);
    }

    // Create mesh writer for VTK output
    VtkMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results_"+time.str(), false);
    mesh_writer.SetParallelFiles(*(this->mpNodesOnlyMesh));

    // Iterate over any cell writers that are present
    for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
         cell_writer_iter != this->mCellWriters.end();
         ++cell_writer_iter)
    {
        // Create vector to store VTK cell data
        std::vector<double> vtk_cell_data(num_nodes);

        // Loop over nodes
        for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
             node_iter != this->mrMesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();

            // If this node is a particle (not a cell), then we set the 'dummy' VTK cell data for this to be -2.0...
            if (this->IsParticle(node_index))
            {
                vtk_cell_data[node_index] = -2.0;
            }
            else
            {
                // ...otherwise we populate the vector of VTK cell data as usual
                CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);
                vtk_cell_data[node_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
            }
        }

        mesh_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    }

    // Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned global_index = this->GetLocationIndexUsingCell(*cell_iter);
        unsigned node_index = this->rGetMesh().SolveNodeMapping(global_index);

        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cell_data[var][node_index] = cell_iter->GetCellData()->GetItem(cell_data_names[var]);
        }

        rank[node_index] = (PetscTools::GetMyRank());
    }

    mesh_writer.AddPointData("Process rank", rank);

    // Loop over nodes
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        particles[node_index] = (double) (this->IsParticle(node_index));
    }

    mesh_writer.AddPointData("Non-particles", particles);

    if (num_cell_data_items > 0)
    {
        for (unsigned var=0; var<cell_data.size(); var++)
        {
            mesh_writer.AddPointData(cell_data_names[var], cell_data[var]);
        }
    }

    mesh_writer.WriteFilesUsingMesh(*(this->mpNodesOnlyMesh));

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    EXCEPT_IF_NOT(PetscTools::IsSequential());
    {
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
/*    {
        // Parallel vtu files  .vtu -> .pvtu
        *(this->mpVtkMetaFile) << ".pvtu\"/>\n";
    }*/
#endif //CHASTE_VTK
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
template class NodeBasedCellPopulationWithParticles<1>;
template class NodeBasedCellPopulationWithParticles<2>;
template class NodeBasedCellPopulationWithParticles<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithParticles)
