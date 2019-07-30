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

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellLocationIndexWriter.hpp"

template<unsigned DIM>
MeshBasedCellPopulationWithGhostNodes<DIM>::MeshBasedCellPopulationWithGhostNodes(
     MutableMesh<DIM, DIM>& rMesh,
     std::vector<CellPtr>& rCells,
     const std::vector<unsigned> locationIndices,
     bool deleteMesh,
     double ghostSpringStiffness)
             : MeshBasedCellPopulation<DIM,DIM>(rMesh, rCells, locationIndices, deleteMesh, false), // do not call the base class Validate()
               mGhostSpringStiffness(ghostSpringStiffness)
{
    if (!locationIndices.empty())
    {
        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<this->GetNumNodes(); i++)
        {
            node_indices.insert(this->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<locationIndices.size(); i++)
        {
            location_indices.insert(locationIndices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices.begin(), location_indices.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        // This method finishes and then calls Validate()
        SetGhostNodes(ghost_node_indices);
    }
    else
    {
        this->mIsGhostNode = std::vector<bool>(this->GetNumNodes(), false);
        Validate();
    }
}

template<unsigned DIM>
MeshBasedCellPopulationWithGhostNodes<DIM>::MeshBasedCellPopulationWithGhostNodes(MutableMesh<DIM, DIM>& rMesh,
                                                                                  double ghostSpringStiffness)
    : MeshBasedCellPopulation<DIM,DIM>(rMesh),
      mGhostSpringStiffness(ghostSpringStiffness)
{
}

template<unsigned DIM>
MeshBasedCellPopulationWithGhostNodes<DIM>::~MeshBasedCellPopulationWithGhostNodes()
{
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* MeshBasedCellPopulationWithGhostNodes<DIM>::GetTetrahedralMeshForPdeModifier()
{
    EXCEPTION("Currently can't solve PDEs on meshes with ghost nodes");
    return static_cast<TetrahedralMesh<DIM, DIM>*>(&(this->mrMesh));
}

template<unsigned DIM>
std::vector<bool>& MeshBasedCellPopulationWithGhostNodes<DIM>::rGetGhostNodes()
{
    return this->mIsGhostNode;
}

template<unsigned DIM>
bool MeshBasedCellPopulationWithGhostNodes<DIM>::IsGhostNode(unsigned index)
{
    return this->mIsGhostNode[index];
}

template<unsigned DIM>
std::set<unsigned> MeshBasedCellPopulationWithGhostNodes<DIM>::GetGhostNodeIndices()
{
    std::set<unsigned> ghost_node_indices;
    for (unsigned i=0; i<this->mIsGhostNode.size(); i++)
    {
        if (this->mIsGhostNode[i])
        {
            ghost_node_indices.insert(i);
        }
    }
    return ghost_node_indices;
}

template<unsigned DIM>
void MeshBasedCellPopulationWithGhostNodes<DIM>::SetGhostNodes(const std::set<unsigned>& rGhostNodeIndices)
{
    // Reinitialise all entries of mIsGhostNode to false
    this->mIsGhostNode = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsGhostNode
    for (std::set<unsigned>::iterator iter=rGhostNodeIndices.begin(); iter!=rGhostNodeIndices.end(); ++iter)
    {
        this->mIsGhostNode[*iter] = true;
    }

    Validate();
}

template<unsigned DIM>
c_vector<double, DIM> MeshBasedCellPopulationWithGhostNodes<DIM>::CalculateForceBetweenGhostNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    assert(rNodeAGlobalIndex != rNodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    const c_vector<double, DIM>& r_node_a_location = this->GetNode(rNodeAGlobalIndex)->rGetLocation();
    const c_vector<double, DIM>& r_node_b_location = this->GetNode(rNodeBGlobalIndex)->rGetLocation();

    // There is reason not to subtract one position from the other (cylindrical meshes)
    unit_difference = this->mrMesh.GetVectorFromAtoB(r_node_a_location, r_node_b_location);

    double distance_between_nodes = norm_2(unit_difference);
    unit_difference /= distance_between_nodes;

    double rest_length = 1.0;

    return mGhostSpringStiffness * unit_difference * (distance_between_nodes - rest_length);
}

template<unsigned DIM>
CellPtr MeshBasedCellPopulationWithGhostNodes<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    // Add new cell to population
    CellPtr p_created_cell = MeshBasedCellPopulation<DIM,DIM>::AddCell(pNewCell, pParentCell);
    assert(p_created_cell == pNewCell);

    // Update size of mIsGhostNode if necessary
    unsigned new_node_index = this->GetLocationIndexUsingCell(p_created_cell);

    if (this->GetNumNodes() > this->mIsGhostNode.size())
    {
        this->mIsGhostNode.resize(this->GetNumNodes());
        this->mIsGhostNode[new_node_index] = false;
    }

    // Return pointer to new cell
    return p_created_cell;
}

template<unsigned DIM>
void MeshBasedCellPopulationWithGhostNodes<DIM>::Validate()
{
    // Get a list of all the nodes that are ghosts
    std::vector<bool> validated_node = mIsGhostNode;
    assert(mIsGhostNode.size()==this->GetNumNodes());

    // Look through all of the cells and record what node they are associated with.
    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell((*cell_iter));

        // If the node attached to this cell is labelled as a ghost node, then throw an error
        if (mIsGhostNode[node_index])
        {
            EXCEPTION("Node " << node_index << " is labelled as a ghost node and has a cell attached");
        }
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            EXCEPTION("Node " << i << " does not appear to be a ghost node or have a cell associated with it");
        }
    }
}

template<unsigned DIM>
void MeshBasedCellPopulationWithGhostNodes<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{
    // Copy mIsGhostNode to a temporary vector
    std::vector<bool> ghost_nodes_before_remesh = mIsGhostNode;

    // Reinitialise mIsGhostNode
    mIsGhostNode.clear();
    mIsGhostNode.resize(this->GetNumNodes());

    // Update mIsGhostNode using the node map
    for (unsigned old_index=0; old_index<rMap.GetSize(); old_index++)
    {
        if (!rMap.IsDeleted(old_index))
        {
            unsigned new_index = rMap.GetNewIndex(old_index);
            mIsGhostNode[new_index] = ghost_nodes_before_remesh[old_index];
        }
    }
}

template<unsigned DIM>
std::set<unsigned> MeshBasedCellPopulationWithGhostNodes<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    unsigned node_index = this->GetLocationIndexUsingCell(pCell);
    std::set<unsigned> neighbour_indices = this->GetNeighbouringNodeIndices(node_index);

    // Remove ghost nodes from the neighbour indices
    for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
         iter != neighbour_indices.end();)
    {
        if (this->IsGhostNode(*iter))
        {
            neighbour_indices.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }

    return neighbour_indices;
}

template<unsigned DIM>
void MeshBasedCellPopulationWithGhostNodes<DIM>::AcceptCellWritersAcrossPopulation()
{
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // If it isn't a ghost node then there might be cell writers attached
        if (! this->IsGhostNode(node_iter->GetIndex()))
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
void MeshBasedCellPopulationWithGhostNodes<DIM>::ApplyGhostForces()
{
    // Initialise vector of forces on ghost nodes
    std::vector<c_vector<double, DIM> > drdt(this->GetNumNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i] = zero_vector<double>(DIM);
    }

    // Calculate forces on ghost nodes
    for (typename MutableMesh<DIM, DIM>::EdgeIterator edge_iterator = static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).EdgesBegin();
        edge_iterator != static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).EdgesEnd();
        ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenGhostNodes(nodeA_global_index, nodeB_global_index);

        if (!this->mIsGhostNode[nodeA_global_index])
        {
            drdt[nodeB_global_index] -= force;
        }
        else
        {
            drdt[nodeA_global_index] += force;

            if (this->mIsGhostNode[nodeB_global_index])
            {
                drdt[nodeB_global_index] -= force;
            }
        }
    }

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        if (this->mIsGhostNode[node_index])
        {
            node_iter->ClearAppliedForce();
            node_iter->AddAppliedForceContribution(drdt[node_index]);
        }
    }
}

template<unsigned DIM>
void MeshBasedCellPopulationWithGhostNodes<DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellLocationIndexWriter>())
        {
            this-> template AddCellWriter<CellLocationIndexWriter>();
        }
    }

    MeshBasedCellPopulation<DIM, DIM>::OpenWritersFiles(rOutputFileHandler);
}

template<unsigned DIM>
void MeshBasedCellPopulationWithGhostNodes<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    if (this->mpVoronoiTessellation != nullptr)
    {
        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
        std::stringstream time;
        time << num_timesteps;

        // Create mesh writer for VTK output
        VertexMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results", false);

        // Iterate over any cell writers that are present
        unsigned num_vtk_cells = this->mpVoronoiTessellation->GetNumElements();
        for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
             cell_writer_iter != this->mCellWriters.end();
             ++cell_writer_iter)
        {
            // Create vector to store VTK cell data
            std::vector<double> vtk_cell_data(num_vtk_cells);

            // Loop over elements of mpVoronoiTessellation
            for (typename VertexMesh<DIM, DIM>::VertexElementIterator elem_iter = this->mpVoronoiTessellation->GetElementIteratorBegin();
                 elem_iter != this->mpVoronoiTessellation->GetElementIteratorEnd();
                 ++elem_iter)
            {
                // Get the indices of this element and the corresponding node in mrMesh
                unsigned elem_index = elem_iter->GetIndex();
                unsigned node_index = this->mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

                // If this node corresponds to a ghost node, set any "cell" data to be -1.0
                if (this->IsGhostNode(node_index))
                {
                    // Populate the vector of VTK cell data
                    vtk_cell_data[elem_index] = -1.0;
                }
                else
                {
                    // Get the cell corresponding to this node
                    CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

                    // Populate the vector of VTK cell data
                    vtk_cell_data[elem_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
                }
            }

            mesh_writer.AddCellData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
        }

        // Next, record which nodes are ghost nodes
        // Note that the cell writer hierarchy can not be used to do this as ghost nodes don't have corresponding cells.
        std::vector<double> ghosts(num_vtk_cells);
        for (typename VertexMesh<DIM, DIM>::VertexElementIterator elem_iter = this->mpVoronoiTessellation->GetElementIteratorBegin();
             elem_iter != this->mpVoronoiTessellation->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            unsigned node_index = this->mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);
            ghosts[elem_index]  = (double) (this->IsGhostNode(node_index));
        }
        mesh_writer.AddCellData("Non-ghosts", ghosts);

        ///\todo #1975 - deal with possibility of information stored in CellData

        mesh_writer.WriteVtkUsingMesh(*(this->mpVoronoiTessellation), time.str());
        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
#endif //CHASTE_VTK
}

template<unsigned DIM>
void MeshBasedCellPopulationWithGhostNodes<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<GhostSpringStiffness>" << mGhostSpringStiffness << "</GhostSpringStiffness>\n";

    // Call method on direct parent class
    MeshBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
template class MeshBasedCellPopulationWithGhostNodes<1>;
template class MeshBasedCellPopulationWithGhostNodes<2>;
template class MeshBasedCellPopulationWithGhostNodes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedCellPopulationWithGhostNodes)
