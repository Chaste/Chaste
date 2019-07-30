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

#include "NodeBasedCellPopulation.hpp"
#include "MathsCustomFunctions.hpp"
#include "VtkMeshWriter.hpp"

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh,
                                      bool validate)
    : AbstractCentreBasedCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDeleteMesh(deleteMesh),
      mUseVariableRadii(false),
      mLoadBalanceMesh(false),
      mLoadBalanceFrequency(100)
{
    mpNodesOnlyMesh = static_cast<NodesOnlyMesh<DIM>* >(&(this->mrMesh));

    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh)
    : AbstractCentreBasedCellPopulation<DIM>(rMesh),
      mDeleteMesh(true),
      mUseVariableRadii(false), // will be set by serialize() method
      mLoadBalanceMesh(false),
      mLoadBalanceFrequency(100)
{
    mpNodesOnlyMesh = static_cast<NodesOnlyMesh<DIM>* >(&(this->mrMesh));
}

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::~NodeBasedCellPopulation()
{
    Clear();
    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
NodesOnlyMesh<DIM>& NodeBasedCellPopulation<DIM>::rGetMesh()
{
    return *mpNodesOnlyMesh;
}

template<unsigned DIM>
const NodesOnlyMesh<DIM>& NodeBasedCellPopulation<DIM>::rGetMesh() const
{
    return *mpNodesOnlyMesh;
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* NodeBasedCellPopulation<DIM>::GetTetrahedralMeshForPdeModifier()
{
    // Note that this code does not yet work in parallel.
    assert(PetscTools::IsSequential());

    std::vector<Node<DIM>*> temp_nodes;

    // Get the nodes of mpNodesOnlyMesh
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = mpNodesOnlyMesh->GetNodeIteratorBegin();
         node_iter != mpNodesOnlyMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        temp_nodes.push_back(new Node<DIM>(node_iter->GetIndex(), node_iter->rGetLocation(), node_iter->IsBoundaryNode()));
    }

    return new MutableMesh<DIM,DIM>(temp_nodes);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Clear()
{
    mNodePairs.clear();
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Validate()
{
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        try
        {
            this->GetCellUsingLocationIndex(node_iter->GetIndex());
        }
        catch (Exception&)
        {
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Node " << node_iter->GetIndex() << " does not appear to have a cell associated with it");
        }
    }
}

template<unsigned DIM>
Node<DIM>* NodeBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return mpNodesOnlyMesh->GetNodeOrHaloNode(index);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mpNodesOnlyMesh->SetNode(nodeIndex, rNewLocation, false);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    UpdateCellProcessLocation();

    mpNodesOnlyMesh->UpdateBoxCollection();

    if (mLoadBalanceMesh)
    {
        if ((SimulationTime::Instance()->GetTimeStepsElapsed() % mLoadBalanceFrequency) == 0)
        {
            mpNodesOnlyMesh->LoadBalanceMesh();

            UpdateCellProcessLocation();

            mpNodesOnlyMesh->UpdateBoxCollection();
        }
    }

    RefreshHaloCells();

    mpNodesOnlyMesh->CalculateInteriorNodePairs(mNodePairs);

    AddReceivedHaloCells();

    mpNodesOnlyMesh->CalculateBoundaryNodePairs(mNodePairs);

    /*
     * Update cell radii based on CellData
     */
    if (mUseVariableRadii)
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
             cell_iter != this->End();
             ++cell_iter)
        {
            double cell_radius = cell_iter->GetCellData()->GetItem("Radius");
            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
            this->GetNode(node_index)->SetRadius(cell_radius);
        }
    }

    // Make sure that everyone exits update together so that all asynchronous communications are complete.
    PetscTools::Barrier("Update");
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::UpdateMapsAfterRemesh(NodeMap& map)
{
    if (!map.IsIdentityMap())
    {
        UpdateParticlesAfterReMesh(map);

        // Update the mappings between cells and location indices
        ///\todo we want to make mCellLocationMap private - we need to find a better way of doing this
        std::map<Cell*, unsigned> old_map = this->mCellLocationMap;

        // Remove any dead pointers from the maps (needed to avoid archiving errors)
        this->mLocationCellMap.clear();
        this->mCellLocationMap.clear();

        for (std::list<CellPtr>::iterator it = this->mCells.begin();
             it != this->mCells.end();
             ++it)
        {
            unsigned old_node_index = old_map[(*it).get()];

            // This shouldn't ever happen, as the cell vector only contains living cells
            assert(!map.IsDeleted(old_node_index));

            unsigned new_node_index = map.GetNewIndex(old_node_index);
            this->SetCellUsingLocationIndex(new_node_index,*it);
        }

        this->Validate();
    }
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         )
    {
        if ((*cell_iter)->IsDead())
        {
            // Remove the node from the mesh
            num_removed++;
            unsigned location_index = mpNodesOnlyMesh->SolveNodeMapping(this->GetLocationIndexUsingCell((*cell_iter)));
            mpNodesOnlyMesh->DeleteNodePriorToReMesh(location_index);

            // Update mappings between cells and location indices
            unsigned location_index_of_removed_node = this->GetLocationIndexUsingCell((*cell_iter));
            this->RemoveCellUsingLocationIndex(location_index_of_removed_node, (*cell_iter));

            // Update vector of cells
            cell_iter = this->mCells.erase(cell_iter);
        }
        else
        {
            ++cell_iter;
        }
    }

    return num_removed;
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mpNodesOnlyMesh->AddNode(pNewNode);
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::GetNumNodes()
{
    return mpNodesOnlyMesh->GetNumNodes();
}

template<unsigned DIM>
CellPtr NodeBasedCellPopulation<DIM>::GetCellUsingLocationIndex(unsigned index)
{
    std::map<unsigned, CellPtr>::iterator iter = mLocationHaloCellMap.find(index);
    if (iter != mLocationHaloCellMap.end())
    {
        return iter->second;
    }
    else
    {
        return AbstractCellPopulation<DIM, DIM>::GetCellUsingLocationIndex(index);
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::UpdateParticlesAfterReMesh(NodeMap& rMap)
{
}

template<unsigned DIM>
std::vector< std::pair<Node<DIM>*, Node<DIM>* > >& NodeBasedCellPopulation<DIM>::rGetNodePairs()
{
    return mNodePairs;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MechanicsCutOffLength>" << mpNodesOnlyMesh->GetMaximumInteractionDistance() << "</MechanicsCutOffLength>\n";
    *rParamsFile << "\t\t<UseVariableRadii>" << mUseVariableRadii << "</UseVariableRadii>\n";

    // Call method on direct parent class
    AbstractCentreBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter)
{
    pPopulationCountWriter->Visit(this);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}

template<unsigned DIM>
double NodeBasedCellPopulation<DIM>::GetMechanicsCutOffLength()
{
    return mpNodesOnlyMesh->GetMaximumInteractionDistance();
}

template<unsigned DIM>
bool NodeBasedCellPopulation<DIM>::GetUseVariableRadii()
{
    return mUseVariableRadii;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetUseVariableRadii(bool useVariableRadii)
{
    mUseVariableRadii = useVariableRadii;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetLoadBalanceMesh(bool loadBalanceMesh)
{
    mLoadBalanceMesh = loadBalanceMesh;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetLoadBalanceFrequency(unsigned loadBalanceFrequency)
{
    mLoadBalanceFrequency = loadBalanceFrequency;
}

template<unsigned DIM>
double NodeBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    return mpNodesOnlyMesh->GetWidth(rDimension);
}

template<unsigned DIM>
c_vector<double, DIM> NodeBasedCellPopulation<DIM>::GetSizeOfCellPopulation()
{
    c_vector<double, DIM> local_size = AbstractCellPopulation<DIM, DIM>::GetSizeOfCellPopulation();
    c_vector<double, DIM> global_size;

    for (unsigned i=0; i<DIM; i++)
    {
        MPI_Allreduce(&local_size[i], &global_size[i], 1, MPI_DOUBLE, MPI_MAX, PetscTools::GetWorld());
    }

    return global_size;
}

template<unsigned DIM>
std::set<unsigned> NodeBasedCellPopulation<DIM>::GetNodesWithinNeighbourhoodRadius(unsigned index, double neighbourhoodRadius)
{
    // Check neighbourhoodRadius is less than the interaction radius. If not you wont return all the correct nodes
    if (neighbourhoodRadius > mpNodesOnlyMesh->GetMaximumInteractionDistance())
    {
        EXCEPTION("neighbourhoodRadius should be less than or equal to the  the maximum interaction radius defined on the NodesOnlyMesh");
    }

    std::set<unsigned> neighbouring_node_indices;

    // Get location
    Node<DIM>* p_node_i = this->GetNode(index);
    const c_vector<double, DIM>& r_node_i_location = p_node_i->rGetLocation();

    // Check the mNodeNeighbours has been set up correctly
    if (!(p_node_i->GetNeighboursSetUp()))
    {
        EXCEPTION("mNodeNeighbours not set up. Call Update() before GetNodesWithinNeighbourhoodRadius()");
    }

    // Get set of 'candidate' neighbours.
    std::vector<unsigned>& near_nodes = p_node_i->rGetNeighbours();

    // Find which ones are actually close
    for (std::vector<unsigned>::iterator iter = near_nodes.begin();
         iter != near_nodes.end();
         ++iter)
    {
        // Be sure not to return the index itself.
        if ((*iter) != index)
        {
            Node<DIM>* p_node_j = this->GetNode((*iter));

            // Get the location of this node
            const c_vector<double, DIM>& r_node_j_location = p_node_j->rGetLocation();

            // Get the vector the two nodes (using GetVectorFromAtoB to catch periodicities etc.)
            c_vector<double, DIM> node_to_node_vector = mpNodesOnlyMesh->GetVectorFromAtoB(r_node_j_location, r_node_i_location);

            // Calculate the distance between the two nodes
            double distance_between_nodes = norm_2(node_to_node_vector);

            // If the cell j is within the neighbourhood of radius neighbourhoodRadius
            // of cell i
            if (distance_between_nodes <= neighbourhoodRadius)// + DBL_EPSILSON)
            {
                // ...then add this node index to the set of neighbouring node indices
                neighbouring_node_indices.insert((*iter));
            }
        }
    }

    return neighbouring_node_indices;
}

template<unsigned DIM>
std::set<unsigned> NodeBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    std::set<unsigned> neighbouring_node_indices;

    // Get location and radius of node
    Node<DIM>* p_node_i = this->GetNode(index);
    const c_vector<double, DIM>& r_node_i_location = p_node_i->rGetLocation();
    double radius_of_cell_i = p_node_i->GetRadius();

    // Check the mNodeNeighbours has been set up correctly
    if (!(p_node_i->GetNeighboursSetUp()))
    {
        EXCEPTION("mNodeNeighbours not set up. Call Update() before GetNeighbouringNodeIndices()");
    }

    // Make sure that the max_interaction distance is smaller than or equal to the box collection size
    if (!(radius_of_cell_i * 2.0 <= mpNodesOnlyMesh->GetMaximumInteractionDistance()))
    {
        EXCEPTION("mpNodesOnlyMesh::mMaxInteractionDistance is smaller than twice the radius of cell " << index << " (" << radius_of_cell_i << ") so interactions may be missed. Make the cut-off larger to avoid errors.");
    }

    // Get set of 'candidate' neighbours
    std::vector<unsigned>& near_nodes = p_node_i->rGetNeighbours();

    // Find which ones are actually close
    for (std::vector<unsigned>::iterator iter = near_nodes.begin();
         iter != near_nodes.end();
         ++iter)
    {
        // Be sure not to return the index itself
        if ((*iter) != index)
        {
            Node<DIM>* p_node_j = this->GetNode((*iter));

            // Get the location of this node
            const c_vector<double, DIM>& r_node_j_location = p_node_j->rGetLocation();

            // Get the vector the two nodes (using GetVectorFromAtoB to catch periodicities etc.)
            c_vector<double, DIM> node_to_node_vector = mpNodesOnlyMesh->GetVectorFromAtoB(r_node_j_location, r_node_i_location);

            // Calculate the distance between the two nodes
            double distance_between_nodes = norm_2(node_to_node_vector);

            // Get the radius of the cell corresponding to this node
            double radius_of_cell_j = p_node_j->GetRadius();

            // If the cells are close enough to exert a force on each other...
            double max_interaction_distance = radius_of_cell_i + radius_of_cell_j;

            // Make sure that the max_interaction distance is smaller than or equal to the box collection size
            if (!(max_interaction_distance <= mpNodesOnlyMesh->GetMaximumInteractionDistance()))
            {
                EXCEPTION("mpNodesOnlyMesh::mMaxInteractionDistance is smaller than the sum of radius of cell " << index << " (" << radius_of_cell_i << ") and cell " << (*iter) << " (" << radius_of_cell_j <<"). Make the cut-off larger to avoid errors.");
            }
            if (distance_between_nodes <= max_interaction_distance)// + DBL_EPSILSON) //Assumes that max_interaction_distance is of order 1
            {
                // ...then add this node index to the set of neighbouring node indices
                neighbouring_node_indices.insert((*iter));
            }
        }
    }

    return neighbouring_node_indices;
}

template<unsigned DIM>
double NodeBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Not implemented or tested in 1D
    assert(DIM==2 ||DIM==3); // LCOV_EXCL_LINE

    // Get node index corresponding to this cell
    unsigned node_index = this->GetLocationIndexUsingCell(pCell);
    Node<DIM>* p_node = this->GetNode(node_index);

    // Get cell radius
    double cell_radius = p_node->GetRadius();

    // Begin code to approximate cell volume
    double averaged_cell_radius = 0.0;
    unsigned num_cells = 0;

    // Get the location of this node
    const c_vector<double, DIM>& r_node_i_location = GetNode(node_index)->rGetLocation();

    // Get the set of node indices corresponding to this cell's neighbours
    std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(node_index);

    // THe number of neighbours in equilibrium configuration, from sphere packing problem
    unsigned num_neighbours_equil;
    if (DIM==2)
    {
        num_neighbours_equil = 6;
    }
    else
    {
        assert(DIM==3);
        num_neighbours_equil = 12;
    }

    // Loop over this set
    for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
         iter != neighbouring_node_indices.end();
         ++iter)
    {
        Node<DIM>* p_node_j = this->GetNode(*iter);

        // Get the location of the neighbouring node
        const c_vector<double, DIM>& r_node_j_location = p_node_j->rGetLocation();

        double neighbouring_cell_radius = p_node_j->GetRadius();

        // If this throws then you may not be considering all cell interactions use a larger cut off length
        assert(cell_radius+neighbouring_cell_radius<mpNodesOnlyMesh->GetMaximumInteractionDistance());

        // Calculate the distance between the two nodes and add to cell radius
        double separation = norm_2(mpNodesOnlyMesh->GetVectorFromAtoB(r_node_j_location, r_node_i_location));

        if (separation < cell_radius+neighbouring_cell_radius)
        {
            // The effective radius is the mid point of the overlap
            averaged_cell_radius = averaged_cell_radius + cell_radius - (cell_radius+neighbouring_cell_radius-separation)/2.0;
            num_cells++;
        }
    }
    if (num_cells < num_neighbours_equil)
    {
        averaged_cell_radius += (num_neighbours_equil-num_cells)*cell_radius;

        averaged_cell_radius /= num_neighbours_equil;
    }
    else
    {
        averaged_cell_radius /= num_cells;
    }
    assert(averaged_cell_radius < mpNodesOnlyMesh->GetMaximumInteractionDistance()/2.0);

    cell_radius = averaged_cell_radius;

    // End code to approximate cell volume

    // Calculate cell volume from radius of cell
    double cell_volume = 0.0;
    if (DIM == 2)
    {
        cell_volume = M_PI*cell_radius*cell_radius;
    }
    else if (DIM == 3)
    {
        cell_volume = (4.0/3.0)*M_PI*cell_radius*cell_radius*cell_radius;
    }

    return cell_volume;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    // Store the present time as a string
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();

    // Make sure the nodes are ordered contiguously in memory
    NodeMap map(1 + this->mpNodesOnlyMesh->GetMaximumNodeIndex());
    this->mpNodesOnlyMesh->ReMesh(map);

    // Create mesh writer for VTK output
    VtkMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results_"+time.str(), false);
    mesh_writer.SetParallelFiles(*mpNodesOnlyMesh);

    auto num_nodes = GetNumNodes();

    // For each cell that this process owns, find the corresponding node index, which we only want to calculate once
    std::vector<unsigned> node_indices_in_cell_order;
    for (auto cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned global_index = this->GetLocationIndexUsingCell(*cell_iter);
        unsigned node_index = this->rGetMesh().SolveNodeMapping(global_index);

        node_indices_in_cell_order.emplace_back(node_index);
    }

    // Iterate over any cell writers that are present.  This is in a separate loop to below, because the writer loop
    // needs to the the outer loop.
    for (auto&& p_cell_writer : this->mCellWriters)
    {
        // Add any scalar data
        if (p_cell_writer->GetOutputScalarData())
        {
            std::vector<double> vtk_cell_data(num_nodes);

            unsigned loop_it = 0;
            for (auto cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter, ++loop_it)
            {
                unsigned node_idx = node_indices_in_cell_order[loop_it];
                vtk_cell_data[node_idx] = p_cell_writer->GetCellDataForVtkOutput(*cell_iter, this);
            }

            mesh_writer.AddPointData(p_cell_writer->GetVtkCellDataName(), vtk_cell_data);
        }

        // Add any vector data
        if (p_cell_writer->GetOutputVectorData())
        {
            std::vector<c_vector<double, DIM>> vtk_cell_data(num_nodes);

            unsigned loop_it = 0;
            for (auto cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter, ++loop_it)
            {
                unsigned node_idx = node_indices_in_cell_order[loop_it];
                vtk_cell_data[node_idx] = p_cell_writer->GetVectorCellDataForVtkOutput(*cell_iter, this);
            }

            mesh_writer.AddPointData(p_cell_writer->GetVtkVectorCellDataName(), vtk_cell_data);
        }
    }

    // Process rank and cell data can be collected on a cell-by-cell basis, and both occur in the following loop
    // We assume the first cell is representative of all cells
    if (this->Begin() != this->End())  // some processes may own no cells, so can't do this->Begin()->GetCellData()
    {
        auto num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
        std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();
        std::vector<std::vector<double>> cell_data(num_cell_data_items, std::vector<double>(num_nodes));

        std::vector<double> rank(num_nodes);

        unsigned loop_it = 0;
        for (auto cell_iter = this->Begin(); cell_iter != this->End(); ++cell_iter, ++loop_it)
        {
            unsigned node_idx = node_indices_in_cell_order[loop_it];

            for (unsigned cell_data_idx = 0; cell_data_idx < num_cell_data_items; ++cell_data_idx)
            {
                cell_data[cell_data_idx][node_idx] = cell_iter->GetCellData()->GetItem(cell_data_names[cell_data_idx]);
            }

            rank[node_idx] = (PetscTools::GetMyRank());
        }

        // Add point data to writers
        mesh_writer.AddPointData("Process rank", rank);
        for (unsigned cell_data_idx = 0; cell_data_idx < num_cell_data_items; ++cell_data_idx)
        {
            mesh_writer.AddPointData(cell_data_names[cell_data_idx], cell_data[cell_data_idx]);
        }
    }

    mesh_writer.WriteFilesUsingMesh(*mpNodesOnlyMesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << R"(" group="" part="0" file="results_)";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    if (PetscTools::IsSequential())
    {
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
    else
    {
        // Parallel vtu files  .vtu -> .pvtu
        *(this->mpVtkMetaFile) << ".pvtu\"/>\n";
    }
#endif //CHASTE_VTK
}

template<unsigned DIM>
CellPtr NodeBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    assert(pNewCell);

    // Add new cell to population
    CellPtr p_created_cell = AbstractCentreBasedCellPopulation<DIM>::AddCell(pNewCell, pParentCell);
    assert(p_created_cell == pNewCell);

    // Then set the new cell radius in the NodesOnlyMesh
    unsigned parent_node_index = this->GetLocationIndexUsingCell(pParentCell);
    Node<DIM>* p_parent_node = this->GetNode(parent_node_index);

    unsigned new_node_index = this->GetLocationIndexUsingCell(p_created_cell);
    Node<DIM>* p_new_node = this->GetNode(new_node_index);

    p_new_node->SetRadius(p_parent_node->GetRadius());

    // Return pointer to new cell
    return p_created_cell;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddMovedCell(CellPtr pCell, boost::shared_ptr<Node<DIM> > pNode)
{
    // Create a new node
    mpNodesOnlyMesh->AddMovedNode(pNode);

    // Update cells vector
    this->mCells.push_back(pCell);

    // Update mappings between cells and location indices
    this->AddCellUsingLocationIndex(pNode->GetIndex(), pCell);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::DeleteMovedCell(unsigned index)
{
    mpNodesOnlyMesh->DeleteMovedNode(index);

    // Update vector of cells
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        if (this->GetLocationIndexUsingCell(*cell_iter) == index)
        {
            // Update mappings between cells and location indices
            this->RemoveCellUsingLocationIndex(index, (*cell_iter));
            cell_iter = this->mCells.erase(cell_iter);

            break;
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SendCellsToNeighbourProcesses()
{
    MPI_Status status;

    if (!PetscTools::AmTopMost())
    {
        boost::shared_ptr<std::vector<std::pair<CellPtr, Node<DIM>* > > > p_cells_right(&mCellsToSendRight, null_deleter());
        mpCellsRecvRight = mRightCommunicator.SendRecvObject(p_cells_right, PetscTools::GetMyRank() + 1, mCellCommunicationTag, PetscTools::GetMyRank() + 1, mCellCommunicationTag, status);
    }
    if (!PetscTools::AmMaster())
    {
        boost::shared_ptr<std::vector<std::pair<CellPtr, Node<DIM>* > > > p_cells_left(&mCellsToSendLeft, null_deleter());
        mpCellsRecvLeft = mLeftCommunicator.SendRecvObject(p_cells_left, PetscTools::GetMyRank() - 1, mCellCommunicationTag, PetscTools::GetMyRank() - 1, mCellCommunicationTag, status);
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::NonBlockingSendCellsToNeighbourProcesses()
{
    if (!PetscTools::AmTopMost())
    {
        boost::shared_ptr<std::vector<std::pair<CellPtr, Node<DIM>* > > > p_cells_right(&mCellsToSendRight, null_deleter());
        int tag = SmallPow(2u, 1+ PetscTools::GetMyRank() ) * SmallPow (3u, 1 + PetscTools::GetMyRank() + 1);
        mRightCommunicator.ISendObject(p_cells_right, PetscTools::GetMyRank() + 1, tag);
    }
    if (!PetscTools::AmMaster())
    {
        int tag = SmallPow (2u, 1 + PetscTools::GetMyRank() ) * SmallPow (3u, 1 + PetscTools::GetMyRank() - 1);
        boost::shared_ptr<std::vector<std::pair<CellPtr, Node<DIM>* > > > p_cells_left(&mCellsToSendLeft, null_deleter());
        mLeftCommunicator.ISendObject(p_cells_left, PetscTools::GetMyRank() - 1, tag);
    }
    // Now post receives to start receiving data before returning.
    if (!PetscTools::AmTopMost())
    {
        int tag = SmallPow (3u, 1 + PetscTools::GetMyRank() ) * SmallPow (2u, 1+ PetscTools::GetMyRank() + 1);
        mRightCommunicator.IRecvObject(PetscTools::GetMyRank() + 1, tag);
    }
    if (!PetscTools::AmMaster())
    {
        int tag = SmallPow (3u, 1 + PetscTools::GetMyRank() ) * SmallPow (2u, 1+ PetscTools::GetMyRank() - 1);
        mLeftCommunicator.IRecvObject(PetscTools::GetMyRank() - 1, tag);
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::GetReceivedCells()
{
    if (!PetscTools::AmTopMost())
    {
        mpCellsRecvRight = mRightCommunicator.GetRecvObject();
    }
    if (!PetscTools::AmMaster())
    {
        mpCellsRecvLeft = mLeftCommunicator.GetRecvObject();
    }
}

template<unsigned DIM>
std::pair<CellPtr, Node<DIM>* > NodeBasedCellPopulation<DIM>::GetCellNodePair(unsigned nodeIndex)
{
    Node<DIM>* p_node = this->GetNode(nodeIndex);

    CellPtr p_cell = this->GetCellUsingLocationIndex(nodeIndex);

    std::pair<CellPtr, Node<DIM>* > new_pair(p_cell, p_node);

    return new_pair;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddNodeAndCellToSendRight(unsigned nodeIndex)
{
    std::pair<CellPtr, Node<DIM>* > pair = GetCellNodePair(nodeIndex);

    mCellsToSendRight.push_back(pair);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddNodeAndCellToSendLeft(unsigned nodeIndex)
{
    std::pair<CellPtr, Node<DIM>* > pair = GetCellNodePair(nodeIndex);

    mCellsToSendLeft.push_back(pair);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddReceivedCells()
{
    if (!PetscTools::AmMaster())
    {
        for (typename std::vector<std::pair<CellPtr, Node<DIM>* > >::iterator iter = mpCellsRecvLeft->begin();
             iter != mpCellsRecvLeft->end();
             ++iter)
        {
            // Make a shared pointer to the node to make sure it is correctly deleted.
            boost::shared_ptr<Node<DIM> > p_node(iter->second);
            AddMovedCell(iter->first, p_node);
        }
    }
    if (!PetscTools::AmTopMost())
    {
        for (typename std::vector<std::pair<CellPtr, Node<DIM>* > >::iterator iter = mpCellsRecvRight->begin();
             iter != mpCellsRecvRight->end();
             ++iter)
        {
            boost::shared_ptr<Node<DIM> > p_node(iter->second);
            AddMovedCell(iter->first, p_node);
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::UpdateCellProcessLocation()
{
    mpNodesOnlyMesh->ResizeBoxCollection();

    mpNodesOnlyMesh->CalculateNodesOutsideLocalDomain();

    std::vector<unsigned> nodes_to_send_right = mpNodesOnlyMesh->rGetNodesToSendRight();
    AddCellsToSendRight(nodes_to_send_right);

    std::vector<unsigned> nodes_to_send_left = mpNodesOnlyMesh->rGetNodesToSendLeft();
    AddCellsToSendLeft(nodes_to_send_left);

    // Post non-blocking send / receives so communication on both sides can start.
    SendCellsToNeighbourProcesses();

    // Post blocking receive calls that wait until communication complete.
    //GetReceivedCells();

    for (std::vector<unsigned>::iterator iter = nodes_to_send_right.begin();
         iter != nodes_to_send_right.end();
         ++iter)
    {
        DeleteMovedCell(*iter);
    }

    for (std::vector<unsigned>::iterator iter = nodes_to_send_left.begin();
         iter != nodes_to_send_left.end();
         ++iter)
    {
        DeleteMovedCell(*iter);
    }

    AddReceivedCells();

    NodeMap map(1 + mpNodesOnlyMesh->GetMaximumNodeIndex());
    mpNodesOnlyMesh->ReMesh(map);
    UpdateMapsAfterRemesh(map);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::RefreshHaloCells()
{
    mpNodesOnlyMesh->ClearHaloNodes();

    mHaloCells.clear();
    mHaloCellLocationMap.clear();
    mLocationHaloCellMap.clear();

    std::vector<unsigned> halos_to_send_right = mpNodesOnlyMesh->rGetHaloNodesToSendRight();
    AddCellsToSendRight(halos_to_send_right);

    std::vector<unsigned> halos_to_send_left = mpNodesOnlyMesh->rGetHaloNodesToSendLeft();
    AddCellsToSendLeft(halos_to_send_left);

    NonBlockingSendCellsToNeighbourProcesses();
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddCellsToSendRight(std::vector<unsigned>& cellLocationIndices)
{
    mCellsToSendRight.clear();

    for (unsigned i=0; i < cellLocationIndices.size(); i++)
    {
        AddNodeAndCellToSendRight(cellLocationIndices[i]);
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddCellsToSendLeft(std::vector<unsigned>& cellLocationIndices)
{
    mCellsToSendLeft.clear();

    for (unsigned i=0; i < cellLocationIndices.size(); i++)
    {
        AddNodeAndCellToSendLeft(cellLocationIndices[i]);
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddReceivedHaloCells()
{
    GetReceivedCells();

    if (!PetscTools::AmMaster())
    {
        for (typename std::vector<std::pair<CellPtr, Node<DIM>* > >::iterator iter = mpCellsRecvLeft->begin();
                iter != mpCellsRecvLeft->end();
                ++iter)
        {
            boost::shared_ptr<Node<DIM> > p_node(iter->second);
            AddHaloCell(iter->first, p_node);

        }
    }
    if (!PetscTools::AmTopMost())
    {
        for (typename std::vector<std::pair<CellPtr, Node<DIM>* > >::iterator iter = mpCellsRecvRight->begin();
                iter != mpCellsRecvRight->end();
                ++iter)
        {
            boost::shared_ptr<Node<DIM> > p_node(iter->second);
            AddHaloCell(iter->first, p_node);
        }
    }

    mpNodesOnlyMesh->AddHaloNodesToBoxes();
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AddHaloCell(CellPtr pCell, boost::shared_ptr<Node<DIM> > pNode)
{
    mHaloCells.push_back(pCell);

    mpNodesOnlyMesh->AddHaloNode(pNode);

    mHaloCellLocationMap[pCell] = pNode->GetIndex();

    mLocationHaloCellMap[pNode->GetIndex()] = pCell;
}

// Explicit instantiation
template class NodeBasedCellPopulation<1>;
template class NodeBasedCellPopulation<2>;
template class NodeBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulation)
