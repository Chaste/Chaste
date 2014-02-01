/*

Copyright (c) 2005-2014, University of Oxford.
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

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellMutationStatesWriter.hpp"

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
            EXCEPTION("Node " << node_iter->GetIndex() << " does not appear to have a cell associated with it");
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
    mpNodesOnlyMesh->GetNode(nodeIndex)->SetPoint(rNewLocation);
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

    mpNodesOnlyMesh->CalculateInteriorNodePairs(mNodePairs, mNodeNeighbours);

    AddReceivedHaloCells();

    mpNodesOnlyMesh->CalculateBoundaryNodePairs(mNodePairs, mNodeNeighbours);

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
    *rParamsFile << "\t\t<UseVariableRadii>" << mUseVariableRadii <<
"</UseVariableRadii>\n";

    // Call method on direct parent class
    AbstractCentreBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
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
std::set<unsigned> NodeBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    // Check the mNodeNeighbours has been set up correctly
    if (mNodeNeighbours.empty())
    {
        EXCEPTION("mNodeNeighbours not set up. Call Update() before GetNeighbouringNodeIndices()");
    }

    std::set<unsigned> neighbouring_node_indices;

    // Get location and radius of node
    Node<DIM>* p_node_i = this->GetNode(index);
    c_vector<double, DIM> node_i_location = p_node_i->rGetLocation();
    double radius_of_cell_i = p_node_i->GetRadius();

    // Make sure that the max_interaction distance is smaller than the box collection size
    if (!(radius_of_cell_i * 2.0 < mpNodesOnlyMesh->GetMaximumInteractionDistance()))
    {
        EXCEPTION("mpNodesOnlyMesh::mMaxInteractionDistance is smaller than twice the radius of cell " << index << " (" << radius_of_cell_i << ") so interactions may be missed. Make the cut-off larger to avoid errors.");
    }


    // Get set of 'candidate' neighbours.
    std::set<unsigned> near_nodes = mNodeNeighbours.find(index)->second;

    // Find which ones are actually close
    for (std::set<unsigned>::iterator iter = near_nodes.begin();
            iter != near_nodes.end();
            ++iter)
    {
        // Be sure not to return the index itself.
        if ((*iter) != index)
        {
            Node<DIM>* p_node_j = this->GetNode((*iter));

            // Get the location of this node
            c_vector<double, DIM> node_j_location = p_node_j->rGetLocation();

            // Get the vector the two nodes (assuming no periodicities etc.)
            c_vector<double, DIM> node_to_node_vector = node_i_location - node_j_location;

            // Calculate the distance between the two nodes
            double distance_between_nodes = norm_2(node_to_node_vector);

            // Get the radius of the cell corresponding to this node
            double radius_of_cell_j = p_node_j->GetRadius();

            // If the cells are close enough to exert a force on each other...
            double max_interaction_distance = radius_of_cell_i + radius_of_cell_j;

            // Make sure that the max_interaction distance is smaller than the box collection size
            if (!(max_interaction_distance < mpNodesOnlyMesh->GetMaximumInteractionDistance()))
            {
                EXCEPTION("mpNodesOnlyMesh::mMaxInteractionDistance is smaller than the sum of radius of cell " << index << " (" << radius_of_cell_i << ") and cell " << (*iter) << " (" << radius_of_cell_j <<"). Make the cut-off larger to avoid errors.");
            }
            if (distance_between_nodes <= max_interaction_distance)// + DBL_EPSILSON) //Assumes that max_interaction_distance is of over 1
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
    // Get node index corresponding to this cell
    unsigned node_index = this->GetLocationIndexUsingCell(pCell);
    Node<DIM>* p_node = this->GetNode(node_index);

    // Get cell radius
    double cell_radius = p_node->GetRadius();

    // Begin code to approximate cell volume
    double averaged_cell_radius = 0.0;
    unsigned num_cells = 0;

    // Get the location of this node
    c_vector<double, DIM> node_i_location = GetNode(node_index)->rGetLocation();

    // Get the set of node indices corresponding to this cell's neighbours
    std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(node_index);

    // Loop over this set
    for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
         iter != neighbouring_node_indices.end();
         ++iter)
    {
        Node<DIM>* p_node_j = this->GetNode(*iter);

        // Get the location of the neighbouring node
        c_vector<double, DIM> node_j_location = p_node_j->rGetLocation();

        double neighbouring_cell_radius = p_node_j->GetRadius();

        // If this throws then you may not be considering all cell interactions use a larger cut off length
        assert(cell_radius+neighbouring_cell_radius<mpNodesOnlyMesh->GetMaximumInteractionDistance());

        // Calculate the distance between the two nodes and add to cell radius
        double separation = norm_2(node_j_location - node_i_location);

        if (separation < cell_radius+neighbouring_cell_radius)
        {
            // The effective radius is the mid point of the overlap
            averaged_cell_radius = averaged_cell_radius + cell_radius - (cell_radius+neighbouring_cell_radius-separation)/2.0;
            num_cells++;
        }
    }
    if (num_cells == 0)
    {
        averaged_cell_radius = cell_radius;
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
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();
    VtkMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results_"+time.str(), false);

    // Make sure the nodes are ordered contiguously in memory.
    NodeMap map(1 + this->mpNodesOnlyMesh->GetMaximumNodeIndex());
    this->mpNodesOnlyMesh->ReMesh(map);

    mesh_writer.SetParallelFiles(*mpNodesOnlyMesh);

    unsigned num_nodes = GetNumNodes();
    std::vector<double> cell_types(num_nodes);
    std::vector<double> cell_ancestors(num_nodes);
    std::vector<double> cell_mutation_states(num_nodes);
    std::vector<double> cell_ages(num_nodes);
    std::vector<double> cell_cycle_phases(num_nodes);
    std::vector<double> cell_radii(num_nodes);
    std::vector<std::vector<double> > cellwise_data;
    std::vector<double> rank(num_nodes);

    unsigned num_cell_data_items = 0;
    std::vector<std::string> cell_data_names;

    // We assume that the first cell is representative of all cells
    if (num_nodes > 0)
    {
        num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
        cell_data_names = this->Begin()->GetCellData()->GetKeys();
    }

    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cellwise_data_var(num_nodes);
        cellwise_data.push_back(cellwise_data_var);
    }

    // Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned global_index = this->GetLocationIndexUsingCell(*cell_iter);

        Node<DIM>* p_node = this->GetNode(global_index);

        unsigned node_index = this->rGetMesh().SolveNodeMapping(global_index);

        if (this-> template HasWriter<CellAncestorWriter>())
        {
            double ancestor_index = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_iter->GetAncestor();
            cell_ancestors[node_index] = ancestor_index;
        }
        if (this-> template HasWriter<CellProliferativeTypesWriter>())
        {
            double cell_type = cell_iter->GetCellProliferativeType()->GetColour();
            cell_types[node_index] = cell_type;
        }
        if (this-> template HasWriter<CellMutationStatesWriter>())
        {
            double mutation_state = cell_iter->GetMutationState()->GetColour();

            CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection();
            CellPropertyCollection label_collection = collection.GetProperties<CellLabel>();

            if (label_collection.GetSize() == 1)
            {
                boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(label_collection.GetProperty());
                mutation_state = p_label->GetColour();
            }

            cell_mutation_states[node_index] = mutation_state;
        }
        if (this-> template HasWriter<CellAgesWriter>())
        {
            double age = cell_iter->GetAge();
            cell_ages[node_index] = age;
        }
        if (this-> template HasWriter<CellProliferativePhasesWriter>())
        {
            double cycle_phase = cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase();
            cell_cycle_phases[node_index] = cycle_phase;
        }
        if (this-> template HasWriter<CellVolumesWriter>())
        {
            double cell_radius = p_node->GetRadius();
            cell_radii[node_index] = cell_radius;
        }

        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cellwise_data[var][node_index] = cell_iter->GetCellData()->GetItem(cell_data_names[var]);
        }

        rank[node_index] = (PetscTools::GetMyRank());
    }

    mesh_writer.AddPointData("Process rank", rank);

    if (this-> template HasWriter<CellProliferativeTypesWriter>())
    {
        mesh_writer.AddPointData("Cell types", cell_types);
    }
    if (this-> template HasWriter<CellAncestorWriter>())
    {
        mesh_writer.AddPointData("Ancestors", cell_ancestors);
    }
    if (this-> template HasWriter<CellMutationStatesWriter>())
    {
        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
    }
    if (this-> template HasWriter<CellAgesWriter>())
    {
        mesh_writer.AddPointData("Ages", cell_ages);
    }
    if (this-> template HasWriter<CellProliferativePhasesWriter>())
    {
        mesh_writer.AddPointData("Cycle phases", cell_cycle_phases);
    }
    if (this-> template HasWriter<CellVolumesWriter>())
    {
        mesh_writer.AddPointData("Cell radii", cell_radii);
    }
    if (num_cell_data_items > 0)
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            mesh_writer.AddPointData(cell_data_names[var], cellwise_data[var]);
        }
    }

    mesh_writer.WriteFilesUsingMesh(*mpNodesOnlyMesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

template<unsigned DIM>
CellPtr NodeBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    assert(pNewCell);

    // Add new cell to cell population
    CellPtr p_created_cell = AbstractCentreBasedCellPopulation<DIM>::AddCell(pNewCell, rCellDivisionVector, pParentCell);
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

/**
 * This null deleter is for the next method, SendCellsToNeighbourProcesses so we can make a
 * shared_ptr copy of mCellsToSendx without it actually being deleted when the pointer goes out of scope.
 * We need a shared pointer to send it because ObjectCommunicator only sends/recvs shared pointers to
 * avoid memory management problems.
 */
struct null_deleter
{
    /**
     * The delete operation that does nothing.
     */
    void operator()(void const *) const
    {
    }
};

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SendCellsToNeighbourProcesses()
{
#if BOOST_VERSION < 103700
    EXCEPTION("Parallel cell-based Chaste requires Boost >= 1.37");
#else // BOOST_VERSION >= 103700
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
#endif
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::NonBlockingSendCellsToNeighbourProcesses()
{
#if BOOST_VERSION < 103700
    EXCEPTION("Parallel cell-based Chaste requires Boost >= 1.37");
#else // BOOST_VERSION >= 103700

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
#endif
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::GetReceivedCells()
{
#if BOOST_VERSION < 103700
    EXCEPTION("Parallel cell-based Chaste requires Boost >= 1.37");
#else // BOOST_VERSION >= 103700

    if (!PetscTools::AmTopMost())
    {
        mpCellsRecvRight = mRightCommunicator.GetRecvObject();
    }
    if (!PetscTools::AmMaster())
    {
        mpCellsRecvLeft = mLeftCommunicator.GetRecvObject();
    }
#endif
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

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class NodeBasedCellPopulation<1>;
template class NodeBasedCellPopulation<2>;
template class NodeBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulation)
