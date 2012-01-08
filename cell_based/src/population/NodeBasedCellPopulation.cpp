/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "NodeBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "VtkMeshWriter.hpp"

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh)
    : AbstractCentreBasedCellPopulation<DIM>(rCells, locationIndices),
      mrMesh(rMesh),
      mpBoxCollection(NULL),
      mDeleteMesh(deleteMesh),
      mMechanicsCutOffLength(DBL_MAX)
{
    Validate();
}

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh)
    : AbstractCentreBasedCellPopulation<DIM>(),
      mrMesh(rMesh),
      mpBoxCollection(NULL),
      mDeleteMesh(true),
      mMechanicsCutOffLength(DBL_MAX) // will be set by serialize() method
{
    // No Validate() because the cells are not associated with the cell population yet in archiving
}

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::~NodeBasedCellPopulation()
{
    Clear();
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}

template<unsigned DIM>
NodesOnlyMesh<DIM>& NodeBasedCellPopulation<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
const NodesOnlyMesh<DIM>& NodeBasedCellPopulation<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Clear()
{
    delete mpBoxCollection;
    mpBoxCollection = NULL;
    mNodePairs.clear();
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Validate()
{
    std::vector<bool> validated_node(GetNumNodes());
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        validated_node[i] = false;
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            EXCEPTION("Node " << i << " does not appear to have a cell associated with it");
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SplitUpIntoBoxes(double cutOffLength, c_vector<double, 2*DIM> domainSize)
{
    mpBoxCollection = new BoxCollection<DIM>(cutOffLength, domainSize);
    mpBoxCollection->SetupLocalBoxesHalfOnly();

    for (unsigned i=0; i<mrMesh.GetNumNodes(); i++)
    {
        unsigned box_index = mpBoxCollection->CalculateContainingBox(this->GetNode(i));
        mpBoxCollection->rGetBox(box_index).AddNode(this->GetNode(i));
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::FindMaxAndMin()
{
    c_vector<double, DIM> min_posn;
    c_vector<double, DIM> max_posn;
    for (unsigned i=0; i<DIM; i++)
    {
        min_posn(i) = DBL_MAX;
        max_posn(i) = -DBL_MAX;
    }

    for (unsigned i=0; i<mrMesh.GetNumNodes(); i++)
    {
        c_vector<double, DIM> posn = this->GetNode(i)->rGetLocation();

        for (unsigned j=0; j<DIM; j++)
        {
            if (posn(j) > max_posn(j))
            {
                max_posn(j) = posn(j);
            }
            if (posn(j) < min_posn(j))
            {
                min_posn(j) = posn(j);
            }
        }
    }

    for (unsigned i=0; i<DIM; i++)
    {
        assert(min_posn(i) != DBL_MAX);
        mMinSpatialPositions(i) = min_posn(i);

        assert(max_posn(i) != -DBL_MAX);
        mMaxSpatialPositions(i) = max_posn(i);
    }
}

template<unsigned DIM>
Node<DIM>* NodeBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mrMesh.GetNode(nodeIndex)->SetPoint(rNewLocation);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    NodeMap map(mrMesh.GetNumAllNodes());
    mrMesh.ReMesh(map);

    if (!map.IsIdentityMap())
    {
        // Update the mappings between cells and location indices
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
            this->mLocationCellMap[new_node_index] = *it;
            this->mCellLocationMap[(*it).get()] = new_node_index;
        }

        this->Validate();
    }

    mrMesh.SetMeshHasChangedSinceLoading();

    if (mpBoxCollection != NULL)
    {
        delete mpBoxCollection;
    }

    FindMaxAndMin();

    // Something here to set up the domain size (max and min of each node position dimension)
    c_vector<double, 2*DIM> domain_size;
    for (unsigned i=0; i<DIM; i++)
    {
        domain_size(2*i) = mMinSpatialPositions(i);
        domain_size(2*i+1) = mMaxSpatialPositions(i);
    }

    if (mMechanicsCutOffLength == DBL_MAX)
    {
        std::string error =  std::string("NodeBasedCellPopulation cannot create boxes if the cut-off length has not been set - ")
                           + std::string("Call SetMechanicsCutOffLength on the CellPopulation ensuring it is larger than GetCutOffLength() on the force law");
        EXCEPTION(error);
    }

    /*
     * Add this parameter and suggest that mechanics systems set it.
     * Allocates memory for mpBoxCollection and does the splitting
     * and putting nodes into boxes.
     */
    SplitUpIntoBoxes(mMechanicsCutOffLength, domain_size);

    std::vector<Node<DIM>*> nodes;
    for (unsigned index=0; index<mrMesh.GetNumNodes(); index++)
    {
        Node<DIM>* p_node = mrMesh.GetNode(index);
        nodes.push_back(p_node);
    }
    mpBoxCollection->CalculateNodePairs(nodes, mNodePairs);
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<CellPtr>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if ((*it)->IsDead())
        {
            // Remove the node from the mesh
            num_removed++;
            mrMesh.DeleteNodePriorToReMesh(this->mCellLocationMap[(*it).get()]);

            // Update mappings between cells and location indices
            unsigned location_index_of_removed_node = this->mCellLocationMap[(*it).get()];
            this->mCellLocationMap.erase((*it).get());
            this->mLocationCellMap.erase(location_index_of_removed_node);

            // Update vector of cells
            it = this->mCells.erase(it);
            --it;
        }
    }

    return num_removed;
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mrMesh.AddNode(pNewNode);
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::GetNumNodes()
{
    return mrMesh.GetNumAllNodes();
}

template<unsigned DIM>
BoxCollection<DIM>* NodeBasedCellPopulation<DIM>::GetBoxCollection()
{
    return mpBoxCollection;
}

template<unsigned DIM>
std::set< std::pair<Node<DIM>*, Node<DIM>* > >& NodeBasedCellPopulation<DIM>::rGetNodePairs()
{
    if (mNodePairs.empty())
    {
        EXCEPTION("No node pairs set up, rGetNodePairs probably called before Update");
    }
    return mNodePairs;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MechanicsCutOffLength>" << mMechanicsCutOffLength << "</MechanicsCutOffLength>\n";

    // Call method on direct parent class
    AbstractCentreBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetMechanicsCutOffLength(double mechanicsCutOffLength)
{
    assert(mechanicsCutOffLength > 0.0);
    mMechanicsCutOffLength = mechanicsCutOffLength;
}

template<unsigned DIM>
double NodeBasedCellPopulation<DIM>::GetMechanicsCutOffLength()
{
    return mMechanicsCutOffLength;
}

template<unsigned DIM>
double NodeBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Update the member variables mMinSpatialPositions and mMaxSpatialPositions
    FindMaxAndMin();

    // Compute the maximum distance between any nodes in this dimension
    double width = mMaxSpatialPositions(rDimension) - mMinSpatialPositions(rDimension);

    return width;
}

template<unsigned DIM>
std::set<unsigned> NodeBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    // Get the location of this node
    c_vector<double, DIM> node_i_location = this->GetNode(index)->rGetLocation();

    // Get the radius of the cell corresponding to this node
    double radius_of_cell_i = mrMesh.GetCellRadius(index);

    // Loop over cells in the population
    std::set<unsigned> neighbouring_node_indices;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned node_j_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Only return the neighbours, not the original node
        if (node_j_index != index)
        {
            // Get the location of this node
            c_vector<double, DIM> node_j_location = this->GetNode(node_j_index)->rGetLocation();

            // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
            c_vector<double, DIM> unit_vector = node_i_location - node_j_location;

            // Calculate the distance between the two nodes
            double distance_between_nodes = norm_2(unit_vector);

            // Get the radius of the cell corresponding to this node
            double radius_of_cell_j = mrMesh.GetCellRadius(node_j_index);

            // If the cells are close enough to exert a force on each other...
            double max_interaction_distance = radius_of_cell_i + radius_of_cell_j;
            if (distance_between_nodes < max_interaction_distance)
            {
                // ...then add this node index to the set of neighbouring node indices
                neighbouring_node_indices.insert(node_j_index);
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

    // Get cell radius
    double cell_radius = mrMesh.GetCellRadius(node_index);

    // Get cell volume from radius
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
void NodeBasedCellPopulation<DIM>::WriteCellVolumeResultsToFile()
{
    assert(DIM==2 || DIM==3);

    // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    // Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the index of the corresponding node in mrMesh
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Write node index to file
        *(this->mpCellVolumesFile) << node_index << " ";

        // Write cell ID to file
        *(this->mpCellVolumesFile) << cell_iter->GetCellId() << " ";

        // Write node location to file
        c_vector<double, DIM> node_location = this->GetNode(node_index)->rGetLocation();
        for (unsigned i=0; i<DIM; i++)
        {
            *(this->mpCellVolumesFile) << node_location[i] << " ";
        }

        // Write cell volume (in 3D) or area (in 2D) to file
        double cell_volume = this->GetVolumeOfCell(*cell_iter);
        *(this->mpCellVolumesFile) << cell_volume << " ";
    }

    *(this->mpCellVolumesFile) << "\n";
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();
    VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);

    unsigned num_nodes = GetNumNodes();
    std::vector<double> cell_types(num_nodes);
    std::vector<double> cell_ancestors(num_nodes);
    std::vector<double> cell_mutation_states(num_nodes);
    std::vector<double> cell_ages(num_nodes);
    std::vector<double> cell_cycle_phases(num_nodes);
    std::vector<double> cell_radii(num_nodes);
    std::vector<std::vector<double> > cellwise_data;

    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
        unsigned num_variables = p_data->GetNumVariables();
        for (unsigned var=0; var<num_variables; var++)
        {
            std::vector<double> cellwise_data_var(num_nodes);
            cellwise_data.push_back(cellwise_data_var);
        }
    }

    // Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

        if (this->mOutputCellAncestors)
        {
            double ancestor_index = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_iter->GetAncestor();
            cell_ancestors[node_index] = ancestor_index;
        }
        if (this->mOutputCellProliferativeTypes)
        {
            double cell_type = cell_iter->GetCellCycleModel()->GetCellProliferativeType();
            cell_types[node_index] = cell_type;
        }
        if (this->mOutputCellMutationStates)
        {
            double mutation_state = cell_iter->GetMutationState()->GetColour();
            cell_mutation_states[node_index] = mutation_state;
        }
        if (this->mOutputCellAges)
        {
            double age = cell_iter->GetAge();
            cell_ages[node_index] = age;
        }
        if (this->mOutputCellCyclePhases)
        {
            double cycle_phase = cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase();
            cell_cycle_phases[node_index] = cycle_phase;
        }
        if (this->mOutputCellVolumes)
        {
            double cell_radius = mrMesh.GetCellRadius(node_index);
            cell_radii[node_index] = cell_radius;
        }
        if (CellwiseData<DIM>::Instance()->IsSetUp())
        {
            CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
            unsigned num_variables = p_data->GetNumVariables();
            for (unsigned var=0; var<num_variables; var++)
            {
                cellwise_data[var][node_index] = p_data->GetValue(*cell_iter, var);
            }
        }
    }

    if (this->mOutputCellProliferativeTypes)
    {
        mesh_writer.AddPointData("Cell types", cell_types);
    }
    if (this->mOutputCellAncestors)
    {
        mesh_writer.AddPointData("Ancestors", cell_ancestors);
    }
    if (this->mOutputCellMutationStates)
    {
        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
    }
    if (this->mOutputCellAges)
    {
        mesh_writer.AddPointData("Ages", cell_ages);
    }
    if (this->mOutputCellCyclePhases)
    {
        mesh_writer.AddPointData("Cycle phases", cell_cycle_phases);
    }
    if (this->mOutputCellVolumes)
    {
        mesh_writer.AddPointData("Cell radii", cell_radii);
    }
    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            std::stringstream data_name;
            data_name << "Cellwise data " << var;
            std::vector<double> cellwise_data_var = cellwise_data[var];
            mesh_writer.AddPointData(data_name.str(), cellwise_data_var);
        }
    }

    mesh_writer.WriteFilesUsingMesh(mrMesh);

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
    ///\todo set the correct cell radius and properly test this (#1808)
    unsigned node_index = this->GetLocationIndexUsingCell(p_created_cell);
    mrMesh.SetCellRadius(node_index, 1.0);

    // Return pointer to new cell
    return p_created_cell;
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
