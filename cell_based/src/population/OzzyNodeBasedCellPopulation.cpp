/*

Copyright (c) 2005-2012, University of Oxford.
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
#include "VtkMeshWriter.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh,
                                      bool validate)
    : AbstractCentreBasedCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDeleteMesh(deleteMesh),
      mMechanicsCutOffLength(DBL_MAX)
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
      mMechanicsCutOffLength(DBL_MAX) // will be set by serialize() method
{
	mpNodesOnlyMesh = static_cast<NodesOnlyMesh<DIM>* >(&(this->mrMesh));
    // No Validate() because the cells are not associated with the cell population yet in archiving
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
    std::vector<bool> validated_node(GetNumNodes());
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        validated_node[i] = false;
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell((*cell_iter));
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
	mpNodesOnlyMesh->SetUpBoxCollection(cutOffLength, domainSize);
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

    for (unsigned i=0; i<this->mrMesh.GetNumNodes(); i++)
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
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
	mpNodesOnlyMesh->GetNode(nodeIndex)->SetPoint(rNewLocation);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    NodeMap map(this->mrMesh.GetNumAllNodes());
    mpNodesOnlyMesh->ReMesh(map);

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

    mpNodesOnlyMesh->SetMeshHasChangedSinceLoading();

    mpNodesOnlyMesh->ClearBoxCollection();

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

    mpNodesOnlyMesh->SetUpBoxCollection(mMechanicsCutOffLength, domain_size);

    mpNodesOnlyMesh->CalculateNodePairs(mNodePairs, mNodeNeighbours);
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
            mpNodesOnlyMesh->DeleteNodePriorToReMesh(this->GetLocationIndexUsingCell((*it)));

            // Update mappings between cells and location indices
            unsigned location_index_of_removed_node = this->GetLocationIndexUsingCell((*it));
            this->RemoveCellUsingLocationIndex(location_index_of_removed_node, (*it));

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
    return mpNodesOnlyMesh->AddNode(pNewNode);
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::GetNumNodes()
{
    return mpNodesOnlyMesh->GetNumAllNodes();
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::UpdateParticlesAfterReMesh(NodeMap& rMap)
{
}

template<unsigned DIM>
BoxCollection<DIM>* NodeBasedCellPopulation<DIM>::GetBoxCollection()
{
    return mpNodesOnlyMesh->GetBoxCollection();
}

template<unsigned DIM>
std::set< std::pair<Node<DIM>*, Node<DIM>* > >& NodeBasedCellPopulation<DIM>::rGetNodePairs()
{
//    if (mNodePairs.empty())
//    {
//        EXCEPTION("No node pairs set up, rGetNodePairs probably called before Update");
//    }
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
    // Check the mNodeNeighbours has been set up correctly.
    if(mNodeNeighbours.empty())
    {
        EXCEPTION("mNodeNeighbours not set up. Call Update() before GetNeighbouringNodeIndices()");
    }

    std::set<unsigned> neighbouring_node_indices;

    // Get location and radius of node
    c_vector<double, DIM> node_i_location = this->GetNode(index)->rGetLocation();
    double radius_of_cell_i = mpNodesOnlyMesh->GetCellRadius(index);

    // Make sure that the max_interaction distance is smaller than the box collection size
	if(!(radius_of_cell_i * 2.0 < mMechanicsCutOffLength))
	{
		EXCEPTION("mMechanicsCutOffLength is smaller than twice the radius of cell " << index << " (" << radius_of_cell_i << ") so interactions may be missed. Make the cut-off larger to avoid errors.");
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
            // Get the location of this node
            c_vector<double, DIM> node_j_location = this->GetNode((*iter))->rGetLocation();

            // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
            c_vector<double, DIM> unit_vector = node_i_location - node_j_location;

            // Calculate the distance between the two nodes
            double distance_between_nodes = norm_2(unit_vector);

            // Get the radius of the cell corresponding to this node
            double radius_of_cell_j = mpNodesOnlyMesh->GetCellRadius((*iter));

            // If the cells are close enough to exert a force on each other...
            double max_interaction_distance = radius_of_cell_i + radius_of_cell_j;

            // Make sure that the max_interaction distance is smaller than the box collection size
            if(!(max_interaction_distance < mMechanicsCutOffLength))
            {
                EXCEPTION("mMechanicsCutOffLength is smaller than the sum of radius of cell " << index << " (" << radius_of_cell_i << ") and cell " << (*iter) << " (" << radius_of_cell_j <<"). Make the cut-off larger to avoid errors.");
            }
            if (distance_between_nodes < max_interaction_distance)
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

    double cell_radius=0.0;
   unsigned num_cells = 0;

    // Get cell radius
//    double current_cell_radius = 0.5; // mpNodesOnlyMesh->GetCellRadius(node_index); // TODO as currently all springs assume radius = 0.5

    // Get the location of this node
	c_vector<double, DIM> node_i_location = GetNode(node_index)->rGetLocation();

	// Get the set of node indices corresponding to this cell's neighbours
	std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(node_index);

		// Loop over this set
		for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
			 iter != neighbouring_node_indices.end();
			 ++iter)
		{
			// Get the location of this node
			c_vector<double, DIM> node_j_location = GetNode(*iter)->rGetLocation();

			// Calculate the distance between the two nodes and add to cell radius
			double seperation = norm_2(node_j_location - node_i_location);

			if(seperation<1.0)//mMechanicsCutOffLength)
			{
				cell_radius = cell_radius + seperation/2.0;
				num_cells++;
			}
		}
//		PRINT_2_VARIABLES(cell_radius,neighbouring_node_indices.size());
//
//		PRINT_VARIABLE(cell_radius);


	if (num_cells == 0)
	{
		cell_radius = 0.5;
	}
	else
	{
		cell_radius /= num_cells;
	}
	assert(cell_radius<mMechanicsCutOffLength/2.0);

    // Calculate cell volume from radius
    double cell_volume = 0.0;
    if (DIM == 2)
    {
        cell_volume = M_PI*cell_radius*cell_radius;
    }
    else if (DIM == 3)
    {
      	cell_volume = (4.0/3.0)*M_PI*cell_radius*cell_radius*cell_radius;
      	if (cell_radius<0.3)
      	{
      	 PRINT_2_VARIABLES(cell_radius,cell_volume);
      	}
    }

//        ///////////////////////////////////////////////////////////////////////////////
//        //Remove volume due to intersections with other cells Search for Spherical Cap to get formulae.
//        double delta_V_c = 0.0;
//
//        // Get the location of this node
//        c_vector<double, DIM> node_i_location = GetNode(node_index)->rGetLocation();
//
//        // Get the set of node indices corresponding to this cell's neighbours
//        std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(node_index);
//
//        // Loop over this set
//        for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
//             iter != neighbouring_node_indices.end();
//             ++iter)
//        {
//            // Get the location of this node
//            c_vector<double, DIM> node_j_location = GetNode(*iter)->rGetLocation();
//
//            // Calculate the distance between the two nodes
//            double dij = norm_2(node_j_location - node_i_location);
//
//            // Get the radius of the cell corresponding to this node
//            double radius_of_cell_j = 0.5; //mpNodesOnlyMesh->GetCellRadius(*iter);// TODO as currently all springs assume radius = 0.5
//
//            // If the cells are close enough to overlap each other...
//                        if (dij < cell_radius + radius_of_cell_j)
//            {
//                // ...add the volume of the spherical cap
//                 double intersection = (dij*dij - radius_of_cell_j*radius_of_cell_j + cell_radius*cell_radius)/2.0/dij;
//                 assert(fabs(intersection-dij/2.0)<1e-6);
//                 double h = cell_radius - intersection;
//
//                 if (h<=0)
//                 {
//                	 PRINT_3_VARIABLES(node_i_location[0],node_i_location[1],node_i_location[2]);
//                	 PRINT_3_VARIABLES(node_j_location[0],node_j_location[1],node_j_location[2]);
//                	 PRINT_5_VARIABLES(cell_radius,radius_of_cell_j,dij,intersection,h);
//                	 assert (0);
//                 }
//                 delta_V_c += M_PI*h*h*(3.0*cell_radius -  h)/3.0;
//            }
//        }
//
//       cell_volume = cell_volume - delta_V_c;
//
//        if  (cell_volume<0)
//        {
//        	//PRINT_2_VARIABLES(cell_volume,delta_V_c );
//        	cell_volume =0.0;
//        }
//        ////////////////////////////////////////////////////////////////////////////////
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

    unsigned num_cell_data_items = 0;
    std::vector<std::string> cell_data_names;
    //We assume that the first cell is representative of all cells
    num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    cell_data_names = this->Begin()->GetCellData()->GetKeys();

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
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

        if (this->mOutputCellAncestors)
        {
            double ancestor_index = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_iter->GetAncestor();
            cell_ancestors[node_index] = ancestor_index;
        }
        if (this->mOutputCellProliferativeTypes)
        {
            double cell_type = cell_iter->GetCellProliferativeType();
            cell_types[node_index] = cell_type;
        }
        if (this->mOutputCellMutationStates)
        {
            double mutation_state = cell_iter->GetMutationState()->GetColour();

            CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection();
            CellPropertyCollection label_collection = collection.GetProperties<CellLabel>();

            if (label_collection.GetSize()==1 )
			{
				boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(label_collection.GetProperty());
				mutation_state = p_label->GetColour();
			}

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
            double cell_radius = mpNodesOnlyMesh->GetCellRadius(node_index);
            cell_radii[node_index] = cell_radius;
        }

        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cellwise_data[var][node_index] = cell_iter->GetCellData()->GetItem(cell_data_names[var]); 
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
    ///\todo set the correct cell radius and properly test this (#1808)
    unsigned node_index = this->GetLocationIndexUsingCell(p_created_cell);
    mpNodesOnlyMesh->SetCellRadius(node_index, 0.5);//#2100 Default the cell radius of 0.5

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
