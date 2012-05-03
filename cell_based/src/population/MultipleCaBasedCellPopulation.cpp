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

#include "MultipleCaBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"
#include "NodesOnlyMesh.hpp"
#include "Exception.hpp"

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::Validate()
{
    // Check each node has at most one cell associated with it
    std::vector<unsigned> validated_nodes = std::vector<unsigned>(this->GetNumNodes(), 0);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned node_index = GetLocationIndexUsingCell(*cell_iter);
        validated_nodes[node_index]++;
        //Check nodes with cells are marked as non-empty
        assert(mEmptySites[node_index] == false);
    }

    for (unsigned i=0; i<validated_nodes.size(); i++)
    {
        if (validated_nodes[i] > 1)
        {
            NEVER_REACHED;
            ///\todo #2066 - This exception is not covered
            //EXCEPTION("Node " << i << " appears to have " << validated_nodes[i] << " cells associated with it");
        }
        else if (validated_nodes[i] < 1)
        {
            //Check nodes without cells are marked as empty
            assert(mEmptySites[i] == true);
        }
    }
}

template<unsigned DIM>
MultipleCaBasedCellPopulation<DIM>::MultipleCaBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                                        std::vector<CellPtr>& rCells,
                                                        const std::vector<unsigned> locationIndices,
                                                        bool deleteMesh,
                                                        bool validate)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh)
{
    // This must always be true
    assert(this->mCells.size() <= this->mrMesh.GetNumNodes());

    if (!locationIndices.empty())
    {
        // Create a set of node indices corresponding to empty sites
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices;
        std::set<unsigned> empty_site_indices;

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
                            std::inserter(empty_site_indices, empty_site_indices.begin()));

        // This method finishes and then calls Validate()
        SetEmptySites(empty_site_indices);
    }
    else
    {
        NEVER_REACHED;
        ///\todo #2066 - This code is not covered
        //mEmptySites = std::vector<bool>(this->GetNumNodes(), false);
        //Validate();
    }
}

template<unsigned DIM>
MultipleCaBasedCellPopulation<DIM>::MultipleCaBasedCellPopulation(PottsMesh<DIM>& rMesh)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
MultipleCaBasedCellPopulation<DIM>::~MultipleCaBasedCellPopulation()
{
    if (this->mDeleteMesh)
    {
        NEVER_REACHED;
        ///\todo #2066 - This code is not covered
        delete &this->mrMesh;
    }
}


template<unsigned DIM>
std::vector<bool>& MultipleCaBasedCellPopulation<DIM>::rGetEmptySites()
{
    return mEmptySites;
}

template<unsigned DIM>
bool MultipleCaBasedCellPopulation<DIM>::IsEmptySite(unsigned index)
{
    return mEmptySites[index];
}

template<unsigned DIM>
std::set<unsigned> MultipleCaBasedCellPopulation<DIM>::GetEmptySiteIndices()
{
    NEVER_REACHED;
    ///\todo #2066 - This code is not covered   
//    std::set<unsigned> empty_site_indices;
//    for (unsigned i=0; i<mEmptySites.size(); i++)
//    {
//        if (mEmptySites[i])
//        {
//            empty_site_indices.insert(i);
//        }
//    }
//    return empty_site_indices;
    std::set<unsigned> empty_site_indices;
    return (empty_site_indices);
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::SetEmptySites(const std::set<unsigned>& rEmptySiteIndices)
{
    // Reinitialise all entries of mEmptySites to false
    mEmptySites = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mEmptySites
    for (std::set<unsigned>::iterator iter = rEmptySiteIndices.begin();
         iter != rEmptySiteIndices.end();
         ++iter)
    {
        mEmptySites[*iter] = true;
    }

    Validate();
}

template<unsigned DIM>
PottsMesh<DIM>& MultipleCaBasedCellPopulation<DIM>::rGetMesh()
{
    return static_cast<PottsMesh<DIM>& >((this->mrMesh));
}

template<unsigned DIM>
const PottsMesh<DIM>& MultipleCaBasedCellPopulation<DIM>::rGetMesh() const
{
    return static_cast<PottsMesh<DIM>& >((this->mrMesh));
}

template<unsigned DIM>
Node<DIM>* MultipleCaBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned MultipleCaBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
c_vector<double, DIM> MultipleCaBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return this->mrMesh.GetNode(this->GetLocationIndexUsingCell(pCell))->rGetLocation();
}

template<unsigned DIM>
Node<DIM>* MultipleCaBasedCellPopulation<DIM>::GetNodeCorrespondingToCell(CellPtr pCell)
{
    return this->mrMesh.GetNode(this->GetLocationIndexUsingCell(pCell));
}

//template<unsigned DIM>
//double MultipleCaBasedCellPopulation<DIM>::GetProbabilityOfDivisionIntoTargetSite(unsigned parentNodeIndex, unsigned targetIndex, unsigned numNeighbours)
//{
//    return 1.0/((double) numNeighbours);
//}

template<unsigned DIM>
CellPtr MultipleCaBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    // Get node index corresponding to the parent cell
    unsigned parent_node_index = this->GetLocationIndexUsingCell(pParentCell);

    // Get the set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(parent_node_index);

    assert(!neighbouring_node_indices.empty());

    unsigned num_neighbours = neighbouring_node_indices.size();

//    double prob_dividing_into_this_site = GetProbabilityOfDivisionIntoTargetSite(parent_node_index, ??, num_neighbours);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    unsigned chosen_start = p_gen->randMod(num_neighbours);
    std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();

    for (unsigned i=0; i<chosen_start; i++)
    {
        neighbour_iter++;
    }

    unsigned daughter_node_index = UNSIGNED_UNSET;

    unsigned count = 0;
    while (count < num_neighbours)
    {
        bool is_empty_site = IsEmptySite(*neighbour_iter);

        if (is_empty_site)
        {
            daughter_node_index = *neighbour_iter;
            break;
        }
        else
        {
        	neighbour_iter++;

        	if (neighbour_iter == neighbouring_node_indices.end())
        	{
        		neighbour_iter = neighbouring_node_indices.begin();
        	}

        	count++;
        }
    }

    if (daughter_node_index == UNSIGNED_UNSET)
    {
        EXCEPTION("No free space to divide.");
    }

    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(daughter_node_index,p_created_cell);

    //Set node to be occupied
    mEmptySites[daughter_node_index] = false;
    return p_created_cell;
}

template<unsigned DIM>
unsigned MultipleCaBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        if ((*cell_iter)->IsDead())
        {
            NEVER_REACHED;
            ///\todo #2066 - This code is not covered   
//            // Get the index of the node corresponding to this cell
//            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
//
//            // Set this node to be an empty site
//            mEmptySites[node_index] = true;
//            this->RemoveCellUsingLocationIndex(node_index, (*cell_iter))
//
//            // Erase cell and update counter
//            num_removed++;
//            cell_iter = this->mCells.erase(cell_iter);
//            --cell_iter;
        }
    }
    return num_removed;
}
//template<unsigned DIM>
//double MultipleCaBasedCellPopulation<DIM>::CalculateProbabilityOfMoving(unsigned node_index, unsigned node_neighbour_index, double dt)
//{
//    c_vector<double, DIM> node_index_location = this->GetNode(node_index)->rGetLocation();
//    c_vector<double, DIM> node_neighbour_location = this->GetNode(node_neighbour_index)->rGetLocation();
//
//    /*
//     * \todo 0.5 in the formula must be replaced by a parameter D to be read from a file
//     */
//    return(0.01*dt/(2* pow(norm_2(this->rGetMesh().GetVectorFromAtoB(node_index_location, node_neighbour_location)), 2)));
//}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
{
	/*
	 * Iterate over cells \todo make this sweep random
	 */
	for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
	         cell_iter != this->mCells.end();
	         ++cell_iter)
	{

		/*
		 * Loop over neighbours and calculate probability of moving (make sure all probabilities are <1)
		 */
		unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
		assert(!IsEmptySite(node_index));

		// Find a random available neighbouring node to overwrite current site
		std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(node_index);
		std::vector<double> neighbouring_node_propensities;
		std::vector<double> neighbouring_node_indices_vector;

		if (!neighbouring_node_indices.empty())
		{
	        //neighbouring_node_indices_list.sort();
		    unsigned num_neighbours = neighbouring_node_indices.size();

			double probability_of_not_moving = 1.0;
			double probability_of_moving = 0.0;

            for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
                 iter != neighbouring_node_indices.end();
                 ++iter)
            {
                neighbouring_node_indices_vector.push_back(*iter);

                if (IsEmptySite(*iter))
            	{
                	// Iterating over the update rule
                    for (typename std::vector<boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > >::iterator iterRule = mUpdateRuleCollection.begin();
                         iterRule != mUpdateRuleCollection.end();
                         ++iterRule)
                    {
                    	probability_of_moving = (*iterRule)->EvaluateProbability(node_index, *iter, *this, dt, 1);
                        assert(probability_of_moving <= 1.0); // Todo Change to exception warn about timestep and parameters
                        assert(probability_of_moving >= 0.0); // Todo Change to exception warn about timestep and parameters
                    }

                    probability_of_not_moving -= probability_of_moving;
                    neighbouring_node_propensities.push_back(probability_of_moving);
            	}
            	else
            	{
            		neighbouring_node_propensities.push_back(0.0);
            	}
			}
            assert(probability_of_not_moving <= 1.0); // Todo Change to exception warn about timestep and parameters
            assert(probability_of_not_moving >= 0.0); // Todo Change to exception warn about timestep and parameters

        	/*
        	 * \todo Rescale so sum of probabilities = 1 (make sure sum of probabilities <1)
        	 */

        	/*
        	 * Sample random number to specify which move to make
        	 */
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            double random_number = p_gen->ranf();

            double total_probability = neighbouring_node_propensities[0];
            unsigned counter = 1u;
            
            ///\todo #2066 This would be less prone to error if it were a for loop
            while ((total_probability < random_number) && (counter < num_neighbours))
            {
                ///\todo #2066 At this point we need to know that counter < num_neighbours
                total_probability += neighbouring_node_propensities[counter];
                counter++; 
            }

            if (counter < num_neighbours)///\todo #2066 Check this!
            {
				unsigned chosen_neighbour_location_index = neighbouring_node_indices_vector[counter-1];///\todo #2066 Check this!

				/*
				 * Move the cell to new location
				 */
				this->MoveCellInLocationMap((*cell_iter), node_index, chosen_neighbour_location_index);

			    mEmptySites[chosen_neighbour_location_index] = false;
			    mEmptySites[node_index] = true;
            }
            //else stay in the same location
		}
		else
		{
			// Each node in the mesh must have at least one neighbour
			NEVER_REACHED;
		}


	}
}

template<unsigned DIM>
bool MultipleCaBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return false;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractCellPopulation<DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizLocationsFile = output_file_handler.OpenOutputFile("results.vizlocations");
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::CloseOutputFiles()
{
    AbstractCellPopulation<DIM>::CloseOutputFiles();
    mpVizLocationsFile->close();
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::WriteResultsToFiles()
{

    AbstractCellPopulation<DIM>::WriteResultsToFiles();

    SimulationTime* p_time = SimulationTime::Instance();

    // Write location data to file
    *mpVizLocationsFile << p_time->GetTime() << "\t";

    // Loop over cells and find associated nodes so in the same order as the cells in output files
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Hack that covers the case where the element is associated with a cell that has just been killed (#1129)

        bool node_corresponds_to_dead_cell = false;

        if (this->IsCellAttachedToLocationIndex(node_index))
        {
            node_corresponds_to_dead_cell = this->GetCellUsingLocationIndex(node_index)->IsDead();
        }

        // Write node data to file
        if (!node_corresponds_to_dead_cell)
        {
            // Write the index of the of Node the cell is associated with.
            *mpVizLocationsFile << node_index << " ";
        }
    }
    *mpVizLocationsFile << "\n";

}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::WriteCellVolumeResultsToFile()
{
    // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    // Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the index of the corresponding node in mrMesh and write to file
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
        *(this->mpCellVolumesFile) << node_index << " ";

        // Get cell ID and write to file
        unsigned cell_index = cell_iter->GetCellId();
        *(this->mpCellVolumesFile) << cell_index << " ";

        // Get node location and write to file
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
double MultipleCaBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    double cell_volume = 1.0;

    return cell_volume;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::GenerateCellResultsAndWriteToFiles()
{
    // Set up cell type counter
    unsigned num_cell_types = this->mCellProliferativeTypeCount.size();
    std::vector<unsigned> cell_type_counter(num_cell_types);
    for (unsigned i=0; i<num_cell_types; i++)
    {
        cell_type_counter[i] = 0;
    }

    // Set up cell cycle phase counter
    unsigned num_cell_cycle_phases = this->mCellCyclePhaseCount.size();
    std::vector<unsigned> cell_cycle_phase_counter(num_cell_cycle_phases);
    for (unsigned i=0; i<num_cell_cycle_phases; i++)
    {
        cell_cycle_phase_counter[i] = 0;
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        this->GenerateCellResults(this->GetLocationIndexUsingCell(*cell_iter), cell_type_counter, cell_cycle_phase_counter);
    }

    this->WriteCellResultsToFiles(cell_type_counter, cell_cycle_phase_counter);
}

template<unsigned DIM>
double MultipleCaBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    NEVER_REACHED;
    ///\todo #2066 - This code is not covered
    //double width = this->mrMesh.GetWidth(rDimension);

    //return width;
    return 0.0;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > pUpdateRule)
{
    mUpdateRuleCollection.push_back(pUpdateRule);
}

template<unsigned DIM>
const std::vector<boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > >& MultipleCaBasedCellPopulation<DIM>::rGetUpdateRuleCollection() const
{
    return mUpdateRuleCollection;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractOnLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
std::set<unsigned> MultipleCaBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    NEVER_REACHED;
    ///\todo #2066 - This exception is not covered
    //EXCEPTION("Cannot call GetNeighbouringNodeIndices() on a MultipleCaBasedCellPopulation, need to go through the PottsMesh instead");
    std::set<unsigned> neighbouring_node_indices;
    return neighbouring_node_indices;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    ///\todo #2032 Compare with MeshBasedCellPopulation::WriteVtkResultsToFile etc.
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();
    VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);

    unsigned num_nodes = GetNumNodes();
    std::vector<double> cell_types(num_nodes, -1.0);
    std::vector<double> cell_mutation_states(num_nodes, -1.0);
    std::vector<double> cell_labels(num_nodes, -1.0);
    std::vector<double> cell_ids(num_nodes, -1.0);
    std::vector<double> cell_ancestors(num_nodes, -1.0);
    std::vector<double> cell_ages(num_nodes, -1.0);
    std::vector<double> cell_cycle_phases(num_nodes, -1.0);

    for (unsigned node_index = 0; node_index < static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetNumNodes(); node_index++)
    {
        if (mEmptySites[node_index] == false)
        {
            CellPtr cell_ptr = this->GetCellUsingLocationIndex(node_index);
            cell_ids[node_index] = cell_ptr->GetCellId();

            if (this->mOutputCellAncestors)
            {
                double ancestor_index = (cell_ptr->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_ptr->GetAncestor();
                cell_ancestors[node_index] = ancestor_index;
            }
            if (this->mOutputCellProliferativeTypes)
            {
                cell_types[node_index] =  cell_ptr->GetCellCycleModel()->GetCellProliferativeType();
            }
            if (this->mOutputCellMutationStates)
            {
                cell_mutation_states[node_index] = cell_ptr->GetMutationState()->GetColour();
            }
            if (this->mOutputCellAges)
            {
                cell_ages[node_index] = cell_ptr->GetAge();
            }
            if (this->mOutputCellCyclePhases)
            {
                cell_cycle_phases[node_index] = cell_ptr->GetCellCycleModel()->GetCurrentCellCyclePhase();
            }
            ///\todo #1515 #2032 Add CellwiseData:  if (CellwiseData<DIM>::Instance()->IsSetUp())
        }
    }

    //Cell IDs can be used to threshold out the empty lattice sites (which have ID=-1)
    mesh_writer.AddPointData("Cell ids", cell_ids);

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

    if (this->mOutputCellMutationStates)
    {
        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
        mesh_writer.AddPointData("Cell labels", cell_labels);
    }

    /*
     * The current VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK then visualized as glyphs.
     */
    NodesOnlyMesh<DIM> temp_mesh;
    temp_mesh.ConstructNodesWithoutMesh(static_cast<PottsMesh<DIM>& >((this->mrMesh)));
    mesh_writer.WriteFilesUsingMesh(temp_mesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << time.str();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << time.str();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MultipleCaBasedCellPopulation<1>;
template class MultipleCaBasedCellPopulation<2>;
template class MultipleCaBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCaBasedCellPopulation)
