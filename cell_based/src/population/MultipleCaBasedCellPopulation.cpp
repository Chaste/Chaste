/*

Copyright (c) 2005-2013, University of Oxford.
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
    NEVER_REACHED;
}

template<unsigned DIM>
MultipleCaBasedCellPopulation<DIM>::MultipleCaBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                                        std::vector<CellPtr>& rCells,
                                                        const std::vector<unsigned> locationIndices,
                                                        unsigned latticeCarryingCapacity,
                                                        bool deleteMesh,
                                                        bool validate)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh),
      mLatticeCarryingCapacity(latticeCarryingCapacity)
{
    mAvailableSpaces = std::vector<unsigned>(this->GetNumNodes(), latticeCarryingCapacity);

    // This must always be true
    assert(this->mCells.size() <= this->mrMesh.GetNumNodes()*latticeCarryingCapacity);

    if (locationIndices.empty())
    {
        EXCEPTION("No location indices being passed. Specify where cells lie before creating the cell population.");
    }
    else
    {
        // Create a set of node indices corresponding to empty sites
        for (unsigned i=0; i<locationIndices.size(); i++)
        {
            if (!IsSiteAvailable(locationIndices[i]))
            {
                EXCEPTION("One of the lattice sites has more cells than the carrying capacity. Check the initial cell locations.");
            }
            mAvailableSpaces[locationIndices[i]]--;
        }
    }
    if (validate)
    {
        EXCEPTION("There is no validation for MultipleCaBasedCellPopulation.");
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
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
std::vector<unsigned>& MultipleCaBasedCellPopulation<DIM>::rGetAvailableSpaces()
{
    return mAvailableSpaces;
}

template<unsigned DIM>
bool MultipleCaBasedCellPopulation<DIM>::IsSiteAvailable(unsigned index)
{
    return (mAvailableSpaces[index] != 0);
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

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::AddCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    if (!IsSiteAvailable(index))
    {
        EXCEPTION("No available spaces at location index " << index << ".");
    }


    mAvailableSpaces[index]--;
    AbstractCellPopulation<DIM,DIM>::AddCellUsingLocationIndex(index, pCell);
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::RemoveCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    AbstractCellPopulation<DIM,DIM>::RemoveCellUsingLocationIndex(index, pCell);

    mAvailableSpaces[index]++;

    assert(mAvailableSpaces[index] <= mLatticeCarryingCapacity);
}

template<unsigned DIM>
bool MultipleCaBasedCellPopulation<DIM>::IsRoomToDivide(CellPtr pCell)
{
    bool is_room = false;

    // Get node index corresponding to this cell
    unsigned node_index = this->GetLocationIndexUsingCell(pCell);

    // Get the set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(node_index);

    // Iterate through the neighbours to see if there are any available sites
    for (std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
         neighbour_iter != neighbouring_node_indices.end();
         ++neighbour_iter)
    {
        if (IsSiteAvailable(*neighbour_iter))
        {
            is_room = true;
            break;
        }
    }

    return is_room;
}

template<unsigned DIM>
CellPtr MultipleCaBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    // Get node index corresponding to the parent cell
    unsigned parent_node_index = this->GetLocationIndexUsingCell(pParentCell);

    // Get the set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(parent_node_index);
    unsigned num_neighbours = neighbouring_node_indices.size();

    // Each node must have at least one neighbour
    assert(!neighbouring_node_indices.empty());

    // Randomly choose one of the neighbouring node indices
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    unsigned chosen_start = p_gen->randMod(num_neighbours);
    std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
    for (unsigned i=0; i<chosen_start; i++)
    {
        neighbour_iter++;
    }

    /*
     * Iterate through the neighbours until the first available site is found.
     * Note that at least one such site is guaranteed, since this method is called
     * within AbstractCellBasedSimulation::DoCellBirth() only if IsRoomToDivide()
     * returns true for the parent cell.
     */
    unsigned daughter_node_index = UNSIGNED_UNSET;
    unsigned count = 0;
    while (count < num_neighbours)
    {
        bool is_empty_site = IsSiteAvailable(*neighbour_iter);

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

    assert(daughter_node_index != UNSIGNED_UNSET);
    assert(daughter_node_index < this->mrMesh.GetNumNodes());

    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    AddCellUsingLocationIndex(daughter_node_index,p_created_cell);

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
            // Get the index of the node corresponding to this cell
            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

            RemoveCellUsingLocationIndex(node_index, (*cell_iter));

            // Erase cell and update counter
            num_removed++;
            cell_iter = this->mCells.erase(cell_iter);
            --cell_iter;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
{
    // Iterate over cells
    ///\todo make this sweep random
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        // Loop over neighbours and calculate probability of moving (make sure all probabilities are <1)
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Find a random available neighbouring node to overwrite current site
        std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(node_index);
        std::vector<double> neighbouring_node_propensities;
        std::vector<unsigned> neighbouring_node_indices_vector;

        if (!neighbouring_node_indices.empty())
        {
            unsigned num_neighbours = neighbouring_node_indices.size();
            double probability_of_not_moving = 1.0;

            for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
                 iter != neighbouring_node_indices.end();
                 ++iter)
            {
                double probability_of_moving = 0.0;

                neighbouring_node_indices_vector.push_back(*iter);

                if (IsSiteAvailable(*iter))
                {
                    // Iterating over the update rule
                    for (typename std::vector<boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > >::iterator iterRule = mUpdateRuleCollection.begin();
                         iterRule != mUpdateRuleCollection.end();
                         ++iterRule)
                    {
                        probability_of_moving += (*iterRule)->EvaluateProbability(node_index, *iter, *this, dt, 1, *cell_iter);
                        if (probability_of_moving < 0)
                        {
                            EXCEPTION("The probability of cellular movement is smaller than zero. In order to prevent it from happening you should change your time step and parameters");
                        }

                        if (probability_of_moving > 1)
                        {
                            EXCEPTION("The probability of the cellular movement is bigger than one. In order to prevent it from happening you should change your time step and parameters");
                        }
                    }

                    probability_of_not_moving -= probability_of_moving;
                    neighbouring_node_propensities.push_back(probability_of_moving);
                }
                else
                {
                    neighbouring_node_propensities.push_back(0.0);
                }
            }
            if (probability_of_not_moving < 0)
            {
                EXCEPTION("The probability of the cell not moving is smaller than zero. In order to prevent it from happening you should change your time step and parameters");
            }

            // Sample random number to specify which move to make
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            double random_number = p_gen->ranf();

            double total_probability = 0.0;
            for (unsigned counter=0; counter<num_neighbours; counter++)
            {
                total_probability += neighbouring_node_propensities[counter];
                if (total_probability >= random_number)
                {
                    //Move the cell to this neighbour location
                    unsigned chosen_neighbour_location_index = neighbouring_node_indices_vector[counter];
                    this->MoveCellInLocationMap((*cell_iter), node_index, chosen_neighbour_location_index);
                    break;
                }
            }
            // If loop completes with total_probability < random_number then stay in the same location
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

        // Write the index of the of Node the cell is associated with.
        *mpVizLocationsFile << node_index << " ";
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
        this->GenerateCellResults(*cell_iter, cell_cycle_phase_counter);
    }

    this->WriteCellResultsToFiles(cell_cycle_phase_counter);
}

template<unsigned DIM>
double MultipleCaBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    double width = this->mrMesh.GetWidth(rDimension);
    return width;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > pUpdateRule)
{
    mUpdateRuleCollection.push_back(pUpdateRule);
}

template<unsigned DIM>
void MultipleCaBasedCellPopulation<DIM>::RemoveAllUpdateRules()
{
    mUpdateRuleCollection.clear();
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
    EXCEPTION("Cannot call GetNeighbouringNodeIndices() on a MultipleCaBasedCellPopulation, need to go through the PottsMesh instead");
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

    unsigned num_cells = this->GetNumRealCells();
    std::vector<double> cell_types(num_cells, -1.0);
    std::vector<double> cell_mutation_states(num_cells, -1.0);
    std::vector<double> cell_labels(num_cells, -1.0);
    std::vector<double> cell_ids(num_cells, -1.0);
    std::vector<double> cell_ancestors(num_cells, -1.0);
    std::vector<double> cell_ages(num_cells, -1.0);
    std::vector<double> cell_cycle_phases(num_cells, -1.0);
    std::vector<Node<DIM>*> nodes;

    std::vector<std::vector<double> > cellwise_data;

    unsigned num_cell_data_items = 0;
    std::vector<std::string> cell_data_names;

    // We assume that the first cell is representative of all cells
    num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    cell_data_names = this->Begin()->GetCellData()->GetKeys();

    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cellwise_data_var(num_cells);
        cellwise_data.push_back(cellwise_data_var);
    }


    unsigned cell = 0;

    // Counter to keep track of how many cells are at a lattice site
    unsigned num_sites = this->mrMesh.GetNumNodes();
    unsigned number_of_cells_at_site[num_sites];
    for (unsigned i=0; i<num_sites; i++)
    {
        number_of_cells_at_site[i] = 0;
    }

    for (std::list<CellPtr>::iterator iter = this->mCells.begin();
         iter != this->mCells.end();
         ++iter)
    {
        CellPtr cell_ptr = *iter;
        cell_ids[cell] = cell_ptr->GetCellId();

        unsigned location_index = this->GetLocationIndexUsingCell(*iter);

        number_of_cells_at_site[location_index]++;
        assert(number_of_cells_at_site[location_index]<=mLatticeCarryingCapacity);

        c_vector<double, DIM> coords = this->mrMesh.GetNode(location_index)->rGetLocation();

        // Move the coordinate slightly so that we can visualise all cells in a lattice site if there is more than one per site
        if (mLatticeCarryingCapacity > 1)
        {
            c_vector<double, DIM> offset;

            if (DIM == 2)
            {
                double angle = (double)number_of_cells_at_site[location_index]*2.0*M_PI/(double)mLatticeCarryingCapacity;
                offset[0] = 0.2*sin(angle);
                offset[1] = 0.2*cos(angle);
            }
            else
            {
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

                for (unsigned i=0; i<DIM; i++)
                {
                    offset[i] = p_gen->ranf(); // This assumes that all sites are 1 apart
                }
            }

            for (unsigned i=0; i<DIM; i++)
            {
                coords[i] += offset[i];
            }

        }

        nodes.push_back(new Node<DIM>(cell, coords, false));

        if (this->mOutputCellAncestors)
        {
            double ancestor_index = (cell_ptr->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_ptr->GetAncestor();
            cell_ancestors[cell] = ancestor_index;
        }
        if (this->mOutputCellProliferativeTypes)
        {
            cell_types[cell] = cell_ptr->GetCellProliferativeType()->GetColour();
        }
        if (this->mOutputCellMutationStates)
        {
            cell_mutation_states[cell] = cell_ptr->GetMutationState()->GetColour();
        }
        if (this->mOutputCellAges)
        {
            cell_ages[cell] = cell_ptr->GetAge();
        }
        if (this->mOutputCellCyclePhases)
        {
            cell_cycle_phases[cell] = cell_ptr->GetCellCycleModel()->GetCurrentCellCyclePhase();
        }
        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cellwise_data[var][cell] = cell_ptr->GetCellData()->GetItem(cell_data_names[var]);
        }

        cell ++;

    }

    // Cell IDs can be used to threshold out the empty lattice sites (which have ID=-1)
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
    if (num_cell_data_items > 0)
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            mesh_writer.AddPointData(cell_data_names[var], cellwise_data[var]);
        }
    }

    /*
     * The current VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK then visualized as glyphs.
     */
    NodesOnlyMesh<DIM> temp_mesh;
    temp_mesh.ConstructNodesWithoutMesh(nodes, 1.5);  // Arbitrary cut off as connectivity not used.
    mesh_writer.WriteFilesUsingMesh(temp_mesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << time.str();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << time.str();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";

    // Tidy up
    for (unsigned i=0; i<nodes.size(); i++)
    {
        delete nodes[i];
    }
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
