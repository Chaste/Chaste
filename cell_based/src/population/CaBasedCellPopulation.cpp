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

#include <boost/scoped_array.hpp>

#include "CaBasedCellPopulation.hpp"
#include "MutableMesh.hpp"
#include "AbstractCaUpdateRule.hpp"
#include "AbstractCaSwitchingUpdateRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellLocationIndexWriter.hpp"
#include "ExclusionCaBasedDivisionRule.hpp"
#include "NodesOnlyMesh.hpp"
#include "ApoptoticCellProperty.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"

// LCOV_EXCL_START
template<unsigned DIM>
void CaBasedCellPopulation<DIM>::Validate()
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
CaBasedCellPopulation<DIM>::CaBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                                        std::vector<CellPtr>& rCells,
                                                        const std::vector<unsigned> locationIndices,
                                                        unsigned latticeCarryingCapacity,
                                                        bool deleteMesh,
                                                        bool validate)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh),
      mLatticeCarryingCapacity(latticeCarryingCapacity)
{
    mAvailableSpaces = std::vector<unsigned>(this->GetNumNodes(), latticeCarryingCapacity);
    mpCaBasedDivisionRule.reset(new ExclusionCaBasedDivisionRule<DIM>());

    // This must always be true
    assert(this->mCells.size() <= this->mrMesh.GetNumNodes()*latticeCarryingCapacity);

    if (locationIndices.empty())
    {
        EXCEPTION("No location indices being passed. Specify where cells lie before creating the cell population.");
    }
    else
    {
        // Create a set of node indices corresponding to empty sites.
        // Note iterating over mCells is OK as it has the same order as location indices at this point (its just coppied from rCells)
        std::list<CellPtr>::iterator it = this->mCells.begin();
        for (unsigned i=0; it != this->mCells.end(); ++it, ++i)
        {
            assert(i < locationIndices.size());
            if (!IsSiteAvailable(locationIndices[i],*it))
            {
                EXCEPTION("One of the lattice sites has more cells than the carrying capacity. Check the initial cell locations.");
            }
            mAvailableSpaces[locationIndices[i]]--;
        }
    }
    if (validate)
    {
        EXCEPTION("There is no validation for CaBasedCellPopulation.");
    }
}

template<unsigned DIM>
CaBasedCellPopulation<DIM>::CaBasedCellPopulation(PottsMesh<DIM>& rMesh)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
CaBasedCellPopulation<DIM>::~CaBasedCellPopulation()
{
    if (this->mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
std::vector<unsigned>& CaBasedCellPopulation<DIM>::rGetAvailableSpaces()
{
    return mAvailableSpaces;
}

template<unsigned DIM>
bool CaBasedCellPopulation<DIM>::IsSiteAvailable(unsigned index, CellPtr pCell)
{
    ///\todo this is where to deal with carrying capacity
    return (mAvailableSpaces[index] != 0);
}

template<unsigned DIM>
PottsMesh<DIM>& CaBasedCellPopulation<DIM>::rGetMesh()
{
    return static_cast<PottsMesh<DIM>& >((this->mrMesh));
}

template<unsigned DIM>
const PottsMesh<DIM>& CaBasedCellPopulation<DIM>::rGetMesh() const
{
    return static_cast<PottsMesh<DIM>& >((this->mrMesh));
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* CaBasedCellPopulation<DIM>::GetTetrahedralMeshForPdeModifier()
{
    std::vector<Node<DIM>*> temp_nodes;

    // Create nodes at the centre of the cells
    unsigned cell_index = 0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        temp_nodes.push_back(new Node<DIM>(cell_index, this->GetLocationOfCellCentre(*cell_iter)));
        cell_index++;
    }

    return new MutableMesh<DIM,DIM>(temp_nodes);
}

template<unsigned DIM>
Node<DIM>* CaBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned CaBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
std::set<unsigned> CaBasedCellPopulation<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    unsigned index = this->GetLocationIndexUsingCell(pCell);
    std::set<unsigned> candidates = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(index);

    std::set<unsigned> neighbour_indices;
    for (std::set<unsigned>::iterator iter = candidates.begin();
         iter != candidates.end();
         ++iter)
    {
        if (!IsSiteAvailable(*iter, pCell))
        {
            neighbour_indices.insert(*iter);
        }
    }

    return neighbour_indices;
}

template<unsigned DIM>
c_vector<double, DIM> CaBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return this->mrMesh.GetNode(this->GetLocationIndexUsingCell(pCell))->rGetLocation();
}

template<unsigned DIM>
Node<DIM>* CaBasedCellPopulation<DIM>::GetNodeCorrespondingToCell(CellPtr pCell)
{
    return this->mrMesh.GetNode(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::AddCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    if (!IsSiteAvailable(index, pCell))
    {
        EXCEPTION("No available spaces at location index " << index << ".");
    }

    mAvailableSpaces[index]--;
    AbstractCellPopulation<DIM,DIM>::AddCellUsingLocationIndex(index, pCell);
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::RemoveCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    AbstractCellPopulation<DIM,DIM>::RemoveCellUsingLocationIndex(index, pCell);

    mAvailableSpaces[index]++;

    assert(mAvailableSpaces[index] <= mLatticeCarryingCapacity);
}

template<unsigned DIM>
bool CaBasedCellPopulation<DIM>::IsRoomToDivide(CellPtr pCell)
{
    return mpCaBasedDivisionRule->IsRoomToDivide(pCell, *this);
}

template<unsigned DIM>
CellPtr CaBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    unsigned daughter_node_index = mpCaBasedDivisionRule->CalculateDaughterNodeIndex(pNewCell,pParentCell,*this);

    // Associate the new cell with the neighbouring node
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    AddCellUsingLocationIndex(daughter_node_index,p_created_cell);

    return p_created_cell;
}

template<unsigned DIM>
double CaBasedCellPopulation<DIM>:: EvaluateDivisionPropensity(unsigned currentNodeIndex,
                                                               unsigned targetNodeIndex,
                                                               CellPtr pCell)
{
    return 1.0;
}

template<unsigned DIM>
unsigned CaBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         )
    {
        if ((*cell_iter)->IsDead())
        {
            // Get the location index corresponding to this cell
            unsigned location_index = this->GetLocationIndexUsingCell(*cell_iter);

            // Use this to remove the cell from the population
            RemoveCellUsingLocationIndex(location_index, (*cell_iter));

            // Erase cell and update counter
            cell_iter = this->mCells.erase(cell_iter);
            num_removed++;
        }
        else
        {
            ++cell_iter;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
{
    /*
     * Here we loop over the nodes and calculate the probability of moving
     * and then select the node to move to.
     */
    if (!(this->mUpdateRuleCollection.empty()))
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

                    if (IsSiteAvailable(*iter, *cell_iter))
                    {
                        // Iterating over the update rule
                        for (typename std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > >::iterator iter_rule = this->mUpdateRuleCollection.begin();
                             iter_rule != this->mUpdateRuleCollection.end();
                             ++iter_rule)
                        {
                            // This static cast is fine, since we assert the update rule must be a CA update rule in AddUpdateRule()
                            double p = (boost::static_pointer_cast<AbstractCaUpdateRule<DIM> >(*iter_rule))->EvaluateProbability(node_index, *iter, *this, dt, 1, *cell_iter);
                            probability_of_moving += p;
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
                        // Move the cell to this neighbour location
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

    /*
     * Here we loop over the nodes and select a neighbour to test if
     * the cells (associated with the nodes) should swap locations
     * or if a cell should move to an empty node
     * Note this currently only works for latticeCarryingCapacity = 1
     */
    if (!(mSwitchingUpdateRuleCollection.empty()))
    {
        assert(mLatticeCarryingCapacity == 1);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        unsigned num_nodes = this->mrMesh.GetNumNodes();

        // Randomly permute mUpdateRuleCollection if specified
        if (this->mIterateRandomlyOverUpdateRuleCollection)
        {
            // Randomly permute mUpdateRuleCollection
            p_gen->Shuffle(mSwitchingUpdateRuleCollection);
        }

        for (unsigned i=0; i<num_nodes; i++)
        {
            unsigned node_index;

            if (this->mUpdateNodesInRandomOrder)
            {
                node_index = p_gen->randMod(num_nodes);
            }
            else
            {
                // Loop over nodes in index order
                node_index = i%num_nodes;
            }

            // Find a random available neighbouring node to switch cells with the current site
            std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(node_index);

            unsigned neighbour_location_index;

            if (!neighbouring_node_indices.empty())
            {
                unsigned num_neighbours = neighbouring_node_indices.size();
                unsigned chosen_neighbour = p_gen->randMod(num_neighbours);

                std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
                for (unsigned j=0; j<chosen_neighbour; j++)
                {
                    neighbour_iter++;
                }
                neighbour_location_index = *neighbour_iter;

                bool is_cell_on_node_index = mAvailableSpaces[node_index] == 0 ? true : false;
                bool is_cell_on_neighbour_location_index = mAvailableSpaces[neighbour_location_index] == 0 ? true : false;

                if (is_cell_on_node_index || is_cell_on_neighbour_location_index)
                {
                    double probability_of_switch = 0.0;

                    // Now add contributions to the probability from each CA switching update rule
                    for (typename std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > >::iterator iter_rule = mSwitchingUpdateRuleCollection.begin();
                         iter_rule != mSwitchingUpdateRuleCollection.end();
                         ++iter_rule)
                    {
                        // This static cast is fine, since we assert the update rule must be a CA switching update rule in AddUpdateRule()
                        double p = (boost::static_pointer_cast<AbstractCaSwitchingUpdateRule<DIM> >(*iter_rule))->EvaluateSwitchingProbability(node_index, neighbour_location_index, *this, dt, 1);
                        probability_of_switch += p;
                    }

                    assert(probability_of_switch >= 0);
                    assert(probability_of_switch <= 1);

                    // Generate a uniform random number to do the random switch
                    double random_number = p_gen->ranf();

                    if (random_number < probability_of_switch)
                    {
                        if (is_cell_on_node_index && is_cell_on_neighbour_location_index)
                        {
                            // Swap the cells associated with the node and the neighbour node
                            CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);
                            CellPtr p_neighbour_cell = this->GetCellUsingLocationIndex(neighbour_location_index);

                            // Remove the cells from their current location
                            RemoveCellUsingLocationIndex(node_index, p_cell);
                            RemoveCellUsingLocationIndex(neighbour_location_index, p_neighbour_cell);

                            // Add cells to their new locations
                            AddCellUsingLocationIndex(node_index, p_neighbour_cell);
                            AddCellUsingLocationIndex(neighbour_location_index, p_cell);
                        }
                        else if (is_cell_on_node_index && !is_cell_on_neighbour_location_index)
                        {
                            // Move the cells associated with the node to the neighbour node
                            CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);
                            RemoveCellUsingLocationIndex(node_index, p_cell);
                            AddCellUsingLocationIndex(neighbour_location_index, p_cell);
                        }
                        else if (!is_cell_on_node_index && is_cell_on_neighbour_location_index)
                        {
                            // Move the cell associated with the neighbour node onto the node
                            CellPtr p_neighbour_cell = this->GetCellUsingLocationIndex(neighbour_location_index);
                            RemoveCellUsingLocationIndex(neighbour_location_index, p_neighbour_cell);
                            AddCellUsingLocationIndex(node_index, p_neighbour_cell);
                        }
                        else
                        {
                            NEVER_REACHED;
                        }
                    }
                }
            }
        }
    }
}

template<unsigned DIM>
bool CaBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return false;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter)
{
    pPopulationCountWriter->Visit(this);
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellLocationIndexWriter>())
        {
            this-> template AddCellWriter<CellLocationIndexWriter>();
        }
    }

    AbstractCellPopulation<DIM>::OpenWritersFiles(rOutputFileHandler);
}

template<unsigned DIM>
double CaBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    double cell_volume = 1.0;
    return cell_volume;
}

template<unsigned DIM>
double CaBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    double width = this->mrMesh.GetWidth(rDimension);
    return width;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule)
{
    // The update rule must be derived from AbstractCaUpdateRule or AbstractCaSwitchingUpdateRule
    assert(bool(dynamic_cast<AbstractCaUpdateRule<DIM>*>(pUpdateRule.get())) ||
           bool(dynamic_cast<AbstractCaSwitchingUpdateRule<DIM>*>(pUpdateRule.get())));

    if (bool(dynamic_cast<AbstractCaUpdateRule<DIM>*>(pUpdateRule.get())))
    {
        this->mUpdateRuleCollection.push_back(pUpdateRule);
    }
    else
    {
        mSwitchingUpdateRuleCollection.push_back(pUpdateRule);
    }
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::RemoveAllUpdateRules()
{
    // Clear mSwitchingUpdateRuleCollection
    mSwitchingUpdateRuleCollection.clear();

    // Clear mUpdateRuleCollection
    AbstractOnLatticeCellPopulation<DIM>::RemoveAllUpdateRules();
}

template<unsigned DIM>
const std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > CaBasedCellPopulation<DIM>::GetUpdateRuleCollection() const
{
    std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > update_rules;

    for (unsigned i=0; i<this->mUpdateRuleCollection.size(); i++)
    {
        update_rules.push_back(this->mUpdateRuleCollection[i]);
    }
    for (unsigned i=0; i<mSwitchingUpdateRuleCollection.size(); i++)
    {
        update_rules.push_back(mSwitchingUpdateRuleCollection[i]);
    }

    return update_rules;
}

template<unsigned DIM>
boost::shared_ptr<AbstractCaBasedDivisionRule<DIM> > CaBasedCellPopulation<DIM>::GetCaBasedDivisionRule()
{
    return mpCaBasedDivisionRule;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::SetCaBasedDivisionRule(boost::shared_ptr<AbstractCaBasedDivisionRule<DIM> > pCaBasedDivisionRule)
{
    mpCaBasedDivisionRule = pCaBasedDivisionRule;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Add the division rule parameters
    *rParamsFile << "\t\t<CaBasedDivisionRule>\n";
    mpCaBasedDivisionRule->OutputCellCaBasedDivisionRuleInfo(rParamsFile);
    *rParamsFile << "\t\t</CaBasedDivisionRule>\n";

    // Call method on direct parent class
    AbstractOnLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    // Store the present time as a string
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    // Store the number of cells for which to output data to VTK
    unsigned num_cells = this->GetNumRealCells();

    // When outputting any CellData, we assume that the first cell is representative of all cells
    unsigned num_cell_data_items = 0u;
    std::vector<std::string> cell_data_names;
    if (num_cells > 0u)
    {
        num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
        cell_data_names = this->Begin()->GetCellData()->GetKeys();
    }
    std::vector<std::vector<double> > cell_data;
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cell_data_var(num_cells);
        cell_data.push_back(cell_data_var);
    }

    // Create mesh writer for VTK output
    VtkMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results_"+time.str(), false);

    // Create a counter to keep track of how many cells are at a lattice site
    unsigned num_sites = this->mrMesh.GetNumNodes();
    boost::scoped_array<unsigned> number_of_cells_at_site(new unsigned[num_sites]);
    for (unsigned i=0; i<num_sites; i++)
    {
        number_of_cells_at_site[i] = 0;
    }

    // Populate a vector of nodes associated with cell locations, by iterating through the list of cells
    std::vector<Node<DIM>*> nodes;
    unsigned node_index = 0;
    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the location index of this cell and update the counter number_of_cells_at_site
        unsigned location_index = this->GetLocationIndexUsingCell(*cell_iter);
        number_of_cells_at_site[location_index]++;
        assert(number_of_cells_at_site[location_index] <= mLatticeCarryingCapacity);

        // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
        c_vector<double, DIM> coords;
        coords = this->mrMesh.GetNode(location_index)->rGetLocation();

        // Move the coordinates slightly so that we can visualise all cells in a lattice site if there is more than one per site
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
                    offset[i] = p_gen->ranf(); // This assumes that all sites are 1 unit apart
                }
            }

            for (unsigned i=0; i<DIM; i++)
            {
                coords[i] += offset[i];
            }
        }

        nodes.push_back(new Node<DIM>(node_index, coords, false));
        node_index++;
    }

    // Iterate over any cell writers that are present
    for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
         cell_writer_iter != this->mCellWriters.end();
         ++cell_writer_iter)
    {
        // Create vector to store VTK cell data
        // (using a default value of -1 to correspond to an empty lattice site)
        std::vector<double> vtk_cell_data(num_cells, -1.0);

        // Loop over cells
        unsigned cell_index = 0;
        for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin();
             cell_iter != this->End();
             ++cell_iter)
        {
            // Populate the vector of VTK cell data
            vtk_cell_data[cell_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(*cell_iter, this);
            cell_index++;
        }

        mesh_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    }

    // Loop over cells to output the cell data
    unsigned cell_index = 0;
    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cell_data[var][cell_index] = cell_iter->GetCellData()->GetItem(cell_data_names[var]);
        }
        cell_index++;
    }
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        mesh_writer.AddPointData(cell_data_names[var], cell_data[var]);
    }

    /*
     * At present, the VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK, then visualized
     * as glyphs in Paraview.
     */
    NodesOnlyMesh<DIM> temp_mesh;

    // Use an approximation of the node spacing as the interaction distance for the nodes only mesh. This is to
    // avoid rounding errors in distributed box collection.
    double volume = this->mrMesh.GetWidth(0);
    for (unsigned idx=1; idx<DIM; idx++)
    {
        volume *= this->mrMesh.GetWidth(idx);
    }

    double spacing;
    if (this->mrMesh.GetNumNodes() >0 && volume > 0.0)
    {
        spacing = std::pow(volume / double(this->mrMesh.GetNumNodes()), 1.0/double(DIM));
    }
    else
    {
        spacing = 1.0;
    }

    temp_mesh.ConstructNodesWithoutMesh(nodes, spacing * 1.2);
    mesh_writer.WriteFilesUsingMesh(temp_mesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";

    // Tidy up
    for (unsigned i=0; i<nodes.size(); i++)
    {
        delete nodes[i];
    }
#endif //CHASTE_VTK
}

template<unsigned DIM>
double CaBasedCellPopulation<DIM>::GetCellDataItemAtPdeNode(
    unsigned pdeNodeIndex,
    std::string& rVariableName,
    bool dirichletBoundaryConditionApplies,
    double dirichletBoundaryValue)
{
    unsigned counter = 0;
    typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
    while (counter != pdeNodeIndex)
    {
        ++cell_iter;
        counter++;
    }

    double value = cell_iter->GetCellData()->GetItem(rVariableName);

    return value;
}

template<unsigned DIM>
bool CaBasedCellPopulation<DIM>::IsPdeNodeAssociatedWithNonApoptoticCell(unsigned pdeNodeIndex)
{
    // pdeNodeIndex corresponds to the 'position' of the cell to interrogate in the vector of cells
    typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();

    assert(pdeNodeIndex < this->GetNumRealCells());
    for (unsigned i=0; i<pdeNodeIndex; i++)
    {
        ++cell_iter;
    }
    bool is_cell_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();

    return !is_cell_apoptotic;
}

// Explicit instantiation
template class CaBasedCellPopulation<1>;
template class CaBasedCellPopulation<2>;
template class CaBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedCellPopulation)
