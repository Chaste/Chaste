/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "CaBasedCellPopulation.hpp"

#include <boost/scoped_array.hpp>

#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"

// Cell population writers
#include "CellMutationStatesWriter.hpp"

#include "NodesOnlyMesh.hpp"
#include "Exception.hpp"

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::Validate()
{
    NEVER_REACHED;
}

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
            assert(i<locationIndices.size());
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
        if (IsSiteAvailable(*neighbour_iter, pCell))
        {
            is_room = true;
            break;
        }
    }

    return is_room;
}

template<unsigned DIM>
CellPtr CaBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    // Get node index corresponding to the parent cell
    unsigned parent_node_index = this->GetLocationIndexUsingCell(pParentCell);

    // Get the set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mrMesh)).GetMooreNeighbouringNodeIndices(parent_node_index);
    unsigned num_neighbours = neighbouring_node_indices.size();

    // Each node must have at least one neighbour
    assert(!neighbouring_node_indices.empty());

    std::vector<double> neighbouring_node_propensities;
    std::vector<unsigned> neighbouring_node_indices_vector;

    double total_propensity = 0.0;

    for (std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
            neighbour_iter != neighbouring_node_indices.end();
         ++neighbour_iter)
    {
        neighbouring_node_indices_vector.push_back(*neighbour_iter);

        double propensity_dividing_into_neighbour = EvaluateDivisionPropensity(parent_node_index,*neighbour_iter,pParentCell);

        if (!IsSiteAvailable(*neighbour_iter, pParentCell))
        {
            propensity_dividing_into_neighbour = 0.0;
        }
        neighbouring_node_propensities.push_back(propensity_dividing_into_neighbour);
        total_propensity += propensity_dividing_into_neighbour;
    }

    assert(total_propensity>0); // if this trips the cell cant divided so need to include this in the IsSiteAvailable method

    for (unsigned i=0; i<num_neighbours; i++)
    {
        neighbouring_node_propensities[i] /= total_propensity;
    }

     // Sample random number to specify which move to make
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double random_number = p_gen->ranf();

    double total_probability = 0.0;
    unsigned daughter_node_index = UNSIGNED_UNSET;

    unsigned counter;
    for (counter=0; counter < num_neighbours; counter++)
    {
        total_probability += neighbouring_node_propensities[counter];
        if (total_probability >= random_number)
        {
            // Divide the parent cell to this neighbour location
            daughter_node_index = neighbouring_node_indices_vector[counter];
            break;
        }
    }
    // This loop should always break as sum(neighbouring_node_propensities) = 1

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
            // Get the index of the node corresponding to this cell
            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

            RemoveCellUsingLocationIndex(node_index, (*cell_iter));

            // Erase cell and update counter
            num_removed++;
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
void CaBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
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
                    for (typename std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >::iterator iterRule = mUpdateRuleCollection.begin();
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
void CaBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule)
{
    mUpdateRuleCollection.push_back(pUpdateRule);
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::RemoveAllUpdateRules()
{
    mUpdateRuleCollection.clear();
}

template<unsigned DIM>
const std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >& CaBasedCellPopulation<DIM>::rGetUpdateRuleCollection() const
{
    return mUpdateRuleCollection;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
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
    if (num_cells > 0u)
    {
        num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    // Not used here: std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();
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
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
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
        for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            // Populate the vector of VTK cell data
            vtk_cell_data[cell_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(*cell_iter, this);
            cell_index++;
        }

        mesh_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    }

    /*
     * At present, the VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK, then visualized
     * as glyphs in Paraview.
     */
    NodesOnlyMesh<DIM> temp_mesh;
    temp_mesh.ConstructNodesWithoutMesh(nodes, 1.5);  // Arbitrary cut off as connectivity not used.
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

// Explicit instantiation
template class CaBasedCellPopulation<1>;
template class CaBasedCellPopulation<2>;
template class CaBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedCellPopulation)
