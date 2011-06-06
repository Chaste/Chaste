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

#include "LatticeBasedCellPopulation.hpp"

#include <cassert>
#include <algorithm>

#include "Exception.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"


template<unsigned DIM>
LatticeBasedCellPopulation<DIM>::LatticeBasedCellPopulation(TetrahedralMesh<DIM, DIM>& rMesh,
                                            std::vector<CellPtr>& rCells,
                                            const std::vector<unsigned> locationIndices,
                                            bool onlyUseNearestNeighboursForDivision,
                                            bool useVonNeumannNeighbourhoods,
                                            bool deleteMesh,
                                            bool validate)
    : AbstractCellPopulation<DIM>(rCells, locationIndices),
      mrMesh(rMesh),
      mOnlyUseNearestNeighboursForDivision(onlyUseNearestNeighboursForDivision),
      mDeleteMesh(deleteMesh),
      mUseVonNeumannNeighbourhoods(useVonNeumannNeighbourhoods)
{
    // This must always be true
    assert(this->mCells.size() <= mrMesh.GetNumNodes());

    if (!locationIndices.empty())
    {
        // Create a set of node indices corresponding to empty sites
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
        SetEmptySites(ghost_node_indices);
    }
    else
    {
        this->mIsEmptySite = std::vector<bool>(this->GetNumNodes(), false);
        Validate();
    }
}

template<unsigned DIM>
LatticeBasedCellPopulation<DIM>::LatticeBasedCellPopulation(TetrahedralMesh<DIM, DIM>& rMesh)
    : mrMesh(rMesh)
{
    mDeleteMesh = true;
}

template<unsigned DIM>
LatticeBasedCellPopulation<DIM>::~LatticeBasedCellPopulation()
{
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::SetVonNeumannNeighbourhoods(bool useVonNeumannNeighbourhoods)
{
    mUseVonNeumannNeighbourhoods = useVonNeumannNeighbourhoods;
}


template<unsigned DIM>
std::vector<bool>& LatticeBasedCellPopulation<DIM>::rGetEmptySites()
{
    return this->mIsEmptySite;
}

template<unsigned DIM>
bool LatticeBasedCellPopulation<DIM>::IsEmptySite(unsigned index)
{
    return this->mIsEmptySite[index];
}

template<unsigned DIM>
std::set<unsigned> LatticeBasedCellPopulation<DIM>::GetEmptySiteIndices()
{
    std::set<unsigned> ghost_node_indices;
    for (unsigned i=0; i<this->mIsEmptySite.size(); i++)
    {
        if (this->mIsEmptySite[i])
        {
            ghost_node_indices.insert(i);
        }
    }
    return ghost_node_indices;
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::SetEmptySites(const std::set<unsigned>& rEmptySiteIndices)
{
    // Reinitialise all entries of mIsEmptySite to false
    this->mIsEmptySite = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsEmptySite
    for (std::set<unsigned>::iterator iter=rEmptySiteIndices.begin(); iter!=rEmptySiteIndices.end(); ++iter)
    {
        this->mIsEmptySite[*iter] = true;
    }

    Validate();
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>& LatticeBasedCellPopulation<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
const TetrahedralMesh<DIM, DIM>& LatticeBasedCellPopulation<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
Node<DIM>* LatticeBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned LatticeBasedCellPopulation<DIM>::GetNumNodes()
{
    return mrMesh.GetNumAllNodes();
}

template<unsigned DIM>
CellPtr LatticeBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    ///\todo This method could probably be made more efficient (#1411)

    // Add rNewCell to mCells
    this->mCells.push_back(pNewCell);
    CellPtr p_created_cell = this->mCells.back();

    // Temporarily provide a location index for the new cell
    this->mCellLocationMap[p_created_cell.get()] = UINT_MAX;

    // Get the location index corresponding to the parent cell
    unsigned parent_index = this->mCellLocationMap[pParentCell.get()];

    unsigned degree_upper_bound;         // Maximum degree which to explore to

    if (mOnlyUseNearestNeighboursForDivision)
    {
        degree_upper_bound = 1;            // Only want to search the next nearest neighbours
    }
    else
    {
        // Get the maximum neighbour degree of this location index in each direction
        std::vector<unsigned> maximum_degrees = GetMaximumDegreeInEachDirection(parent_index);

        // Find the largest maximum degree, giving an upper bound on our search for a free location index for the new cell
        degree_upper_bound = *(std::max_element(maximum_degrees.begin(), maximum_degrees.end()));
    }

    bool free_neighbour_was_found = false;

    // Loop over increasing neighbour degrees
    for (unsigned degree=1; degree<=degree_upper_bound; degree++)
    {
        // Get the set of nth degree neighbour indices
        std::set<unsigned> nth_degree_neighbours = GetNthDegreeNeighbouringNodeIndices(parent_index, degree);

        // Now find which (if any) of these neighbours are free (i.e. empty sites)
        std::set<unsigned> free_nth_degree_neighbours;
        for (std::set<unsigned>::iterator neighbour_iter = nth_degree_neighbours.begin();
             neighbour_iter != nth_degree_neighbours.end();
             ++neighbour_iter)
        {
            if (IsEmptySite(*neighbour_iter))
            {
                free_nth_degree_neighbours.insert(*neighbour_iter);
            }
        }

        // If there are any free neighbours at this degree...
        if (!free_nth_degree_neighbours.empty())
        {
            // ... then choose a free neighbour at random
            unsigned num_free_neighbours = free_nth_degree_neighbours.size();
            unsigned random_offset = RandomNumberGenerator::Instance()->randMod(num_free_neighbours);

            std::set<unsigned>::iterator free_neighbour_iter = free_nth_degree_neighbours.begin();
            for (unsigned i=0; i<random_offset; i++)
            {
                free_neighbour_iter++;
            }

            unsigned chosen_free_neighbour_index = *free_neighbour_iter;

            // Compute the direction from the parent cell to the chosen free nth degree neighbour
            int direction_index = ((int)chosen_free_neighbour_index - (int)parent_index)/(int)degree;

            // Move each cell between the parent cell and the chosen free nth degree neighbour
            std::vector<unsigned> indices(degree);
            for (unsigned i=0; i<degree; i++)
            {
                indices[i] = parent_index + (i+1)*direction_index;
            }
            for (unsigned i=degree-1; i>0; i--)
            {
                CellPtr p_current_cell = this->mLocationCellMap[indices[i-1]];
                assert(p_current_cell);

                MoveCell(p_current_cell, indices[i]);
            }

            // This creates a free nearest neighbour of the parent cell, which we associate with the new cell
            MoveCell(p_created_cell, parent_index+direction_index);

            // Break out of the loop
            free_neighbour_was_found = true;
            break;
        }
    }

    // If no free neighbour was found, throw an exception
    if (!free_neighbour_was_found)
    {
        EXCEPTION("Cell can not divide as there are no free neighbours at maximum degree in any direction");
    }

    return p_created_cell;
}

template<unsigned DIM>
std::vector<unsigned> LatticeBasedCellPopulation<DIM>::GetMaximumDegreeInEachDirection(unsigned nodeIndex)
{
    double width = this->mrMesh.GetWidth(0);
    unsigned nodes_across = (unsigned)width + 1;
    double height = this->mrMesh.GetWidth(1);
    unsigned nodes_up = (unsigned)height + 1;

    /*
     * Create a temporary vector to store the degrees in each direction,
     * some of which may be zero (if the node is on the boundary). We
     * take care to use the same ordering as is used in the method
     * GetNeighbouringNodeIndicesVector(), namely N, NW, W, SW, S, SE,
     * E, NE.
     */

    unsigned degrees_south = nodeIndex/nodes_up;
    unsigned degrees_north = nodes_up - degrees_south - 1;
    unsigned degrees_west = (nodeIndex < nodes_across) ? nodeIndex : nodeIndex%nodes_across;
    unsigned degrees_east = nodes_across - degrees_west - 1;

    std::vector<unsigned> degrees(8);
    degrees[0] = degrees_north; // N
    degrees[1] = std::min(degrees_north, degrees_west); // NW
    degrees[2] = degrees_west; // W
    degrees[3] = std::min(degrees_south, degrees_west); // SW
    degrees[4] = degrees_south; // S
    degrees[5] = std::min(degrees_south, degrees_east); // SE
    degrees[6] = degrees_east; // E
    degrees[7] = std::min(degrees_north, degrees_east); // NE

    /*
     * Now use this vector to populate another vector, which stores
     * only the non-zero degrees (we need to do this because this
     * method is used alongside GetNeighbouringNodeIndicesVector()
     * in the method GetNthDegreeNeighbouringNodeIndices(), so should
     * be the same length.
     */

    std::vector<unsigned> non_zero_degrees_in_each_direction;
    for (unsigned i=0; i<8; i++)
    {
        if (degrees[i] > 0)
        {
            non_zero_degrees_in_each_direction.push_back(degrees[i]);
        }
    }

    return non_zero_degrees_in_each_direction;
}


template<unsigned DIM>
std::set<unsigned> LatticeBasedCellPopulation<DIM>::GetNthDegreeNeighbouringNodeIndices(unsigned nodeIndex, unsigned degree)
{
    std::set<unsigned> nth_degree_neighbours;
    std::vector<unsigned> nearest_neighbours = this->GetNeighbouringNodeIndicesVector(nodeIndex);

    int next_neighbour;

    std::vector<unsigned> max_degrees = GetMaximumDegreeInEachDirection(nodeIndex);

    // Loop over all directions
    unsigned current_node;

    for (unsigned i=0; i<nearest_neighbours.size(); i++)
    {
        current_node = nearest_neighbours[i];

        // Find the n-th degree next nearest neighbour in the current direction
        if (degree <= max_degrees[i])
        {
            next_neighbour = (int) nodeIndex + (int) degree * ((int) current_node - (int) nodeIndex);
            nth_degree_neighbours.insert(next_neighbour);
        }
    }
    return nth_degree_neighbours;
}


template<unsigned DIM>
std::set<unsigned> LatticeBasedCellPopulation<DIM>::GetFreeNeighbouringNodeIndices(unsigned nodeIndex)
{
    std::set<unsigned> free_neighbouring_nodes;

    // Get all neighbouring node indices
    std::set<unsigned> all_neighbours = GetNeighbouringNodeIndices(nodeIndex);

    // Now find which neighbours are free (i.e. empty sites)
    for (std::set<unsigned>::iterator neighbour_iter = all_neighbours.begin();
         neighbour_iter != all_neighbours.end();
         ++neighbour_iter)
    {
        if (IsEmptySite(*neighbour_iter))
        {
            free_neighbouring_nodes.insert(*neighbour_iter);
        }
    }
    return free_neighbouring_nodes;
}


template<unsigned DIM>
std::set<unsigned> LatticeBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Get the neighbouring node indices as a vector
    std::vector<unsigned> neighbouring_node_indices_vector = GetNeighbouringNodeIndicesVector(nodeIndex);

    // Now populate a std::set with these node indices
    std::set<unsigned> all_neighbours;
    for (unsigned i=0; i<neighbouring_node_indices_vector.size(); i++)
    {
        all_neighbours.insert(neighbouring_node_indices_vector[i]);
    }

    return all_neighbours;
}


template<unsigned DIM>
std::vector<unsigned> LatticeBasedCellPopulation<DIM>::GetNeighbouringNodeIndicesVector(unsigned nodeIndex)
{
    std::vector<unsigned> all_neighbours;

    switch (DIM)
    {
        case 1:
        {
            double width = this->mrMesh.GetWidth(0);
            unsigned nodes_across = (unsigned)width + 1; // Getting the number of nodes along x axis of mesh

            /*
             * This stores the available neighbours using the following numbering:
             *
             * 0------x------1
             *
             */

            // Create a vector of possible neighbouring node indices
            std::vector<unsigned> neighbour_indices(2, nodeIndex);
            neighbour_indices[0] -= 1;
            neighbour_indices[1] += 1;

            // Work out whether this node lies on any edge of the mesh
            bool on_west_edge = (nodeIndex == 0);
            bool on_east_edge = (nodeIndex == nodes_across-1);

            // Create a vector of booleans for which neighbours are available
            std::vector<bool> available_neighbours = std::vector<bool>(8, true);
            available_neighbours[0] = !on_west_edge;
            available_neighbours[1] = !on_east_edge;

            // Using neighbour_indices and available_neighbours, store the indices of all available neighbours to the set all_neighbours
            for (unsigned i=0; i<2; i++)
            {
                if (available_neighbours[i])
                {
                    all_neighbours.push_back(neighbour_indices[i]);
                }
            }
            break;
        }
        case 2:
        {
            double width = this->mrMesh.GetWidth(0);
            unsigned nodes_across = (unsigned)width + 1; // Getting the number of nodes along x axis of mesh

            double height = this->mrMesh.GetWidth(1);
            unsigned nodes_up = (unsigned)height + 1;

            /*
             * This stores the available neighbours using the following numbering:
             *
             *  1-----0-----7
             *  |     |     |
             *  |     |     |
             *  2-----x-----6
             *  |     |     |
             *  |     |     |
             *  3-----4-----5
             */

            // Create a vector of possible neighbouring node indices
            std::vector<unsigned> neighbour_indices(8, nodeIndex);
            neighbour_indices[0] += nodes_across;
            neighbour_indices[1] += nodes_across - 1;
            neighbour_indices[2] -= 1;
            neighbour_indices[3] -= nodes_across + 1;
            neighbour_indices[4] -= nodes_across;
            neighbour_indices[5] -= nodes_across - 1;
            neighbour_indices[6] += 1;
            neighbour_indices[7] += nodes_across + 1;

            // Work out whether this node lies on any edge of the mesh
            bool on_south_edge = (nodeIndex < nodes_across);
            bool on_north_edge = (nodeIndex > nodes_up*(nodes_across - 1)-1);
            bool on_west_edge = (nodeIndex%nodes_across == 0);
            bool on_east_edge = (nodeIndex%nodes_across == nodes_across - 1);

            // Create a vector of booleans for which neighbours are available
            // Use the order N, NW, W, SW, S, SE, E, NE
            std::vector<bool> available_neighbours = std::vector<bool>(8, true);
            available_neighbours[0] = !on_north_edge;
            available_neighbours[1] = !(on_north_edge || on_west_edge || mUseVonNeumannNeighbourhoods);
            available_neighbours[2] = !on_west_edge;
            available_neighbours[3] = !(on_south_edge || on_west_edge || mUseVonNeumannNeighbourhoods);
            available_neighbours[4] = !on_south_edge;
            available_neighbours[5] = !(on_south_edge || on_east_edge || mUseVonNeumannNeighbourhoods);
            available_neighbours[6] = !on_east_edge;
            available_neighbours[7] = !(on_north_edge || on_east_edge || mUseVonNeumannNeighbourhoods);

            // Using neighbour_indices and available_neighbours, store the indices of all available neighbours to the set all_neighbours
            for (unsigned i=0; i<8; i++)
            {
                if (available_neighbours[i])
                {
                    all_neighbours.push_back(neighbour_indices[i]);
                }
            }
            break;
        }
        case 3:
        {
            double width = this->mrMesh.GetWidth(0);
            unsigned nodes_across = (unsigned)width + 1; // Getting the number of nodes along x axis of mesh

            double height = this->mrMesh.GetWidth(1);
            unsigned nodes_up = (unsigned)height + 1; // Getting the number of nodes along y axis of mesh

            double depth = this->mrMesh.GetWidth(2);
            unsigned nodes_deep = (unsigned)depth + 1; // Getting the number of nodes along z axis of mesh

            unsigned nodes_layer = nodes_across*nodes_up; // Number of nodes in each layer

            /*
             * This stores the available neighbours using the following numbering:
             *
             *     6------7------8           14-----15-----16            23-----24-----25
             *     |      |      |           |      |      |             |      |      |
             *     |      |      |           |      |      |             |      |      |
             *     3------4------5           12-----x------13            20-----21-----22
             *     |      |      |           |      |      |             |      |      |
             *     |      |      |           |      |      |             |      |      |
             *     0------1------2           9------10-----11            17-----18-----19
             *
             */

            // Create a vector of possible neighbouring node indices
            std::vector<unsigned> neighbour_indices(26, nodeIndex);
            for (unsigned i=0; i<9; i++)
            {
                neighbour_indices[i] -= nodes_layer;
                neighbour_indices[25-i] += nodes_layer;
            }
            for (unsigned i=0; i<3; i++)
            {
                neighbour_indices[i] -= nodes_across;
                neighbour_indices[9+i] -= nodes_across;
                neighbour_indices[17+i] -= nodes_across;
                neighbour_indices[6+i] += nodes_across;
                neighbour_indices[14+i] += nodes_across;
                neighbour_indices[23+i] += nodes_across;
            }
            for (unsigned i=0; i<4; i++)
            {
                neighbour_indices[3*i] -= 1;
                neighbour_indices[25-3*i] += 1;
                neighbour_indices[2+3*i] += 1;
                neighbour_indices[23-3*i] -= 1;
            }
            neighbour_indices[12] -= 1;
            neighbour_indices[13] += 1;

            // Work out whether this node lies on any face of the mesh
            bool on_bottom_face = (nodeIndex < nodes_layer);
            bool on_top_face = (nodeIndex > (nodes_deep-1)*nodes_layer-1);
            bool on_south_face = (nodeIndex%nodes_layer < nodes_across);
            bool on_north_face = (nodeIndex%nodes_layer > nodes_across*(nodes_up-1)-1);
            bool on_west_face = (nodeIndex%nodes_across == 0);
            bool on_east_face = (nodeIndex%nodes_across == nodes_across - 1);

            // Create a vector of booleans for which neighbours are available
            std::vector<bool> available_neighbours = std::vector<bool>(26, true);
            available_neighbours[0] = !(on_bottom_face || on_west_face || on_south_face || mUseVonNeumannNeighbourhoods);
            available_neighbours[1] = !(on_bottom_face || on_south_face || mUseVonNeumannNeighbourhoods);
            available_neighbours[2] = !(on_bottom_face || on_east_face || on_south_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[3] = !(on_bottom_face || on_west_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[4] = !on_bottom_face;
            available_neighbours[5] = !(on_bottom_face || on_east_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[6] = !(on_bottom_face || on_west_face || on_north_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[7] = !(on_bottom_face || on_north_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[8] = !(on_bottom_face || on_east_face || on_north_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[9] = !(on_west_face || on_south_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[10] = !on_south_face;
            available_neighbours[11] = !(on_east_face || on_south_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[12] = !on_west_face;
            available_neighbours[13] = !on_east_face;
            available_neighbours[14] = !(on_west_face || on_north_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[15] = !on_north_face;
            available_neighbours[16] = !(on_east_face || on_north_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[17] = !(on_top_face || on_west_face || on_south_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[18] = !(on_top_face || on_south_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[19] = !(on_top_face || on_east_face || on_south_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[20] = !(on_top_face || on_west_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[21] = !on_top_face;
            available_neighbours[22] = !(on_top_face || on_east_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[23] = !(on_top_face || on_west_face || on_north_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[24] = !(on_top_face || on_north_face|| mUseVonNeumannNeighbourhoods);
            available_neighbours[25] = !(on_top_face || on_east_face || on_north_face|| mUseVonNeumannNeighbourhoods);

            // Using neighbour_indices and available_neighbours, store the indices of all available neighbours to the set all_neighbours
            for (unsigned i=0; i<26; i++)
            {
                if (available_neighbours[i])
                {
                    all_neighbours.push_back(neighbour_indices[i]);
                }
            }
            break;
        }
        default:
            NEVER_REACHED;
    }
    return all_neighbours;
}


template<unsigned DIM>
unsigned LatticeBasedCellPopulation<DIM>::RemoveDeadCells()
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

            // Set this node to be an empty site
            this->mIsEmptySite[node_index] = true;
            this->mCellLocationMap.erase((*cell_iter).get());
            this->mLocationCellMap.erase(node_index);

            // Erase cell and update counter
            num_removed++;
            cell_iter = this->mCells.erase(cell_iter);
            --cell_iter;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    Validate();
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::Validate()
{
    // Get a list of all the nodes that are ghosts
    std::vector<bool> validated_node = mIsEmptySite;

    assert(mIsEmptySite.size()==this->GetNumNodes());

    // Look through all of the cells and record what node they are associated with.
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];

        // If the node attached to this cell is labelled as an empty site, then throw an error
        if (mIsEmptySite[node_index])
        {
            std::stringstream ss;
            ss << "Node " << node_index << " is labelled as an empty site and has a cell attached";
            EXCEPTION(ss.str());
        }
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to be an empty site or have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractCellPopulation<DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::CloseOutputFiles()
{
    AbstractCellPopulation<DIM>::CloseOutputFiles();
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::WriteResultsToFiles()
{
    double time = SimulationTime::Instance()->GetTime();

    *(this->mpVizNodesFile) << time << "\t";
    *(this->mpVizBoundaryNodesFile) << time << "\t";

    // Write node data to file
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        ///\todo we may need a hack to cover the case where the node is associated with a cell that has just been killed (#1129)

        Node<DIM>* p_node = GetNode(node_index);

        // Write node data to file
        if (!(p_node->IsDeleted()))
        {
            const c_vector<double, DIM>& position = p_node->rGetLocation();

            for (unsigned i=0; i<DIM; i++)
            {
                *(this->mpVizNodesFile) << position[i] << " ";
            }
            *(this->mpVizBoundaryNodesFile) << p_node->IsBoundaryNode() << " ";
        }
    }
    *(this->mpVizNodesFile) << "\n";
    *(this->mpVizBoundaryNodesFile) << "\n";

    *(this->mpVizCellProliferativeTypesFile) << time << "\t";

    if (this->mOutputCellAncestors)
    {
        *(this->mpVizCellAncestorsFile) << time << "\t";
    }
    if (this->mOutputCellMutationStates)
    {
        *(this->mpCellMutationStatesFile) << time << "\t";
    }
    if (this->mOutputCellProliferativeTypes)
    {
        *(this->mpCellProliferativeTypesFile) << time << "\t";
    }
    if (this->mOutputCellVariables)
    {
        *(this->mpCellVariablesFile) << time << "\t";
    }
    if (this->mOutputCellCyclePhases)
    {
        *(this->mpCellCyclePhasesFile) << time << "\t";
    }
    if (this->mOutputCellAges)
    {
        *(this->mpCellAgesFile) << time << "\t";
    }

    GenerateCellResultsAndWriteToFiles();

    // Write logged cell data if required
    if (this->mOutputCellIdData)
    {
        this->WriteCellIdDataToFile();
    }
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::WriteCellVolumeResultsToFile()
{   
    // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    // Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); 
         cell_iter!=this->End(); ++cell_iter)
    {
        // Get the index of the corresponding node in mrMesh
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Cell volumes all one (equal-sized lattice sites)
        double cell_volume = 1.0;
           
        // Write node index to file
        *(this->mpCellVolumesFile) << node_index << " ";

        // Write cell ID to file
        unsigned cell_index = cell_iter->GetCellId();
        *(this->mpCellVolumesFile) << cell_index << " ";

        // Write node location to file
        c_vector<double, DIM> node_location = this->GetNode(node_index)->rGetLocation();
        for (unsigned i=0; i<DIM; i++)
        {
            *(this->mpCellVolumesFile) << node_location[i] << " ";
        }

        // Write cell volume (in 3D) or area (in 2D) to file
        *(this->mpCellVolumesFile) << cell_volume << " ";
    }
    *(this->mpCellVolumesFile) << "\n";
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::GenerateCellResults(unsigned locationIndex,
                                                  std::vector<unsigned>& rCellProliferativeTypeCounter,
                                                  std::vector<unsigned>& rCellCyclePhaseCounter)
{
    if (IsEmptySite(locationIndex) == true)
    {
        *(this->mpVizCellProliferativeTypesFile) << INVISIBLE_COLOUR << " ";
    }
    else
    {
        AbstractCellPopulation<DIM>::GenerateCellResults(locationIndex,
                                                 rCellProliferativeTypeCounter,
                                                 rCellCyclePhaseCounter);
    }
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::GenerateCellResultsAndWriteToFiles()
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

    // Write cell data to file
    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        ///\todo we may need a hack that covers the case where the node is associated with a cell that has just been killed (#1129)

        // Write cell data to file
        if (!(this->GetNode(node_index)->IsDeleted()))
        {
            this->GenerateCellResults(node_index, cell_type_counter, cell_cycle_phase_counter);
        }
    }

    this->WriteCellResultsToFiles(cell_type_counter, cell_cycle_phase_counter);
}

template<unsigned DIM>
bool LatticeBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return false;
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::MoveCell(CellPtr pCell, unsigned newLocationIndex)
{
    // Get the current location index corresponding to this cell
    unsigned current_location_index = this->mCellLocationMap[pCell.get()];

    if (current_location_index != newLocationIndex)
    {
        // The new location index must not already correspond to a cell
        assert(IsEmptySite(newLocationIndex));

        // Update the new location index to correspond to this cell
        this->mIsEmptySite[newLocationIndex] = false;
        this->mLocationCellMap[newLocationIndex] = pCell;

        // Update this cell to correspond to the new location index
        this->mCellLocationMap[pCell.get()] = newLocationIndex;

        // If this cell is not a new cell...
        if (current_location_index != UINT_MAX)
        {
            // ...update the current location index to correspond to an empty site
            this->mLocationCellMap.erase(current_location_index);
            this->mIsEmptySite[current_location_index] = true;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//                             Methods Not used                             //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    EXCEPTION("UpdateNodeLocations() cannot be called on a LatticeBasedCellPopulation");
}

template<unsigned DIM>
c_vector<double, DIM> LatticeBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    unsigned node_index = this->mCellLocationMap[pCell.get()];
    c_vector<double, DIM> node_location = this->GetNode(node_index)->rGetLocation();
    return node_location;
}

template<unsigned DIM>
unsigned LatticeBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    EXCEPTION("AddNode() cannot be called on a LatticeBasedCellPopulation");
    return 0;
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    EXCEPTION("SetNode() cannot be called on a LatticeBasedCellPopulation");
}

template<unsigned DIM>
double LatticeBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    EXCEPTION("GetDampingConstant() cannot be called on a LatticeBasedCellPopulation");
    return 0.0;
}

template<unsigned DIM>
void LatticeBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Currently this is not called from LatticeBasedCellBasedSimulation; see #1453 for a discussion on this for centre- and vertex-based cell populations.
    EXCEPTION("OutputCellPopulationParameters() is not yet implemented for LatticeBasedCellPopulation see #1453");
}

template<unsigned DIM>
double LatticeBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = mrMesh.GetWidth(rDimension);

    return width;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class LatticeBasedCellPopulation<1>;
template class LatticeBasedCellPopulation<2>;
template class LatticeBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LatticeBasedCellPopulation)
