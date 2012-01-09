/*

Copyright (C) University of Oxford, 2005-2012

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

#include <cassert>
#include <algorithm>

#include "CaBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"
#include "NodesOnlyMesh.hpp"
#include "Exception.hpp"

template<unsigned DIM>
CaBasedCellPopulation<DIM>::CaBasedCellPopulation(TetrahedralMesh<DIM, DIM>& rMesh,
                                            std::vector<CellPtr>& rCells,
                                            const std::vector<unsigned> locationIndices,
                                            bool deleteMesh,
                                            bool validate)
    : AbstractOnLatticeCellPopulation<DIM>(rCells, locationIndices),
      mrMesh(rMesh),
      mOnlyUseNearestNeighboursForDivision(false),
      mUseVonNeumannNeighbourhoods(false)
{
    // This must always be true
    assert(this->mCells.size() <= mrMesh.GetNumNodes());

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
        mIsEmptySite = std::vector<bool>(this->GetNumNodes(), false);
        Validate();
    }
}

template<unsigned DIM>
CaBasedCellPopulation<DIM>::CaBasedCellPopulation(TetrahedralMesh<DIM, DIM>& rMesh)
    : AbstractOnLatticeCellPopulation<DIM>(),
      mrMesh(rMesh),
      mOnlyUseNearestNeighboursForDivision(false),
      mUseVonNeumannNeighbourhoods(false)
{
}

template<unsigned DIM>
CaBasedCellPopulation<DIM>::~CaBasedCellPopulation()
{
    if (this->mDeleteMesh)
    {
        delete &mrMesh;
    }
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
{
    // Iterate over contributions from each UpdateRule
    if (this->mIterateRandomlyOverUpdateRuleCollection)
    {
        // Randomly permute mUpdateRuleCollection
        /// \todo #1942 this call will invalidate the state of the Random Number Generator
        std::random_shuffle(mUpdateRuleCollection.begin(), mUpdateRuleCollection.end());
    }

    for (typename std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >::iterator update_iter = mUpdateRuleCollection.begin();
         update_iter != mUpdateRuleCollection.end();
         ++update_iter)
    {
        // Randomly permute cells
        if (this->mUpdateNodesInRandomOrder)
        {
            std::vector<CellPtr> cells_vector;
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
                 cell_iter != this->End();
                 ++cell_iter)
            {
                cells_vector.push_back(*cell_iter);
            }
            /// \todo #1942 this call will invalidate the state of the Random Number Generator
            std::random_shuffle(cells_vector.begin(), cells_vector.end());

            for (unsigned i=0; i<cells_vector.size(); i++)
            {
                // Get index of the node associated with cell
                unsigned current_location_index = this->GetLocationIndexUsingCell(cells_vector[i]);

                assert(!this->IsEmptySite(current_location_index));

                // Get index of node the cell is to move to
                unsigned new_location_index = (*update_iter)->GetNewLocationOfCell(current_location_index, *this, dt);

                // Update the location index of the cell and free the old site
                this->MoveCell(cells_vector[i], new_location_index);
            }
        }
        else
        {
            // Iterate over all cells and update their positions
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
                 cell_iter != this->End();
                 ++cell_iter)
            {
                // Get index of the node associated with cell
                unsigned current_location_index = this->GetLocationIndexUsingCell(*cell_iter);

                assert(!this->IsEmptySite(current_location_index));

                // Get index of node the cell is to move to
                unsigned new_location_index = (*update_iter)->GetNewLocationOfCell(current_location_index, *this, dt);

                // Update the location index of the cell and free the old site
                this->MoveCell(*cell_iter, new_location_index);
            }
        }
    }
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::SetOnlyUseNearestNeighboursForDivision(bool onlyUseNearestNeighboursForDivision)
{
    mOnlyUseNearestNeighboursForDivision = onlyUseNearestNeighboursForDivision;
}

template<unsigned DIM>
bool CaBasedCellPopulation<DIM>::GetOnlyUseNearestNeighboursForDivision()
{
    return mOnlyUseNearestNeighboursForDivision;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule)
{
    mUpdateRuleCollection.push_back(pUpdateRule);
}

template<unsigned DIM>
const std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >& CaBasedCellPopulation<DIM>::rGetUpdateRuleCollection() const
{
    return mUpdateRuleCollection;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::SetUseVonNeumannNeighbourhoods(bool useVonNeumannNeighbourhoods)
{
    mUseVonNeumannNeighbourhoods = useVonNeumannNeighbourhoods;
}

template<unsigned DIM>
bool CaBasedCellPopulation<DIM>::GetUseVonNeumannNeighbourhoods()
{
    return mUseVonNeumannNeighbourhoods;
}

template<unsigned DIM>
std::vector<bool>& CaBasedCellPopulation<DIM>::rGetEmptySites()
{
    return mIsEmptySite;
}

template<unsigned DIM>
bool CaBasedCellPopulation<DIM>::IsEmptySite(unsigned index)
{
    return mIsEmptySite[index];
}

template<unsigned DIM>
std::set<unsigned> CaBasedCellPopulation<DIM>::GetEmptySiteIndices()
{
    std::set<unsigned> empty_site_indices;
    for (unsigned i=0; i<mIsEmptySite.size(); i++)
    {
        if (mIsEmptySite[i])
        {
            empty_site_indices.insert(i);
        }
    }
    return empty_site_indices;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::SetEmptySites(const std::set<unsigned>& rEmptySiteIndices)
{
    // Reinitialise all entries of mIsEmptySite to false
    mIsEmptySite = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsEmptySite
    for (std::set<unsigned>::iterator iter = rEmptySiteIndices.begin();
         iter != rEmptySiteIndices.end();
         ++iter)
    {
        mIsEmptySite[*iter] = true;
    }

    Validate();
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>& CaBasedCellPopulation<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
const TetrahedralMesh<DIM, DIM>& CaBasedCellPopulation<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
Node<DIM>* CaBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned CaBasedCellPopulation<DIM>::GetNumNodes()
{
    return mrMesh.GetNumAllNodes();
}

template<unsigned DIM>
CellPtr CaBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    /*
     * In the case of a CaBasedCellPopulation, we must provided a parent cell when calling AddCell().
     * This is because the location of the parent cell is used as a starting point in the search
     * among neighbours for an empty site in which to locate the new cell. Therefore, if no parent
     * cell is provided, we throw the following exception.
     */
    if (pParentCell == CellPtr())
    {
        EXCEPTION("A parent cell must be provided when calling AddCell() on a CaBasedCellPopulation.");
    }

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
        degree_upper_bound = 1;
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
        EXCEPTION("Cell can not divide as there are no free neighbours at maximum degree in any direction.");
    }

    return p_created_cell;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::WriteVtkResultsToFile()
{
    ///\todo Implement VTK output for CA simulations (see #1914)
//#ifdef CHASTE_VTK
//    std::stringstream time;
//    time << SimulationTime::Instance()->GetTimeStepsElapsed();
//    VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);
//
//    unsigned num_nodes = GetNumNodes();
//    std::vector<double> cell_types;
//    std::vector<double> cell_mutation_states;
//    std::vector<double> cell_labels;
//    cell_types.reserve(num_nodes);
//    cell_mutation_states.reserve(num_nodes);
//    cell_labels.reserve(num_nodes);
//
//    for (unsigned node_index=0; node_index<num_nodes; node_index++)
//    {
//        if (node_index)
//        {
//            // No cell is associated with this node
//            cell_types.push_back(-1.0);
//            if (this->mOutputCellMutationStates)
//            {
//                cell_mutation_states.push_back(-1.0);
//                cell_labels.push_back(-1.0);
//            }
//        }
//        else
//        {
//            CellPtr p_cell = this->mLocationCellMap[node_index];
//            double cell_type = p_cell->GetCellCycleModel()->GetCellProliferativeType();
//            cell_types.push_back(cell_type);
//
//            if (this->mOutputCellMutationStates)
//            {
//                double cell_mutation_state = p_cell->GetMutationState()->GetColour();
//                cell_mutation_states.push_back(cell_mutation_state);
//
//                double cell_label = 0.0;
//                if (p_cell->HasCellProperty<CellLabel>())
//                {
//                    CellPropertyCollection collection = p_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
//                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
//                    cell_label = p_label->GetColour();
//                }
//                cell_labels.push_back(cell_label);
//            }
//        }
//    }
//
//    assert(cell_types.size() == num_nodes);
//
//    mesh_writer.AddPointData("Cell types", cell_types);
//
//    if (this->mOutputCellMutationStates)
//    {
//        assert(cell_mutation_states.size() == num_nodes);
//        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
//        assert(cell_labels.size() == num_nodes);
//        mesh_writer.AddPointData("Cell labels", cell_labels);
//    }
//
//    /*
//     * The current VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
//     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK then visualized as glyphs.
//     */
//    NodesOnlyMesh<DIM> temp_mesh;
//    temp_mesh.ConstructNodesWithoutMesh(mrMesh);
//    mesh_writer.WriteFilesUsingMesh(temp_mesh);
//
//    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
//    *(this->mpVtkMetaFile) << time.str();
//    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
//    *(this->mpVtkMetaFile) << time.str();
//    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
//#endif //CHASTE_VTK
}

template<unsigned DIM>
std::vector<unsigned> CaBasedCellPopulation<DIM>::GetMaximumDegreeInEachDirection(unsigned nodeIndex)
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
std::set<unsigned> CaBasedCellPopulation<DIM>::GetNthDegreeNeighbouringNodeIndices(unsigned nodeIndex, unsigned degree)
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
std::set<unsigned> CaBasedCellPopulation<DIM>::GetFreeNeighbouringNodeIndices(unsigned nodeIndex)
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
std::set<unsigned> CaBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
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
std::vector<unsigned> CaBasedCellPopulation<DIM>::GetNeighbouringNodeIndicesVector(unsigned nodeIndex)
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
unsigned CaBasedCellPopulation<DIM>::RemoveDeadCells()
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
            mIsEmptySite[node_index] = true;
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
void CaBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    Validate();
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::Validate()
{
    // Get a list of all the sites that are empty
    std::vector<bool> validated_node = mIsEmptySite;

    assert(mIsEmptySite.size() == this->GetNumNodes());

    // Look through all of the cells and record what node they are associated with
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];

        // If the node attached to this cell is labelled as an empty site, then throw an error
        if (mIsEmptySite[node_index])
        {
            EXCEPTION("Node " << node_index << " is labelled as an empty site and has a cell attached");
        }
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            EXCEPTION("Node " << i << " does not appear to be an empty site or have a cell associated with it");
        }
    }
}

template<unsigned DIM>
double CaBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{    
    // Cell volumes all one (equal-sized lattice sites)
    ///\todo modify this method to account for more general lattices
    double cell_volume = 1.0;

    return cell_volume;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::WriteCellVolumeResultsToFile()
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
void CaBasedCellPopulation<DIM>::GenerateCellResults(unsigned locationIndex,
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
void CaBasedCellPopulation<DIM>::GenerateCellResultsAndWriteToFiles()
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
bool CaBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return false;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::MoveCell(CellPtr pCell, unsigned newLocationIndex)
{
    // Get the current location index corresponding to this cell
    unsigned current_location_index = this->mCellLocationMap[pCell.get()];

    if (current_location_index != newLocationIndex)
    {
        // The new location index must not already correspond to a cell
        assert(IsEmptySite(newLocationIndex));

        // Update the new location index to correspond to this cell
        mIsEmptySite[newLocationIndex] = false;
        this->mLocationCellMap[newLocationIndex] = pCell;

        // Update this cell to correspond to the new location index
        this->mCellLocationMap[pCell.get()] = newLocationIndex;

        // If this cell is not a new cell...
        if (current_location_index != UINT_MAX)
        {
            // ...update the current location index to correspond to an empty site
            this->mLocationCellMap.erase(current_location_index);
            mIsEmptySite[current_location_index] = true;
        }
    }
}

template<unsigned DIM>
c_vector<double, DIM> CaBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    unsigned node_index = this->mCellLocationMap[pCell.get()];
    c_vector<double, DIM> node_location = this->GetNode(node_index)->rGetLocation();
    return node_location;
}

template<unsigned DIM>
void CaBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OnlyUseNearestNeighboursForDivision>" << mOnlyUseNearestNeighboursForDivision << "</OnlyUseNearestNeighboursForDivision>\n";
    *rParamsFile << "\t\t<UseVonNeumannNeighbourhoods>" << mUseVonNeumannNeighbourhoods << "</UseVonNeumannNeighbourhoods>\n";

    // Call method on direct parent class
    AbstractCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
double CaBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = mrMesh.GetWidth(rDimension);

    return width;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CaBasedCellPopulation<1>;
template class CaBasedCellPopulation<2>;
template class CaBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedCellPopulation)
