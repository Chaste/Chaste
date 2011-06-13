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

#include "PottsBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

void PottsBasedCellPopulation::Validate()
{
    // Check each element has only one cell associated with it
    std::vector<unsigned> validated_element = std::vector<unsigned>(this->GetNumElements(), 0);

    for (AbstractCellPopulation<2>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = GetLocationIndexUsingCell(*cell_iter);
        validated_element[elem_index]++;
    }

    for (unsigned i=0; i<validated_element.size(); i++)
    {
        if (validated_element[i] == 0)
        {
            std::stringstream ss;
            ss << "Element " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str());
        }

        if (validated_element[i] > 1)
        {
            std::stringstream ss;
            ss << "Element " << i << " appears to have " << validated_element[i] << " cells associated with it";
            EXCEPTION(ss.str());
        }
    }
}

PottsBasedCellPopulation::PottsBasedCellPopulation(PottsMesh& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractCellPopulation<2>(rCells, locationIndices),
      mrMesh(rMesh),
      mDeleteMesh(deleteMesh)
{
    // Check each element has only one cell associated with it
    if (validate)
    {
        Validate();
    }
}

PottsBasedCellPopulation::~PottsBasedCellPopulation()
{
    if (mDeleteMesh)
    {
        // Not used until archiving is implemented
        #define COVERAGE_IGNORE
        delete &mrMesh;
        #undef COVERAGE_IGNORE
    }
}

PottsMesh& PottsBasedCellPopulation::rGetMesh()
{
    return mrMesh;
}

const PottsMesh& PottsBasedCellPopulation::rGetMesh() const
{
    return mrMesh;
}

PottsElement* PottsBasedCellPopulation::GetElement(unsigned elementIndex)
{
    return mrMesh.GetElement(elementIndex);
}

unsigned PottsBasedCellPopulation::GetNumElements()
{
    return mrMesh.GetNumElements();
}

Node<2>* PottsBasedCellPopulation::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}

unsigned PottsBasedCellPopulation::GetNumNodes()
{
    return mrMesh.GetNumNodes();
}

c_vector<double, 2> PottsBasedCellPopulation::GetLocationOfCellCentre(CellPtr pCell)
{
    return mrMesh.GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

PottsElement* PottsBasedCellPopulation::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mrMesh.GetElement(this->mCellLocationMap[pCell.get()]);
}

CellPtr PottsBasedCellPopulation::AddCell(CellPtr pNewCell, const c_vector<double,2>& rCellDivisionVector, CellPtr pParentCell)
{
    // Get the element associated with this cell
    PottsElement* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element
    unsigned new_element_index = mrMesh.DivideElement(p_element, true); // new element will be below the existing element

    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->mLocationCellMap[new_element_index] = p_created_cell;
    this->mCellLocationMap[p_created_cell.get()] = new_element_index;
    return p_created_cell;
}

unsigned PottsBasedCellPopulation::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if ((*it)->IsDead())
        {
            // Remove the element from the mesh
            num_removed++;
            mrMesh.DeleteElement(this->mCellLocationMap[(*it).get()]);
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

void PottsBasedCellPopulation::UpdateNodeLocations(const std::vector< c_vector<double, 2> >& rNodeForces, double dt)
{
    /*
     * This is where we currently perform the Monte Carlo simulations.
     * It will eventually be moved to the Simulation class.
     */

    // Loop over nodes and exchange
    ///\todo make this a random sweep (see #1664 and associated tickets)
    for (AbstractMesh<2,2>::NodeIterator node_iter = mrMesh.GetNodeIteratorBegin();
         node_iter != mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {

        // Helper bool to save computational time if there's no point swapping elements
        bool swap_elements=true;

        // All nodes should be in at most one element.
        assert(node_iter->GetNumContainingElements() <= 1);

        // Find a random site from all of the available neighbouring nodes to extend the element/medium into
        std::set<unsigned> neighboring_node_indices = mrMesh.GetNeighbouringNodeIndices(node_iter->GetIndex());

        unsigned new_location_index;

        if (!neighboring_node_indices.empty())
        {
            unsigned num_neighbours = neighboring_node_indices.size();
            unsigned chosen_neighbour = RandomNumberGenerator::Instance()->randMod(num_neighbours);

            std::set<unsigned>::iterator neighbour_iter = neighboring_node_indices.begin();
            for (unsigned i=0; i<chosen_neighbour; i++)
            {
                neighbour_iter++;
            }

            new_location_index = *neighbour_iter;
        }
        else  //(neighboring_node_indices.empty())
        {
            // Every node in your mesh should have at least one neighbour.
            NEVER_REACHED;
        }

        std::set<unsigned> containing_elements = node_iter->rGetContainingElementIndices();
        std::set<unsigned> new_location_containing_elements = GetNode(new_location_index)->rGetContainingElementIndices();

        // All nodes should be in at most one element.
        assert(new_location_containing_elements.size() <= 1);

        double H_0 = 0.0; // Hamiltonian before swap
        double H_1 = 0.0; // Hamiltonian after swap

        // Hamiltonian parameters
        double lambda_volume = 0.1;
        double target_volume = 16.0;
        double lambda_contact = 0.1;

        if (containing_elements.size() == 1) // current node is in an element
        {
            if (new_location_containing_elements.size() == 1) // target node is in an element
            {
                /*
                 * Here the two nodes are in both contained in elements.
                 */

                // Check if the elements are different otherwise no point doing anything
                if ((*containing_elements.begin()) != (*new_location_containing_elements.begin()))
                {
                    /*
                     * Here the two nodes are in different elements, so we should calculate
                     * the Hamiltonian to see whether or not we make the replacement
                     */
                    unsigned current_element = (*containing_elements.begin());
                    unsigned target_element = (*new_location_containing_elements.begin());

                    // This is the Hamiltonian

                    // Add the volume constraint
                    H_0 = lambda_volume*pow(mrMesh.GetVolumeOfElement(current_element)-target_volume, 2.0)+
                          lambda_volume*pow(mrMesh.GetVolumeOfElement(target_element)-target_volume, 2.0);
                    H_1 = lambda_volume*pow(mrMesh.GetVolumeOfElement(current_element)+1.0-target_volume, 2.0)+
                          lambda_volume*pow(mrMesh.GetVolumeOfElement(target_element)-1.0-target_volume, 2.0);

                    // Iterate over nodes neighbouring the target node to work out the contact energy contribution
                    std::set<unsigned> target_neighboring_node_indices = mrMesh.GetNeighbouringNodeIndices(new_location_index);

                    for (std::set<unsigned>::iterator iter = target_neighboring_node_indices.begin();
                         iter != target_neighboring_node_indices.end();
                         ++iter)
                    {
                        std::set<unsigned> neighboring_node_containing_elements = mrMesh.GetNode(*iter)->rGetContainingElementIndices();
                        unsigned neighbour_element = (*neighboring_node_containing_elements.begin());

                        if ( neighboring_node_containing_elements.size() != 0u )
                        {
                            // If the nodes are currently from different elements
                            if ( target_element != neighbour_element )
                            {
                                H_0 += lambda_contact;
                            }

                            // If the nodes will be in different elements after swap
                            if ( current_element != neighbour_element )
                            {
                                 H_1 += lambda_contact;
                            }
                        }
                    }
                }
                else //elements the same so no swap needed
                {
                    // Don't need to do any swapping for this case to save computations.
                    swap_elements = false;
                }
            }
            else // (new_location_containing_elements.size()==0) // target node is not in an element
            {
                assert(new_location_containing_elements.size()==0);

                /*
                 * Here the target node not in an element but the current node is. So we should calculate
                 * the Hamiltonian to see whether or not we make the replacement
                 */

                unsigned current_element = (*containing_elements.begin());

                // This is the Hamiltonian

                // Add the volume constraint
                H_0 = lambda_volume*pow(mrMesh.GetVolumeOfElement(current_element)-target_volume, 2.0);
                H_1 = lambda_volume*pow(mrMesh.GetVolumeOfElement(current_element)+1.0-target_volume, 2.0);

                // Iterate over nodes neighbouring the target node to work out the contact energy contribution
                std::set<unsigned> target_neighboring_node_indices = mrMesh.GetNeighbouringNodeIndices(new_location_index);

                for (std::set<unsigned>::iterator iter = target_neighboring_node_indices.begin();
                     iter != target_neighboring_node_indices.end();
                     ++iter)
                {
                    std::set<unsigned> neighboring_node_containing_elements = mrMesh.GetNode(*iter)->rGetContainingElementIndices();
                    unsigned neighbour_element = (*neighboring_node_containing_elements.begin());

                    // If the neighboring node is currently not in elements
                    if ( neighboring_node_containing_elements.size() != 0u )
                    {
                        H_0 += lambda_contact;
                    }

                    // If the nodes will be in different elements after swap
                    if ( ( current_element != neighbour_element ) && ( neighboring_node_containing_elements.size() != 0u ) )
                    {
                         H_1 += lambda_contact;
                    }
                }
            }
        }
        else // (containing_elements.size() == 0) // current node is not in an element
        {
            if (new_location_containing_elements.size()==1) // target node is in an element
            {
                /*
                 * Here the target node is in an element and the current node is not. So we should calculate
                 * the Hamiltonian to see whether or not we make the replacement
                 */

                unsigned target_element = (*new_location_containing_elements.begin());

                // This is the Hamiltonian

                // Add the volume constraint
                H_0 = lambda_volume*pow(mrMesh.GetVolumeOfElement(target_element)-target_volume, 2.0);
                H_1 = lambda_volume*pow(mrMesh.GetVolumeOfElement(target_element)-1.0-target_volume, 2.0);

                // Iterate over nodes neighbouring the target node to work out the contact energy contribution
                std::set<unsigned> target_neighboring_node_indices = mrMesh.GetNeighbouringNodeIndices(new_location_index);

                for (std::set<unsigned>::iterator iter = target_neighboring_node_indices.begin();
                     iter != target_neighboring_node_indices.end();
                     ++iter)
                {
                    std::set<unsigned> neighboring_node_containing_elements = mrMesh.GetNode(*iter)->rGetContainingElementIndices();
                    unsigned neighbour_element = (*neighboring_node_containing_elements.begin());

                    // If the nodes are currently from different elements
                    if ( ( target_element != neighbour_element ) && ( neighboring_node_containing_elements.size() != 0u ) )
                    {
                        H_0 += lambda_contact;
                    }

                    // If the neighbouring node is not in an element (like the target node after the swap).
                    if ( neighboring_node_containing_elements.size() != 0u )
                    {
                        H_1 += lambda_contact;
                    }
                }

            }
            else // (new_location_containing_elements.size()==0) // target node is not in an element
            {
                // Don't need to do any swapping for this case to save computations.
                swap_elements = false;
            }
        }

        double delta_H = H_1 - H_0;
        double T = 0.1;

        // Generate a uniform random number to do the random motion.
        double random_number = RandomNumberGenerator::Instance()->ranf();
        double p = exp(-delta_H/T);

        if ( ( delta_H <= 0 || random_number < p) && swap_elements )
        {
            // Do swap

            // Iterate over the elements containing the target node to remove node, this should be at most one element.
            for (std::set<unsigned>::iterator iter = new_location_containing_elements.begin();
                 iter != new_location_containing_elements.end();
                 ++iter)
            {
                GetElement(*iter)->DeleteNode(GetElement(*iter)->GetNodeLocalIndex(new_location_index));

                ///\todo If this causes the element to have no nodes then flag the element and cell to be deleted
            }

            // Now iterate over the elements containing the current node to add the target node, this should be at most one element.
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                GetElement(*iter)->AddNode(mrMesh.GetNode(new_location_index));
            }
        }
    }
}

bool PottsBasedCellPopulation::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();;
}

void PottsBasedCellPopulation::Update(bool hasHadBirthsOrDeaths)
{
    ///\todo implement this method
}

void PottsBasedCellPopulation::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractCellPopulation<2>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizElementsFile = output_file_handler.OpenOutputFile("results.vizelements");
}

void PottsBasedCellPopulation::CloseOutputFiles()
{
    AbstractCellPopulation<2>::CloseOutputFiles();
    mpVizElementsFile->close();
}

void PottsBasedCellPopulation::WriteResultsToFiles()
{
    AbstractCellPopulation<2>::WriteResultsToFiles();

    CreateElementTessellation(); //To be used to output to the visualiser

    SimulationTime* p_time = SimulationTime::Instance();

    // Write element data to file
    *mpVizElementsFile << p_time->GetTime() << "\t";

    // Loop over cells and find associated elements so in the same order as the cells in output files
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Hack that covers the case where the element is associated with a cell that has just been killed (#1129)
        bool elem_corresponds_to_dead_cell = false;

        if (this->mLocationCellMap[elem_index])
        {
            elem_corresponds_to_dead_cell = this->mLocationCellMap[elem_index]->IsDead();
        }

        // Write node data to file
        if ( !(GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
            PottsElement* p_element = mrMesh.GetElement(elem_index);

            unsigned num_nodes_in_element = p_element->GetNumNodes();

            // First write the number of Nodes belonging to this PottsElement
            *mpVizElementsFile << num_nodes_in_element << " ";

            // Then write the global index of each Node in this element
            for (unsigned i=0; i<num_nodes_in_element; i++)
            {
                *mpVizElementsFile << p_element->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
    *mpVizElementsFile << "\n";
}

void PottsBasedCellPopulation::WriteCellVolumeResultsToFile()
{    
    // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    // Loop over cells and find associated elements so in the same order as the cells in output files
     for (AbstractCellPopulation<2>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Hack that covers the case where the element is associated with a cell that has just been killed (#1129)
        bool elem_corresponds_to_dead_cell = false;

        if (this->mLocationCellMap[elem_index])
        {
            elem_corresponds_to_dead_cell = this->mLocationCellMap[elem_index]->IsDead();
        }

        // Write node data to file
        if ( !(GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
           // Write element index to file
            *(this->mpCellVolumesFile) << elem_index << " ";

            // Write cell ID to file
            unsigned cell_index = cell_iter->GetCellId();
            *(this->mpCellVolumesFile) << cell_index << " ";

            // Write centroid location to file            
            c_vector<double, 2> centroid_location = mrMesh.GetCentroidOfElement(elem_index); 
            
            *(this->mpCellVolumesFile) << centroid_location[0] << " ";
            *(this->mpCellVolumesFile) << centroid_location[1] << " ";

            // Write cell volume (in 3D) or area (in 2D) to file
            double cell_volume = mrMesh.GetVolumeOfElement(elem_index);
            *(this->mpCellVolumesFile) << cell_volume << " ";
        }
    }
    *(this->mpCellVolumesFile) << "\n";
}

void PottsBasedCellPopulation::GenerateCellResultsAndWriteToFiles()
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

    for (AbstractCellPopulation<2>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        this->GenerateCellResults(this->GetLocationIndexUsingCell(*cell_iter), cell_type_counter, cell_cycle_phase_counter);
    }

    this->WriteCellResultsToFiles(cell_type_counter, cell_cycle_phase_counter);
}

double PottsBasedCellPopulation::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = mrMesh.GetWidth(rDimension);

    return width;
}

void PottsBasedCellPopulation::CreateElementTessellation()
{
    //delete mpElementTessellation;

    // Calculate Voronoi tesselation of the nodes (only needs to be done once)
    //mpVoronoiTessellation = new VertexMesh<2, 2>(tesselation_nodes, is_mesh_periodic);;

    // Create the ElementTeselation
    //mpElementTessellation = mpVoronoiTesselation

    /*
     * For each PottsElement find the Voronoi cells of the contained nodes that share
     * two common voronoi nodes and merge the voronoi cells. Keep searching until
     * there are no such pairs. This gives the ElementTesselation for the PottsElement
     */

//    for (PottsMesh::PottsElementIterator elem_iter = mrMesh.GetElementIteratorBegin();
//         elem_iter != mrMesh.GetElementIteratorEnd();
//         ++elem_iter)
//    {
//        unsigned element_index = elem_iter->GetIndex();
//        unsigned num_nodes = elem_iter->GetNumNodes();
//
//        bool are_shared_edges = true;
//
//        while (are_shared_edges)
//        {
//            for (unsigned i=0; i<num_nodes; i++)
//            {
//                for (unsigned j=0; j<num_nodes; j++)
//                {
//                    // Look for shared voronoi nodes between cell i and j
//                    // if find shared edge perform voronoi cell merge
//                    //mpElementTessellation.MergeElements(mpElementTessellation.GetElement(elem_iter->GetNodeGlobalIndex(i),
//                    //                                    mpElementTessellation.GetElement(elem_iter->GetNodeGlobalIndex(j));
//                    // and break;
//                }
//            }
//        }
//
//    }

    // Make Vertex Mesh
    //mpElementTessellation = new VertexMesh<2, 2>(tesselation_nodes, is_mesh_periodic);
}

VertexMesh<2, 2>* PottsBasedCellPopulation::GetElementTessellation()
{
    return mpElementTessellation;
}

void PottsBasedCellPopulation::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellPopulation<2>::OutputCellPopulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
///\todo Unused Methods to be refactored out of the AbstractCellPopulation?
/////////////////////////////////////////////////////////////////////////////

unsigned PottsBasedCellPopulation::AddNode(Node<2>* pNewNode)
{
    ///\todo Method not needed for this population type; need to refactor out?
    //return mrMesh.AddNode(pNewNode);
    return 0;
}

void PottsBasedCellPopulation::SetNode(unsigned nodeIndex, ChastePoint<2>& rNewLocation)
{
    ///\todo Method not needed for this population type; need to refactor out?
    //mrMesh.SetNode(nodeIndex, rNewLocation);
}

double PottsBasedCellPopulation::GetDampingConstant(unsigned nodeIndex)
{
    ///\todo Method not needed for this population type; need to refactor out?
    #define COVERAGE_IGNORE
    assert(0);
    #undef COVERAGE_IGNORE
    return 0.0;
}

std::set<unsigned> PottsBasedCellPopulation::GetNeighbouringNodeIndices(unsigned index)
{
    // Get pointer to this node
    Node<2>* p_node = mrMesh.GetNode(index);

    // Loop over containing elements
    std::set<unsigned> neighbouring_node_indices;
    for (Node<2>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
         elem_iter != p_node->ContainingElementsEnd();
         ++elem_iter)
    {
        // Get pointer to this containing element
        PottsElement* p_element = mrMesh.GetElement(*elem_iter);

        // Loop over nodes contained in this element
        for (unsigned i=0; i<p_element->GetNumNodes(); i++)
        {
            // Get index of this node and add its index to the set if not the original node
            unsigned node_index = p_element->GetNodeGlobalIndex(i);
            if (node_index != index)
            {
                neighbouring_node_indices.insert(node_index);
            }
        }
    }
    return neighbouring_node_indices;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PottsBasedCellPopulation)
