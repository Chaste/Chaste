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

//template<unsigned DIM>
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

//template<unsigned DIM>
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

//template<unsigned DIM>
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

//template<unsigned DIM>
PottsMesh& PottsBasedCellPopulation::rGetMesh()
{
    return mrMesh;
}

//template<unsigned DIM>
const PottsMesh& PottsBasedCellPopulation::rGetMesh() const
{
    return mrMesh;
}

//template<unsigned DIM>
PottsElement* PottsBasedCellPopulation::GetElement(unsigned elementIndex)
{
    return mrMesh.GetElement(elementIndex);
}

//template<unsigned DIM>
unsigned PottsBasedCellPopulation::GetNumElements()
{
    return mrMesh.GetNumElements();
}

//template<unsigned DIM>
Node<2>* PottsBasedCellPopulation::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}

//template<unsigned DIM>
unsigned PottsBasedCellPopulation::GetNumNodes()
{
    return mrMesh.GetNumNodes();
}

//template<unsigned DIM>
c_vector<double, 2> PottsBasedCellPopulation::GetLocationOfCellCentre(CellPtr pCell)
{
    return mrMesh.GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

//template<unsigned DIM>
PottsElement* PottsBasedCellPopulation::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mrMesh.GetElement(this->mCellLocationMap[pCell.get()]);
}

//template<unsigned DIM>
CellPtr PottsBasedCellPopulation::AddCell(CellPtr pNewCell, const c_vector<double,2>& rCellDivisionVector, CellPtr pParentCell)
{
    // Method Not Written Yet
    #define COVERAGE_IGNORE
    assert(0);
    #undef COVERAGE_IGNORE

    // Get the element associated with this cell
    //PottsElement* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
//    this->mLocationCellMap[new_element_index] = p_created_cell;
//    this->mCellLocationMap[p_created_cell.get()] = new_element_index;
    return p_created_cell;
}

//template<unsigned DIM>
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

//template<unsigned DIM>
void PottsBasedCellPopulation::UpdateNodeLocations(const std::vector< c_vector<double, 2> >& rNodeForces, double dt)
{
    /*
     * This is where we currently perform the Monte Carlo simulations.
     * It will eventually be moved to the Simulation class.
     */

    /*
     * Currently we just randomly move elements by looping over all the nodes
     * finding a random neighbour and assigning this neighbour to be the same
     * element as the original node.
     */

    // Loop over nodes and exchange
    ///\todo make this a random sweep (see #1664 and associated tickets)
    for (AbstractMesh<2,2>::NodeIterator node_iter = mrMesh.GetNodeIteratorBegin();
         node_iter != mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        // Node should be in at most one element
        assert(node_iter->GetNumContainingElements() <= 1);

        if (node_iter->GetNumContainingElements() == 1)
        {
            // Find a random site from all of the available neighbouring nodes to extend the element into
            std::set<unsigned> neighboring_node_indices = mrMesh.GetNeighbouringNodeIndices(node_iter->GetIndex());

            if (!neighboring_node_indices.empty())
            {
                unsigned num_neighbours = neighboring_node_indices.size();
                unsigned chosen_neighbour = RandomNumberGenerator::Instance()->randMod(num_neighbours);

                std::set<unsigned>::iterator neighbour_iter = neighboring_node_indices.begin();
                for (unsigned i=0; i<chosen_neighbour; i++)
                {
                    neighbour_iter++;
                }

                double new_location_index = *neighbour_iter;

                std::set<unsigned> containing_elements = node_iter->rGetContainingElementIndices();
                std::set<unsigned> new_location_containing_elements = GetNode(new_location_index)->rGetContainingElementIndices();

                if (containing_elements.size() == 1)
                {
                    std::set<unsigned>::iterator elem_iter = containing_elements.begin();

                    if (new_location_containing_elements.size()>0)
                    {
                        assert(new_location_containing_elements.size()==1);

                        // Check if the elements are different otherwise no point doing anything

                        if( (*elem_iter) != (*new_location_containing_elements.begin()))
                        {

                            /*
                             * Here the two nodes are in different elements, so we should calculate
                             * the Hamiltonian to see whether or not we make the replacement
                             */

                            unsigned curent_element = (*elem_iter);
                            unsigned target_element = (*new_location_containing_elements.begin());

                            // This is the Hamiltonian

                            // Add the volume constraint
                            double lambda_volume = 0.1;
                            double target_volume = 16.0;
                            double H_0 = lambda_volume*pow(mrMesh.GetVolumeOfElement(curent_element)-target_volume, 2.0)+
                                         lambda_volume*pow(mrMesh.GetVolumeOfElement(target_element)-target_volume, 2.0);
                            double H_1 = lambda_volume*pow(mrMesh.GetVolumeOfElement(curent_element)+1.0-target_volume, 2.0)+
                                         lambda_volume*pow(mrMesh.GetVolumeOfElement(target_element)-1.0-target_volume, 2.0);

                            double lambda_contact = 0.1;

                            // Iterate over nodes neighbouring the target node to work out the contact energy contribution
                            std::set<unsigned> target_neighboring_node_indices = mrMesh.GetNeighbouringNodeIndices(new_location_index);

                            for (std::set<unsigned>::iterator iter = target_neighboring_node_indices.begin();
                                 iter != target_neighboring_node_indices.end();
                                 ++iter)
                            {
                                std::set<unsigned> neighboring_node_containing_elements = mrMesh.GetNode(*iter)->rGetContainingElementIndices();
                                unsigned neighbour_element = (*neighboring_node_containing_elements.begin());

                                // If the nodes are currently from different elements
                                if ( target_element != neighbour_element )
                                {
                                    H_0 += lambda_contact;
                                }

                                // If the nodes will be different elements after swap
                                if ( curent_element != neighbour_element )
                                {
                                     H_1 += lambda_contact;
                                }
                            }

                            double delta_H = H_1 - H_0;
                            double T = 0.1;

                            // Uniform random number
                            double random_number = RandomNumberGenerator::Instance()->ranf();
                            double p = exp(-delta_H/T);

                            if (delta_H <= 0 || random_number < p)
                            {
                                // Do swap

                                // Iterate over the elements containing the target node to remove node should be at most one element
                                for (std::set<unsigned>::iterator iter = new_location_containing_elements.begin();
                                     iter != new_location_containing_elements.end();
                                     ++iter)
                                {
                                    GetElement(*iter)->DeleteNode(GetElement(*iter)->GetNodeLocalIndex(new_location_index));

                                    // If this causes the element to have no nodes then flag the element and cell to be deleted
                                }

                                // Now add node to original element
                                GetElement(*elem_iter)->AddNode(mrMesh.GetNode(new_location_index));
                            }
                        }
                    }
                    else // new_location_containing_elements.size()==0
                    {
                        assert(new_location_containing_elements.size()==0);

                        /*
                         * Here the target node not in an element, so we should calculate
                         * the Hamiltonian to see whether or not we make the replacement
                         */

                        unsigned curent_element = (*elem_iter);

                        // This is the Hamiltonian

                        // Add the volume constraint
                        double lambda_volume = 0.1;
                        double target_volume = 16.0;
                        double H_0 = lambda_volume*pow(mrMesh.GetVolumeOfElement(curent_element)-target_volume, 2.0);
                        double H_1 = lambda_volume*pow(mrMesh.GetVolumeOfElement(curent_element)+1.0-target_volume, 2.0);

                        double lambda_contact = 0.1;

                        // Iterate over nodes neighbouring the target node to work out the contact energy contribution
                        std::set<unsigned> target_neighboring_node_indices = mrMesh.GetNeighbouringNodeIndices(new_location_index);

                        for (std::set<unsigned>::iterator iter = target_neighboring_node_indices.begin();
                             iter != target_neighboring_node_indices.end();
                             ++iter)
                        {
                            std::set<unsigned> neighboring_node_containing_elements = mrMesh.GetNode(*iter)->rGetContainingElementIndices();
                            unsigned neighbour_element = (*neighboring_node_containing_elements.begin());

                            // If the nodes are currently from different elements
                            if ( neighboring_node_containing_elements.size() != 0u )
                            {
                                H_0 += lambda_contact;
                            }

                            // If the nodes will be different elements after swap
                            if ( curent_element != neighbour_element )
                            {
                                 H_1 += lambda_contact;
                            }
                        }

                        double delta_H = H_1 - H_0;
                        double T = 0.1;

                        // Uniform random number
                        double random_number = RandomNumberGenerator::Instance()->ranf();
                        double p = exp(-delta_H/T);

                        if (delta_H <= 0 || random_number < p)
                        {
                            // Do swap

                            // Iterate over the elements containing the target node to remove node should be at most one element
                            for (std::set<unsigned>::iterator iter = new_location_containing_elements.begin();
                                 iter != new_location_containing_elements.end();
                                 ++iter)
                            {
                                GetElement(*iter)->DeleteNode(GetElement(*iter)->GetNodeLocalIndex(new_location_index));

                                // If this causes the element to have no nodes then flag the element and cell to be deleted
                            }

                            // Now add node to original element
                            GetElement(*elem_iter)->AddNode(mrMesh.GetNode(new_location_index));
                        }
                    }
                }
                else  //(neighboring_node_indices.empty())
                {
                    NEVER_REACHED;
                }
            }
        }
    }
}

//template<unsigned DIM>
bool PottsBasedCellPopulation::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();;
}

//template<unsigned DIM>
void PottsBasedCellPopulation::Update(bool hasHadBirthsOrDeaths)
{
//    PottsElementMap element_map(mrMesh.GetNumAllElements());
//
//    mrMesh.ReMesh(element_map);
//
//    if (!element_map.IsIdentityMap())
//    {
//    	// Fix up the mappings between CellPtrs and PottsElements
//        std::map<Cell*, unsigned> old_map = this->mCellLocationMap;
//
//        this->mCellLocationMap.clear();
//        this->mLocationCellMap.clear();
//
//        for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
//             cell_iter != this->mCells.end();
//             ++cell_iter)
//        {
//            // This shouldn't ever happen, as the cell vector only contains living cells
//            unsigned old_elem_index = old_map[(*cell_iter).get()];
//
//            if (element_map.IsDeleted(old_elem_index))
//            {
//            	/*\todo this is a kludge to remove the cell once a T2Swap occurs this is not included in the dead cells counter.
//            	 * This should be included in the RemoveDeadCells method so the death is counted
//            	 */
//            	WARNING("Cell removed due to T2Swap this is not counted in the dead cells counter");
//            	cell_iter = this->mCells.erase(cell_iter);
//            	--cell_iter;
//            }
//            else
//            {
//				unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);
//
//				this->mLocationCellMap[new_elem_index] = *cell_iter;
//				this->mCellLocationMap[(*cell_iter).get()] = new_elem_index;
//            }
//        }
//
//        // Check that each PottsElement has only one CellPtr associated with it in the updated cell population
//        Validate();
//    }
//
//    element_map.ResetToIdentity();
}

//template<unsigned DIM>
void PottsBasedCellPopulation::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractCellPopulation<2>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizElementsFile = output_file_handler.OpenOutputFile("results.vizelements");
}

//template<unsigned DIM>
void PottsBasedCellPopulation::CloseOutputFiles()
{
    AbstractCellPopulation<2>::CloseOutputFiles();
    mpVizElementsFile->close();
}

//template<unsigned DIM>
void PottsBasedCellPopulation::WriteResultsToFiles()
{
    // Only works for 2D at present
    //assert(DIM ==2);

    AbstractCellPopulation<2>::WriteResultsToFiles();

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

//template<unsigned DIM>
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

//template<unsigned DIM>
double PottsBasedCellPopulation::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = mrMesh.GetWidth(rDimension);

    return width;
}

//template<unsigned DIM>
void PottsBasedCellPopulation::OutputCellPopulationParameters(out_stream& rParamsFile)
{
	// Call direct parent class method
	AbstractCellPopulation<2>::OutputCellPopulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
///\todo Unused Methods to be refactored out of the AbstractCellPopulation?
/////////////////////////////////////////////////////////////////////////////

//template<unsigned DIM>
unsigned PottsBasedCellPopulation::AddNode(Node<2>* pNewNode)
{
    ///\todo Method not needed for this population type; need to refactor out?
    //return mrMesh.AddNode(pNewNode);
    return 0;
}

//template<unsigned DIM>
void PottsBasedCellPopulation::SetNode(unsigned nodeIndex, ChastePoint<2>& rNewLocation)
{
    ///\todo Method not needed for this population type; need to refactor out?
    //mrMesh.SetNode(nodeIndex, rNewLocation);
}

//template<unsigned DIM>
double PottsBasedCellPopulation::GetDampingConstant(unsigned nodeIndex)
{
    ///\todo Method not needed for this population type; need to refactor out?
    #define COVERAGE_IGNORE
    assert(0);
    #undef COVERAGE_IGNORE
    return 0.0;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class PottsBasedCellPopulation<1>;
//template class PottsBasedCellPopulation<2>;
//template class PottsBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PottsBasedCellPopulation)
