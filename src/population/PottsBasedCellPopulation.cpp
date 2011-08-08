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
            EXCEPTION("Element " << i << " does not appear to have a cell associated with it");
        }

        if (validated_element[i] > 1)
        {
            EXCEPTION("Element " << i << " appears to have " << validated_element[i] << " cells associated with it");
        }
    }
}

PottsBasedCellPopulation::PottsBasedCellPopulation(PottsMesh<2>& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractCellPopulation<2>(rCells, locationIndices),
      mrMesh(rMesh),
      mDeleteMesh(deleteMesh),
      mTemperature(0.1)
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

PottsMesh<2>& PottsBasedCellPopulation::rGetMesh()
{
    return mrMesh;
}

const PottsMesh<2>& PottsBasedCellPopulation::rGetMesh() const
{
    return mrMesh;
}

PottsElement<2>* PottsBasedCellPopulation::GetElement(unsigned elementIndex)
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

PottsElement<2>* PottsBasedCellPopulation::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mrMesh.GetElement(this->mCellLocationMap[pCell.get()]);
}

CellPtr PottsBasedCellPopulation::AddCell(CellPtr pNewCell, const c_vector<double,2>& rCellDivisionVector, CellPtr pParentCell)
{
    // Get the element associated with this cell
    PottsElement<2>* p_element = GetElementCorrespondingToCell(pParentCell);

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
     * This is where we perform the Monte Carlo simulations.
     */

    // Loop over nodes and exchange
    ///\todo make this a random sweep (see #1664 and associated tickets)
    for (AbstractMesh<2,2>::NodeIterator node_iter = mrMesh.GetNodeIteratorBegin();
         node_iter != mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        // Each node in the mesh must be in at most one element
        assert(node_iter->GetNumContainingElements() <= 1);

        // Find a random available neighbouring node to extend the element/medium into
        std::set<unsigned> neighbouring_node_indices = mrMesh.GetMooreNeighbouringNodeIndices(node_iter->GetIndex());
        unsigned new_location_index;
        if (!neighbouring_node_indices.empty())
        {
            unsigned num_neighbours = neighbouring_node_indices.size();
            unsigned chosen_neighbour = RandomNumberGenerator::Instance()->randMod(num_neighbours);

            std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
            for (unsigned i=0; i<chosen_neighbour; i++)
            {
                neighbour_iter++;
            }

            new_location_index = *neighbour_iter;
        }
        else
        {
            // Each node in the mesh must have at least one neighbour
            NEVER_REACHED;
        }

        std::set<unsigned> containing_elements = node_iter->rGetContainingElementIndices();
        std::set<unsigned> new_location_containing_elements = GetNode(new_location_index)->rGetContainingElementIndices();

        // Only calculate Hamiltonian and update elements if the nodes are from different elements
        if (   ( *containing_elements.begin() != *new_location_containing_elements.begin() )
            && ( !containing_elements.empty() || !new_location_containing_elements.empty() ) )
        {
        	double delta_H = 0.0;

            // Now add contributions to the Hamiltonian from each AbstractPottsUpdateRule
            for (std::vector<AbstractPottsUpdateRule<2>*>::iterator iter = mUpdateRuleCollection.begin();
                 iter != mUpdateRuleCollection.end();
                 ++iter)
            {
                delta_H += (*iter)->EvaluateHamiltonianContribution(node_iter->GetIndex(), new_location_index, *this);
            }

			// Generate a uniform random number to do the random motion
			double random_number = RandomNumberGenerator::Instance()->ranf();
			double p = exp(-delta_H/mTemperature);

			if (delta_H <= 0 || random_number < p)
			{
				// Do swap

			    // Remove the target node from any elements containing it (there should be at most one such element)
				for (std::set<unsigned>::iterator iter = new_location_containing_elements.begin();
					 iter != new_location_containing_elements.end();
					 ++iter)
				{
					GetElement(*iter)->DeleteNode(GetElement(*iter)->GetNodeLocalIndex(new_location_index));

					///\todo If this causes the element to have no nodes then flag the element and cell to be deleted
				}

                // Next add the target node to any elements containing the current node (there should be at most one such element)
				for (std::set<unsigned>::iterator iter = containing_elements.begin();
					 iter != containing_elements.end();
					 ++iter)
				{
					GetElement(*iter)->AddNode(mrMesh.GetNode(new_location_index));
				}
			}
        }
    }
}

bool PottsBasedCellPopulation::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();
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
        if (!(GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
            PottsElement<2>* p_element = mrMesh.GetElement(elem_index);

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
        if (!(GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
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

void PottsBasedCellPopulation::AddUpdateRule(AbstractPottsUpdateRule<2>* pUpdateRule)
{
	mUpdateRuleCollection.push_back(pUpdateRule);
}

void PottsBasedCellPopulation::CreateElementTessellation()
{
	///\todo create a Potts tessellation here to enable VTK output
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
        PottsElement<2>* p_element = mrMesh.GetElement(*elem_iter);

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

void PottsBasedCellPopulation::SetTemperature(double temperature)
{
    mTemperature = temperature;
}

double PottsBasedCellPopulation::GetTemperature()
{
    return mTemperature;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PottsBasedCellPopulation)
