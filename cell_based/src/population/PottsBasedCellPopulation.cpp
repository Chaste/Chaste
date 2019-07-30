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

#include "PottsBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractPottsUpdateRule.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellPopulationElementWriter.hpp"
#include "CellIdWriter.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::Validate()
{
    // Check each element has only one cell associated with it
    std::vector<unsigned> validated_element = std::vector<unsigned>(this->GetNumElements(), 0);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);
        validated_element[elem_index]++;
    }

    for (unsigned i=0; i<validated_element.size(); i++)
    {
        if (validated_element[i] == 0)
        {
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Element " << i << " does not appear to have a cell associated with it");
        }

        if (validated_element[i] > 1)
        {
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Element " << i << " appears to have " << validated_element[i] << " cells associated with it");
        }
    }
}

template<unsigned DIM>
PottsBasedCellPopulation<DIM>::PottsBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                                        std::vector<CellPtr>& rCells,
                                                        bool deleteMesh,
                                                        bool validate,
                                                        const std::vector<unsigned> locationIndices)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh),
      mpElementTessellation(nullptr),
      mpMutableMesh(nullptr),
      mTemperature(0.1),
      mNumSweepsPerTimestep(1)
{
    mpPottsMesh = static_cast<PottsMesh<DIM>* >(&(this->mrMesh));
    // Check each element has only one cell associated with it
    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
PottsBasedCellPopulation<DIM>::PottsBasedCellPopulation(PottsMesh<DIM>& rMesh)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh),
      mpElementTessellation(nullptr),
      mpMutableMesh(nullptr),
      mTemperature(0.1),
      mNumSweepsPerTimestep(1)
{
    mpPottsMesh = static_cast<PottsMesh<DIM>* >(&(this->mrMesh));
}

template<unsigned DIM>
PottsBasedCellPopulation<DIM>::~PottsBasedCellPopulation()
{
    delete mpElementTessellation;

    delete mpMutableMesh;

    if (this->mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
PottsMesh<DIM>& PottsBasedCellPopulation<DIM>::rGetMesh()
{
    return *mpPottsMesh;
}

template<unsigned DIM>
const PottsMesh<DIM>& PottsBasedCellPopulation<DIM>::rGetMesh() const
{
    return *mpPottsMesh;
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* PottsBasedCellPopulation<DIM>::GetTetrahedralMeshForPdeModifier()
{
    std::vector<Node<DIM>*> temp_nodes;

    // Create nodes at the centre of the cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned index = this->GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, DIM> location = this->GetLocationOfCellCentre(*cell_iter);
        temp_nodes.push_back(new Node<DIM>(index, location));
    }

    return new MutableMesh<DIM, DIM>(temp_nodes);
}

template<unsigned DIM>
PottsElement<DIM>* PottsBasedCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mpPottsMesh->GetElement(elementIndex);
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::GetNumElements()
{
    return mpPottsMesh->GetNumElements();
}

template<unsigned DIM>
Node<DIM>* PottsBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
std::set<unsigned> PottsBasedCellPopulation<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);
    return mpPottsMesh->GetNeighbouringElementIndices(elem_index);
}

template<unsigned DIM>
c_vector<double, DIM> PottsBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mpPottsMesh->GetCentroidOfElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
PottsElement<DIM>* PottsBasedCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mpPottsMesh->GetElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
CellPtr PottsBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    // Get the element associated with this cell
    PottsElement<DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element
    unsigned new_element_index = mpPottsMesh->DivideElement(p_element, true); // new element will be below the existing element

    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(new_element_index,p_created_cell);
    return p_created_cell;
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::RemoveDeadCells()
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
            mpPottsMesh->DeleteElement(location_index);

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
void PottsBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
{
    /*
     * This method implements a Monte Carlo method to update the cell population.
     * We sample randomly from all nodes in the mesh. Once we have selected a target
     * node we randomly select a neighbour. The Hamiltonian is evaluated in the
     * current configuration (H_0) and with the target node added to the same
     * element as the neighbour (H_1). Based on the vale of deltaH = H_1 - H_0,
     * the switch is either made or not.
     *
     * For each time step (i.e. each time this method is called) we sample
     * mrMesh.GetNumNodes() nodes. This is known as a Monte Carlo Step (MCS).
     */

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    unsigned num_nodes = this->mrMesh.GetNumNodes();

    // Randomly permute mUpdateRuleCollection if specified
    if (this->mIterateRandomlyOverUpdateRuleCollection)
    {
        // Randomly permute mUpdateRuleCollection
        p_gen->Shuffle(this->mUpdateRuleCollection);
    }

    for (unsigned i=0; i<num_nodes*mNumSweepsPerTimestep; i++)
    {
        unsigned node_index;

        if (this->mUpdateNodesInRandomOrder)
        {
            node_index = p_gen->randMod(num_nodes);
        }
        else
        {
            // Loop over nodes in index order.
            node_index = i%num_nodes;
        }

        Node<DIM>* p_node = this->mrMesh.GetNode(node_index);

        // Each node in the mesh must be in at most one element
        assert(p_node->GetNumContainingElements() <= 1);

        // Find a random available neighbouring node to overwrite current site
        std::set<unsigned> neighbouring_node_indices = mpPottsMesh->GetMooreNeighbouringNodeIndices(node_index);
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

            std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            std::set<unsigned> neighbour_containing_elements = GetNode(neighbour_location_index)->rGetContainingElementIndices();
            // Only calculate Hamiltonian and update elements if the nodes are from different elements, or one is from the medium
            if ((!containing_elements.empty() && neighbour_containing_elements.empty())
                || (containing_elements.empty() && !neighbour_containing_elements.empty())
                || (!containing_elements.empty() && !neighbour_containing_elements.empty() && *containing_elements.begin() != *neighbour_containing_elements.begin()))
            {
                double delta_H = 0.0; // This is H_1-H_0.

                // Now add contributions to the Hamiltonian from each AbstractPottsUpdateRule
                for (typename std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > >::iterator iter = this->mUpdateRuleCollection.begin();
                     iter != this->mUpdateRuleCollection.end();
                     ++iter)
                {
                    // This static cast is fine, since we assert the update rule must be a Potts update rule in AddUpdateRule()
                    double dH = (boost::static_pointer_cast<AbstractPottsUpdateRule<DIM> >(*iter))->EvaluateHamiltonianContribution(neighbour_location_index, p_node->GetIndex(), *this);
                    delta_H += dH;
                }

                // Generate a uniform random number to do the random motion
                double random_number = p_gen->ranf();

                double p = exp(-delta_H/mTemperature);
                if (delta_H <= 0 || random_number < p)
                {
                    // Do swap

                    // Remove the current node from any elements containing it (there should be at most one such element)
                    for (std::set<unsigned>::iterator iter = containing_elements.begin();
                         iter != containing_elements.end();
                         ++iter)
                    {
                        GetElement(*iter)->DeleteNode(GetElement(*iter)->GetNodeLocalIndex(node_index));

                        ///\todo If this causes the element to have no nodes then flag the element and cell to be deleted
                    }

                    // Next add the current node to any elements containing the neighbouring node (there should be at most one such element)
                    for (std::set<unsigned>::iterator iter = neighbour_containing_elements.begin();
                         iter != neighbour_containing_elements.end();
                         ++iter)
                    {
                        GetElement(*iter)->AddNode(this->mrMesh.GetNode(node_index));
                    }
                }
            }
        }
    }
}

template<unsigned DIM>
bool PottsBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellPopulationElementWriter>())
        {
            this-> template AddPopulationWriter<CellPopulationElementWriter>();
        }
    }
    // Add a CellID writer so that a VTK file will contain IDs for visualisation.  (It will also dump a "loggedcell.dat" file as a side-effect.)
    this-> template AddCellWriter<CellIdWriter>();

    AbstractCellPopulation<DIM>::OpenWritersFiles(rOutputFileHandler);
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::WriteResultsToFiles(const std::string& rDirectory)
{
    CreateElementTessellation(); // To be used to output to the visualizer

    AbstractCellPopulation<DIM>::WriteResultsToFiles(rDirectory);
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter)
{
    pPopulationCountWriter->Visit(this);
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}

template<unsigned DIM>
double PottsBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Get element index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get volume of this element in the Potts mesh
    double cell_volume = mpPottsMesh->GetVolumeOfElement(elem_index);

    return cell_volume;
}

template<unsigned DIM>
double PottsBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = this->mrMesh.GetWidth(rDimension);

    return width;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule)
{
    assert(bool(dynamic_cast<AbstractPottsUpdateRule<DIM>*>(pUpdateRule.get())));
    this->mUpdateRuleCollection.push_back(pUpdateRule);
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::CreateElementTessellation()
{
    ///\todo implement this method (#1666)
//  delete mpElementTessellation;
//
//    ///\todo this code would need to be extended if the domain were required to be periodic
//
//  std::vector<Node<2>*> nodes;
//  for (unsigned node_index=0; node_index<mrMesh.GetNumNodes(); node_index++)
//  {
//      Node<2>* p_temp_node = mrMesh.GetNode(node_index);
//      nodes.push_back(p_temp_node);
//  }
//  MutableMesh<2,2> mesh(nodes);
//    mpElementTessellation = new VertexMesh<2,2>(mesh);
}

template<unsigned DIM>
VertexMesh<DIM,DIM>* PottsBasedCellPopulation<DIM>::GetElementTessellation()
{
//    assert(mpElementTessellation != NULL);
    return mpElementTessellation;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::CreateMutableMesh()
{
    delete mpMutableMesh;

    // Get the nodes of the PottsMesh
    std::vector<Node<DIM>*> nodes;
    for (unsigned node_index=0; node_index<this->mrMesh.GetNumNodes(); node_index++)
    {
      const c_vector<double, DIM>& r_location = this->mrMesh.GetNode(node_index)->rGetLocation();
      nodes.push_back(new Node<DIM>(node_index, r_location));
    }

    mpMutableMesh = new MutableMesh<DIM,DIM>(nodes);
}

template<unsigned DIM>
MutableMesh<DIM,DIM>* PottsBasedCellPopulation<DIM>::GetMutableMesh()
{
    assert(mpMutableMesh);
    return mpMutableMesh;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<Temperature>" << mTemperature << "</Temperature>\n";
    *rParamsFile << "\t\t<NumSweepsPerTimestep>" << mNumSweepsPerTimestep << "</NumSweepsPerTimestep>\n";

    // Call method on direct parent class
    AbstractOnLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::SetTemperature(double temperature)
{
    mTemperature = temperature;
}

template<unsigned DIM>
double PottsBasedCellPopulation<DIM>::GetTemperature()
{
    return mTemperature;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::SetNumSweepsPerTimestep(unsigned numSweepsPerTimestep)
{
    mNumSweepsPerTimestep = numSweepsPerTimestep;
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::GetNumSweepsPerTimestep()
{
    return mNumSweepsPerTimestep;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    // Create mesh writer for VTK output
    VtkMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results_"+time.str(), false);

    // Iterate over any cell writers that are present
    unsigned num_nodes = GetNumNodes();
    for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
         cell_writer_iter != this->mCellWriters.end();
         ++cell_writer_iter)
    {
        // Create vector to store VTK cell data
        std::vector<double> vtk_cell_data(num_nodes);

        // Iterate over nodes in the mesh
        for (typename AbstractMesh<DIM,DIM>::NodeIterator iter = mpPottsMesh->GetNodeIteratorBegin();
             iter != mpPottsMesh->GetNodeIteratorEnd();
             ++iter)
        {
            // Get the index of this node in the mesh and those elements (i.e. cells) that contain this node
            unsigned node_index = iter->GetIndex();
            std::set<unsigned> element_indices = iter->rGetContainingElementIndices();

            // If there are no elements associated with this node, then we set the value of any VTK cell data to be -1 at this node...
            if (element_indices.empty())
            {
                // Populate the vector of VTK cell data
                vtk_cell_data[node_index] = -1.0;
            }
            else
            {
                // ... otherwise there should be exactly one element (i.e. cell) containing this node
                assert(element_indices.size() == 1);
                unsigned elem_index = *(element_indices.begin());
                CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);

                // Populate the vector of VTK cell data
                vtk_cell_data[node_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
            }
        }

        mesh_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    }

    // When outputting any CellData, we assume that the first cell is representative of all cells
    unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

    std::vector<std::vector<double> > cell_data;
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cell_data_var(num_nodes);
        cell_data.push_back(cell_data_var);
    }

    for (typename AbstractMesh<DIM,DIM>::NodeIterator iter = mpPottsMesh->GetNodeIteratorBegin();
         iter != mpPottsMesh->GetNodeIteratorEnd();
         ++iter)
    {
        // Get the index of this node in the mesh and those elements (i.e. cells) that contain this node
        unsigned node_index = iter->GetIndex();
        std::set<unsigned> element_indices = iter->rGetContainingElementIndices();

        // If there are no elements associated with this node, then we set the value of any VTK cell data to be -1 at this node...
        if (element_indices.empty())
        {
            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cell_data[var][node_index] = -1.0;
            }
        }
        else
        {
            // ... otherwise there should be exactly one element (i.e. cell) containing this node
            assert(element_indices.size() == 1);
            unsigned elem_index = *(element_indices.begin());
            CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);

            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cell_data[var][node_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
            }
        }
    }
    for (unsigned var=0; var<cell_data.size(); var++)
    {
        mesh_writer.AddPointData(cell_data_names[var], cell_data[var]);
    }

    /*
     * At present, the VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK, then visualized
     * as glyphs in Paraview.
     */
    NodesOnlyMesh<DIM> temp_mesh;
    temp_mesh.ConstructNodesWithoutMesh(*mpPottsMesh, 1.5); // This cut-off is arbitrary, as node connectivity is not used here
    mesh_writer.WriteFilesUsingMesh(temp_mesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";

    // Extra part to output the outlines of cells
    if (DIM ==2)
    {
        std::vector<Node<2>*> outline_nodes;
        std::vector<VertexElement<2,2>*>  outline_elements;

        unsigned outline_node_index = 0;
        unsigned outline_element_index = 0;
        for (typename AbstractMesh<DIM,DIM>::NodeIterator iter = mpPottsMesh->GetNodeIteratorBegin();
             iter != mpPottsMesh->GetNodeIteratorEnd();
             ++iter)
        {
            // Get the index of this node in the mesh and those elements (i.e. cells) that contain this node
            unsigned node_index = iter->GetIndex();
            std::set<unsigned> element_indices = iter->rGetContainingElementIndices();

            std::set<unsigned> target_neighbouring_node_indices = this->rGetMesh().GetVonNeumannNeighbouringNodeIndices(node_index);

            for (std::set<unsigned>::iterator neighbour_iter = target_neighbouring_node_indices.begin();
                 neighbour_iter != target_neighbouring_node_indices.end();
                 ++neighbour_iter)
            {
                std::set<unsigned> neighbouring_element_indices = this->rGetMesh().GetNode(*neighbour_iter)->rGetContainingElementIndices();

                // If different cells add a line
                if (element_indices != neighbouring_element_indices)
                {
                    std::vector<Node<2>*> element_nodes;

                    const c_vector<double, 2>& r_node_location = this->mrMesh.GetNode(node_index)->rGetLocation();
                    const c_vector<double, 2>& r_neighbour_node_location = this->mrMesh.GetNode(*neighbour_iter)->rGetLocation();

                    c_vector<double, 2> unit_tangent = r_neighbour_node_location - r_node_location;

                    if (norm_2(unit_tangent) > 1.0) // It's a periodic neighbour
                    {
                        c_vector<double, 2> mid_point = 0.5*(r_node_location + r_neighbour_node_location);
                        if (unit_tangent(0) == 0)
                        {
                            if (r_node_location(1) < r_neighbour_node_location(1))
                            {
                                mid_point(1) = r_node_location(1) - 0.5;
                            }
                            else
                            {
                                mid_point(1) = r_node_location(1) + 0.5;
                            }
                        }
                        else
                        {
                            assert(unit_tangent(1) == 0);

                            if (r_node_location(0) < r_neighbour_node_location(0))
                            {
                                mid_point(0) = r_node_location(0) - 0.5;
                            }
                            else
                            {
                                mid_point(0) = r_node_location(0) + 0.5;
                            }
                        }
                        assert(norm_2(unit_tangent) > 0);
                        unit_tangent = unit_tangent/norm_2(unit_tangent);

                        c_vector<double, DIM> unit_normal;
                        unit_normal(0) = -unit_tangent(1);
                        unit_normal(1) = unit_tangent(0);

                        std::vector<Node<2>*> element_nodes;

                        // Need at least three points to visualise in Paraview
                        Node<2>* p_node_1 = new Node<2>(outline_node_index, mid_point - 0.5*unit_normal);
                        outline_nodes.push_back(p_node_1);
                        element_nodes.push_back(outline_nodes[outline_node_index]);
                        outline_node_index++;

                        Node<2>* p_node_2 = new Node<2>(outline_node_index, mid_point);
                        outline_nodes.push_back(p_node_2);
                        element_nodes.push_back(outline_nodes[outline_node_index]);
                        outline_node_index++;

                        Node<2>* p_node_3 = new Node<2>(outline_node_index, mid_point + 0.5*unit_normal);
                        outline_nodes.push_back(p_node_3);
                        element_nodes.push_back(outline_nodes[outline_node_index]);
                        outline_node_index++;

                        VertexElement<2,2>* p_element = new VertexElement<2,2>(outline_element_index, element_nodes);
                        outline_elements.push_back(p_element);
                        outline_element_index++;

                    }
                    else // Standard Neighbour
                    {
                        c_vector<double, 2> mid_point = 0.5*(r_node_location + r_neighbour_node_location);

                        assert(norm_2(unit_tangent) > 0);
                        unit_tangent = unit_tangent/norm_2(unit_tangent);

                        c_vector<double, 2> unit_normal;
                        unit_normal(0) = -unit_tangent(1);
                        unit_normal(1) = unit_tangent(0);

                        std::vector<Node<2>*> element_nodes;

                        // Need at least three points to visualise in Paraview
                        Node<2>* p_node_1 = new Node<2>(outline_node_index, mid_point - 0.5*unit_normal);
                        outline_nodes.push_back(p_node_1);
                        element_nodes.push_back(outline_nodes[outline_node_index]);
                        outline_node_index++;

                        Node<2>* p_node_2 = new Node<2>(outline_node_index, mid_point);
                        outline_nodes.push_back(p_node_2);
                        element_nodes.push_back(outline_nodes[outline_node_index]);
                        outline_node_index++;

                        Node<2>* p_node_3 = new Node<2>(outline_node_index, mid_point + 0.5*unit_normal);
                        outline_nodes.push_back(p_node_3);
                        element_nodes.push_back(outline_nodes[outline_node_index]);
                        outline_node_index++;

                        VertexElement<2,2>* p_element = new VertexElement<2,2>(outline_element_index, element_nodes);
                        outline_elements.push_back(p_element);
                        outline_element_index++;
                    }
                }
            }
        }
        VertexMesh<2,2> cell_outline_mesh(outline_nodes,outline_elements);

        VertexMeshWriter<2, 2> outline_mesh_writer(rDirectory, "outlines", false);
        outline_mesh_writer.WriteVtkUsingMesh(cell_outline_mesh, time.str());
        outline_mesh_writer.WriteFilesUsingMesh(cell_outline_mesh);
    }
#endif //CHASTE_VTK
}

template<unsigned DIM>
double PottsBasedCellPopulation<DIM>::GetCellDataItemAtPdeNode(
    unsigned pdeNodeIndex,
    std::string& rVariableName,
    bool dirichletBoundaryConditionApplies,
    double dirichletBoundaryValue)
{
    CellPtr p_cell = this->GetCellUsingLocationIndex(pdeNodeIndex);
    double value = p_cell->GetCellData()->GetItem(rVariableName);
    return value;
}

// Explicit instantiation
template class PottsBasedCellPopulation<1>;
template class PottsBasedCellPopulation<2>;
template class PottsBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsBasedCellPopulation)
