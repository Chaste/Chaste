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

#include "VertexBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "VertexMeshWriter.hpp"
#include "Warnings.hpp"

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<DIM>(rCells, locationIndices),
      mrMesh(rMesh),
      mDeleteMesh(deleteMesh)
{
    // Check the mesh contains boundary nodes
    bool contains_boundary_nodes = false;
    for (typename MutableVertexMesh<DIM,DIM>::NodeIterator node_iter = mrMesh.GetNodeIteratorBegin();
         node_iter != mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        if (node_iter->IsBoundaryNode())
        {
            contains_boundary_nodes = true;
        }
    }
    //if (mrMesh.GetNumBoundaryNodes() == 0)
    ///\todo we should be able to do this, but mBoundaryNodes is not used in vertex meshes (#1558)
    if (!contains_boundary_nodes)
    {
        EXCEPTION("No boundary nodes are defined in the supplied vertex mesh which are needed for vertex based simulations.");
    }

    // Check each element has only one cell attached
    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh)
    : AbstractOffLatticeCellPopulation<DIM>(),
      mrMesh(rMesh),
      mDeleteMesh(true)
{
}

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::~VertexBasedCellPopulation()
{
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    // Take the average of the cells containing this vertex
    double average_damping_constant = 0.0;

    std::set<unsigned> containing_elements = GetNode(nodeIndex)->rGetContainingElementIndices();
    double temp = 1.0/((double) containing_elements.size());

    for (std::set<unsigned>::iterator iter = containing_elements.begin();
         iter != containing_elements.end();
         ++iter)
    {
        CellPtr p_cell = this->GetCellUsingLocationIndex(*iter);
        bool cell_is_wild_type = p_cell->GetMutationState()->IsType<WildTypeCellMutationState>();
        bool cell_is_labelled = p_cell->HasCellProperty<CellLabel>();

        if (cell_is_wild_type && !cell_is_labelled)
        {
            average_damping_constant += this->GetDampingConstantNormal()*temp;
        }
        else
        {
            average_damping_constant += this->GetDampingConstantMutant()*temp;
        }
    }

    return average_damping_constant;
}

template<unsigned DIM>
MutableVertexMesh<DIM, DIM>& VertexBasedCellPopulation<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
const MutableVertexMesh<DIM, DIM>& VertexBasedCellPopulation<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mrMesh.GetElement(elementIndex);
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::GetNumNodes()
{
    return mrMesh.GetNumNodes();
}

template<unsigned DIM>
c_vector<double, DIM> VertexBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mrMesh.GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

template<unsigned DIM>
Node<DIM>* VertexBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mrMesh.AddNode(pNewNode);
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mrMesh.SetNode(nodeIndex, rNewLocation);
}

template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mrMesh.GetElement(this->mCellLocationMap[pCell.get()]);
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::GetNumElements()
{
    return mrMesh.GetNumElements();
}

template<unsigned DIM>
CellPtr VertexBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    // Get the element associated with this cell
    VertexElement<DIM, DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element
    unsigned new_element_index;
    if (norm_2(rCellDivisionVector) < DBL_EPSILON)
    {
        // If the cell division vector is the default zero vector, divide the element along the short axis
        new_element_index = mrMesh.DivideElementAlongShortAxis(p_element, true);
    }
    else
    {
        // If the cell division vector has any non-zero component, divide the element along this axis
        new_element_index = mrMesh.DivideElementAlongGivenAxis(p_element, rCellDivisionVector, true);
    }

    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->mLocationCellMap[new_element_index] = p_created_cell;
    this->mCellLocationMap[p_created_cell.get()] = new_element_index;
    return p_created_cell;
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::RemoveDeadCells()
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
            mrMesh.DeleteElementPriorToReMesh(this->mCellLocationMap[(*it).get()]);
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        // Get the damping constant for this node
        double damping_const = this->GetDampingConstant(node_index);

        // Compute the displacement of this node
        c_vector<double, DIM> displacement = dt*rNodeForces[node_index]/damping_const;

        /*
         * If the displacement of this node is greater than half the cell rearrangement threshold,
         * this could result in nodes moving into the interior of other elements, which should not
         * be possible. Therefore in this case we restrict the displacement of the node to the cell
         * rearrangement threshold and warn the user that a smaller timestep should be used. This
         * restriction ensures that vertex elements remain well defined (see #1376).
         */
        if (norm_2(displacement) > 0.5*mrMesh.GetCellRearrangementThreshold())
        {
            WARN_ONCE_ONLY("Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
            displacement *= 0.5*mrMesh.GetCellRearrangementThreshold()/norm_2(displacement);
        }

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + displacement;

        for (unsigned i=0; i<DIM; i++)
        {
            assert(!std::isnan(new_node_location(i)));
        }

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}

template<unsigned DIM>
bool VertexBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    VertexElementMap element_map(mrMesh.GetNumAllElements());

    mrMesh.ReMesh(element_map);

    if (!element_map.IsIdentityMap())
    {
        // Fix up the mappings between CellPtrs and VertexElements
        std::map<Cell*, unsigned> old_map = this->mCellLocationMap;

        this->mCellLocationMap.clear();
        this->mLocationCellMap.clear();

        for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            // This shouldn't ever happen, as the cell vector only contains living cells
            unsigned old_elem_index = old_map[(*cell_iter).get()];

            if (element_map.IsDeleted(old_elem_index))
            {
                /**
                 * \todo this is a kludge to remove the cell once a T2Swap occurs this is not included in the dead cells counter.
                 * This should be included in the RemoveDeadCells method so the death is counted
                 */
                WARNING("Cell removed due to T2Swap this is not counted in the dead cells counter");
                cell_iter = this->mCells.erase(cell_iter);
                --cell_iter;
            }
            else
            {
                unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);

                this->mLocationCellMap[new_elem_index] = *cell_iter;
                this->mCellLocationMap[(*cell_iter).get()] = new_elem_index;
            }
        }

        // Check that each VertexElement has only one CellPtr associated with it in the updated cell population
        Validate();
    }

    element_map.ResetToIdentity();
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::Validate()
{
    // Check each element has only one cell attached
    std::vector<unsigned> validated_element = std::vector<unsigned>(this->GetNumElements(), 0);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
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

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::WriteResultsToFiles()
{
    AbstractOffLatticeCellPopulation<DIM>::WriteResultsToFiles();

    SimulationTime* p_time = SimulationTime::Instance();

    // Write Locations of T1Swaps to file
    *mpT1SwapLocationsFile << p_time->GetTime() << "\t";
    std::vector< c_vector<double, DIM> > t1_swap_locations = mrMesh.GetLocationsOfT1Swaps();
    *mpT1SwapLocationsFile << t1_swap_locations.size() << "\t";
    for (unsigned index = 0;  index < t1_swap_locations.size(); index++)
    {
        for (unsigned i=0; i<DIM; i++)
        {
            *mpT1SwapLocationsFile << t1_swap_locations[index][i] << "\t";
        }
    }
    *mpT1SwapLocationsFile << "\n";
    mrMesh.ClearLocationsOfT1Swaps();

    // Write Locations of T3Swaps to file
    *mpT3SwapLocationsFile << p_time->GetTime() << "\t";
    std::vector< c_vector<double, DIM> > t3_swap_locations = mrMesh.GetLocationsOfT3Swaps();
    *mpT3SwapLocationsFile << t3_swap_locations.size() << "\t";
    for (unsigned index = 0;  index < t3_swap_locations.size(); index++)
    {
        for (unsigned i=0; i<DIM; i++)
        {
            *mpT3SwapLocationsFile << t3_swap_locations[index][i] << "\t";
        }
    }
    *mpT3SwapLocationsFile << "\n";
    mrMesh.ClearLocationsOfT3Swaps();

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

        // Write element data to file
        if (!(GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
            VertexElement<DIM, DIM>* p_element = mrMesh.GetElement(elem_index);

            unsigned num_nodes_in_element = p_element->GetNumNodes();

            // First write the number of Nodes belonging to this VertexElement
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

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::WriteCellVolumeResultsToFile()
{
    assert(DIM==2);

     // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    // Loop over cells and find associated elements so in the same order as the cells in output files
     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
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
            unsigned cell_index = (*cell_iter)->GetCellId();
            *(this->mpCellVolumesFile) << cell_index << " ";

            // Write location of element centroid to file
            c_vector<double, DIM> centre_location = GetLocationOfCellCentre(*cell_iter);
            for (unsigned i=0; i<DIM; i++)
            {
                *(this->mpCellVolumesFile) << centre_location[i] << " ";
            }

            // Write cell volume (in 3D) or area (in 2D) to file
            double cell_volume = mrMesh.GetVolumeOfElement(elem_index);
            *(this->mpCellVolumesFile) << cell_volume << " ";
        }
    }
    *(this->mpCellVolumesFile) << "\n";
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    SimulationTime* p_time = SimulationTime::Instance();

    VertexMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results", false);
    std::stringstream time;
    time << p_time->GetTimeStepsElapsed();

    unsigned num_elements = mrMesh.GetNumElements();
    std::vector<double> cell_types(num_elements);
    std::vector<double> cell_ancestors(num_elements);
    std::vector<double> cell_mutation_states(num_elements);
    std::vector<double> cell_ages(num_elements);
    std::vector<double> cell_cycle_phases(num_elements);
    std::vector<double> cell_volumes(num_elements);
    std::vector<std::vector<double> > cellwise_data;

    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
        unsigned num_variables = p_data->GetNumVariables();
        for (unsigned var=0; var<num_variables; var++)
        {
            std::vector<double> cellwise_data_var(num_elements);
            cellwise_data.push_back(cellwise_data_var);
        }
    }

    // Loop over vertex elements
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mrMesh.GetElementIteratorBegin();
         elem_iter != mrMesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in the vertex mesh
        unsigned elem_index = elem_iter->GetIndex();

        // Get the cell corresponding to this element
        CellPtr p_cell = this->mLocationCellMap[elem_index];
        assert(p_cell);

        if (this->mOutputCellAncestors)
        {
            double ancestor_index = (p_cell->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)p_cell->GetAncestor();
            cell_ancestors[elem_index] = ancestor_index;
        }
        if (this->mOutputCellProliferativeTypes)
        {
            double cell_type = p_cell->GetCellCycleModel()->GetCellProliferativeType();
            cell_types[elem_index] = cell_type;
        }
        if (this->mOutputCellMutationStates)
        {
            double mutation_state = p_cell->GetMutationState()->GetColour();
            cell_mutation_states[elem_index] = mutation_state;
        }
        if (this->mOutputCellAges)
        {
            double age = p_cell->GetAge();
            cell_ages[elem_index] = age;
        }
        if (this->mOutputCellCyclePhases)
        {
            double cycle_phase = p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase();
            cell_cycle_phases[elem_index] = cycle_phase;
        }
        if (this->mOutputCellVolumes)
        {
            double cell_volume = mrMesh.GetVolumeOfElement(elem_index);
            cell_volumes[elem_index] = cell_volume;
        }
        if (CellwiseData<DIM>::Instance()->IsSetUp())
        {
            CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
            unsigned num_variables = p_data->GetNumVariables();

            for (unsigned var=0; var<num_variables; var++)
            {
                cellwise_data[var][elem_index] = p_data->GetValue(p_cell, var);
            }
        }
    }

    if (this->mOutputCellProliferativeTypes)
    {
        mesh_writer.AddCellData("Cell types", cell_types);
    }
    if (this->mOutputCellAncestors)
    {
        mesh_writer.AddCellData("Ancestors", cell_ancestors);
    }
    if (this->mOutputCellMutationStates)
    {
        mesh_writer.AddCellData("Mutation states", cell_mutation_states);
    }
    if (this->mOutputCellAges)
    {
        mesh_writer.AddCellData("Ages", cell_ages);
    }
    if (this->mOutputCellCyclePhases)
    {
        mesh_writer.AddCellData("Cycle phases", cell_cycle_phases);
    }
    if (this->mOutputCellVolumes)
    {
        mesh_writer.AddCellData("Cell volumes", cell_volumes);
    }
    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            std::stringstream data_name;
            data_name << "Cellwise data " << var;
            std::vector<double> cellwise_data_var = cellwise_data[var];
            mesh_writer.AddCellData(data_name.str(), cellwise_data_var);
        }
    }

    mesh_writer.WriteVtkUsingMesh(mrMesh, time.str());
    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << p_time->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << p_time->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractOffLatticeCellPopulation<DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizElementsFile = output_file_handler.OpenOutputFile("results.vizelements");
    mpT1SwapLocationsFile = output_file_handler.OpenOutputFile("T1SwapLocations.dat");
    mpT3SwapLocationsFile = output_file_handler.OpenOutputFile("T3SwapLocations.dat");
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::CloseOutputFiles()
{
    AbstractOffLatticeCellPopulation<DIM>::CloseOutputFiles();
    mpVizElementsFile->close();
    mpT1SwapLocationsFile->close();
    mpT3SwapLocationsFile->close();
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::GenerateCellResultsAndWriteToFiles()
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
void VertexBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<CellRearrangementThreshold>" << mrMesh.GetCellRearrangementThreshold() << "</CellRearrangementThreshold>\n";
    *rParamsFile << "\t\t<T2Threshold>" <<  mrMesh.GetT2Threshold() << "</T2Threshold>\n";
    *rParamsFile << "\t\t<CellRearrangementRatio>" << mrMesh.GetCellRearrangementRatio() << "</CellRearrangementRatio>\n";

    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = mrMesh.GetWidth(rDimension);

    return width;
}

template<unsigned DIM>
std::set<unsigned> VertexBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    return mrMesh.GetNeighbouringNodeIndices(index);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VertexBasedCellPopulation<1>;
template class VertexBasedCellPopulation<2>;
template class VertexBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedCellPopulation)
