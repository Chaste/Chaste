/*

Copyright (C) University of Oxford, 2005-2010

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

#include "VertexBasedTissue.hpp"
#include "CellwiseData.hpp"
#include "VertexMeshWriter.hpp"
#include "Warnings.hpp"

template<unsigned DIM>
VertexBasedTissue<DIM>::VertexBasedTissue(MutableVertexMesh<DIM, DIM>& rMesh,
                                          std::vector<TissueCellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractTissue<DIM>(rCells, locationIndices),
      mrMesh(rMesh),
      mDeleteMesh(deleteMesh)
{
    this->mTissueContainsMesh = true;

    if (validate)
    {
        Validate();
    }
}


template<unsigned DIM>
VertexBasedTissue<DIM>::VertexBasedTissue(MutableVertexMesh<DIM, DIM>& rMesh)
             : mrMesh(rMesh)
{
    this->mTissueContainsMesh = true;
    mDeleteMesh = true;
}


template<unsigned DIM>
VertexBasedTissue<DIM>::~VertexBasedTissue()
{
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}


template<unsigned DIM>
double VertexBasedTissue<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    // Take the average of the cells containing this vertex
    double average_damping_constant = 0.0;

    std::set<unsigned> containing_elements = GetNode(nodeIndex)->rGetContainingElementIndices();
    double temp = 1.0/((double) containing_elements.size());

    for (std::set<unsigned>::iterator iter = containing_elements.begin();
         iter != containing_elements.end();
         ++iter)
    {
        boost::shared_ptr<AbstractCellMutationState> p_state = this->GetCellUsingLocationIndex(*iter)->GetMutationState();
        if (p_state->IsType<WildTypeCellMutationState>())
        {
            average_damping_constant += TissueConfig::Instance()->GetDampingConstantNormal()*temp;
        }
        else
        {
            average_damping_constant += TissueConfig::Instance()->GetDampingConstantMutant()*temp;
        }
    }

    return average_damping_constant;
}


template<unsigned DIM>
MutableVertexMesh<DIM, DIM>& VertexBasedTissue<DIM>::rGetMesh()
{
    return mrMesh;
}


template<unsigned DIM>
const MutableVertexMesh<DIM, DIM>& VertexBasedTissue<DIM>::rGetMesh() const
{
    return mrMesh;
}


template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::GetElement(unsigned elementIndex)
{
    return mrMesh.GetElement(elementIndex);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::GetNumNodes()
{
    return mrMesh.GetNumNodes();
}


template<unsigned DIM>
c_vector<double, DIM> VertexBasedTissue<DIM>::GetLocationOfCellCentre(TissueCellPtr pCell)
{
    return mrMesh.GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}


template<unsigned DIM>
Node<DIM>* VertexBasedTissue<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mrMesh.AddNode(pNewNode);
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mrMesh.SetNode(nodeIndex, rNewLocation);
}


template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::GetElementCorrespondingToCell(TissueCellPtr pCell)
{
    return mrMesh.GetElement(this->mCellLocationMap[pCell.get()]);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::GetNumElements()
{
    return mrMesh.GetNumElements();
}


template<unsigned DIM>
TissueCellPtr VertexBasedTissue<DIM>::AddCell(TissueCellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCellPtr pParentCell)
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
    TissueCellPtr p_created_cell = this->mCells.back();
    this->mLocationCellMap[new_element_index] = p_created_cell;
    this->mCellLocationMap[p_created_cell.get()] = new_element_index;
    return p_created_cell;
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<TissueCellPtr>::iterator it = this->mCells.begin();
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
void VertexBasedTissue<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        // Get damping constant for node
        double damping_const = this->GetDampingConstant(node_index);

        //Get displacement
        c_vector<double, DIM> displacement = dt*rNodeForces[node_index]/damping_const;

        /*
         * If this displacement is greater than the cell rearrangement threshold then this could
         * result in nodes appearing in the middle of other elements which would be bad. If this occurs
         * we restrict the motion to the cell rearrangement threshold and warn the user to use a smaller timesstep? see #1376
         */
        if (norm_2(displacement)>0.5*mrMesh.GetCellRearrangementThreshold())
        {
        	//WARNING("Use a smaller timestep");
        	//displacement *= 0.5*mrMesh.GetCellRearrangementThreshold()/norm_2(displacement);
        }

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + displacement;

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}


template<unsigned DIM>
bool VertexBasedTissue<DIM>::IsCellAssociatedWithADeletedLocation(TissueCellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();;
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    VertexElementMap element_map(mrMesh.GetNumAllElements());

    mrMesh.ReMesh(element_map);

    if (!element_map.IsIdentityMap())
    {
    	// Fix up the mappings between TissueCellPtrs and VertexElements
        std::map<TissueCell*, unsigned> old_map = this->mCellLocationMap;

        this->mCellLocationMap.clear();
        this->mLocationCellMap.clear();

        for (std::list<TissueCellPtr>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            // This shouldn't ever happen, as the cell vector only contains living cells
            unsigned old_elem_index = old_map[(*cell_iter).get()];

            if (element_map.IsDeleted(old_elem_index))
            {
            	/*\todo this is a kludge to remove the cell once a T2Swap occurs this is not included in the dead cells counter.
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
    }

    element_map.ResetToIdentity();

    // Check that each VertexElement has only one TissueCellPtr associated with it in the updated tissue
    Validate();
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::Validate()
{
	// Check each element has only one cell attached.
    std::vector<unsigned> validated_element = std::vector<unsigned>(this->GetNumElements(), 0);

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
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


template<unsigned DIM>
double VertexBasedTissue<DIM>::GetTargetAreaOfCell(const TissueCellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = TissueConfig::Instance()->GetMatureCellTargetArea();

    double cell_age = pCell->GetAge();
    double g1_duration = pCell->GetCellCycleModel()->GetG1Duration();

    // If the cell is differentiated then its G1 duration is infinite
    if (g1_duration == DBL_MAX) // don't use magic number, compare to DBL_MAX
    {
        // This is just for fixed cell cycle models, need to work out how to find the g1 duration
        g1_duration = TissueConfig::Instance()->GetTransitCellG1Duration();
    }

    if (pCell->GetMutationState()->IsType<ApoptoticCellMutationState>())
    {
        // Age of cell when apoptosis begins
        if (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime() < g1_duration)
        {
            cell_target_area *= 0.5*(1 + (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime())/g1_duration);
        }

        // The target area of an apoptotic cell decreases linearly to zero (and past it negative)
        cell_target_area = cell_target_area - 0.5*cell_target_area/(TissueConfig::Instance()->GetApoptosisTime())*(SimulationTime::Instance()->GetTime()-pCell->GetStartOfApoptosisTime());

        // Don't allow a negative target area
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else
    {
        // The target area of a proliferating cell increases linearly from A/2 to A over the course of the G1 phase
        if (cell_age < g1_duration)
        {
            cell_target_area *= 0.5*(1 + cell_age/g1_duration);
        }
    }

    return cell_target_area;
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::WriteResultsToFiles()
{
    AbstractTissue<DIM>::WriteResultsToFiles();

    SimulationTime* p_time = SimulationTime::Instance();

    // Write element data to file
    *mpVizElementsFile << p_time->GetTime() << "\t";

    // Loop over cells and find associated elements so in the same order as the cells in output files
    for (std::list<TissueCellPtr>::iterator cell_iter = this->mCells.begin();
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

#ifdef CHASTE_VTK
    VertexMeshWriter<DIM, DIM> mesh_writer(mDirPath, "results", false);
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
        TissueCellPtr p_cell = this->mLocationCellMap[elem_index];
        assert(p_cell);

        if (TissueConfig::Instance()->GetOutputCellAncestors())
        {
            double ancestor_index = (p_cell->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)p_cell->GetAncestor();
            cell_ancestors[elem_index] = ancestor_index;
        }
        if (TissueConfig::Instance()->GetOutputCellProliferativeTypes())
        {
            double cell_type = p_cell->GetCellCycleModel()->GetCellProliferativeType();
            cell_types[elem_index] = cell_type;
        }
        if (TissueConfig::Instance()->GetOutputCellMutationStates())
        {
            double mutation_state = p_cell->GetMutationState()->GetColour();
            cell_mutation_states[elem_index] = mutation_state;
        }
        if (TissueConfig::Instance()->GetOutputCellAges())
        {
            double age = p_cell->GetAge();
            cell_ages[elem_index] = age;
        }
        if (TissueConfig::Instance()->GetOutputCellCyclePhases())
        {
            double cycle_phase = p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase();
            cell_cycle_phases[elem_index] = cycle_phase;
        }
        if (TissueConfig::Instance()->GetOutputCellVolumes())
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

    if (TissueConfig::Instance()->GetOutputCellProliferativeTypes())
    {
        mesh_writer.AddCellData("Cell types", cell_types);
    }
    if (TissueConfig::Instance()->GetOutputCellAncestors())
    {
        mesh_writer.AddCellData("Ancestors", cell_ancestors);
    }
    if (TissueConfig::Instance()->GetOutputCellMutationStates())
    {
        mesh_writer.AddCellData("Mutation states", cell_mutation_states);
    }
    if (TissueConfig::Instance()->GetOutputCellAges())
    {
        mesh_writer.AddCellData("Ages", cell_ages);
    }
    if (TissueConfig::Instance()->GetOutputCellCyclePhases())
    {
        mesh_writer.AddCellData("Cycle phases", cell_cycle_phases);
    }
    if (TissueConfig::Instance()->GetOutputCellVolumes())
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
    *mpVtkMetaFile << "        <DataSet timestep=\"";
    *mpVtkMetaFile << p_time->GetTimeStepsElapsed();
    *mpVtkMetaFile << "\" group=\"\" part=\"0\" file=\"results_";
    *mpVtkMetaFile << p_time->GetTimeStepsElapsed();
    *mpVtkMetaFile << ".vtu\"/>\n";
#endif //CHASTE_VTK
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractTissue<DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizElementsFile = output_file_handler.OpenOutputFile("results.vizelements");
    mDirPath = rDirectory;
#ifdef CHASTE_VTK
    mpVtkMetaFile = output_file_handler.OpenOutputFile("results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";
#endif //CHASTE_VTK
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::CloseOutputFiles()
{
    AbstractTissue<DIM>::CloseOutputFiles();
    mpVizElementsFile->close();
#ifdef CHASTE_VTK
    *mpVtkMetaFile << "    </Collection>\n";
    *mpVtkMetaFile << "</VTKFile>\n";
    mpVtkMetaFile->close();
#endif //CHASTE_VTK
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::GenerateCellResultsAndWriteToFiles()
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

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        this->GenerateCellResults(this->GetLocationIndexUsingCell(*cell_iter), cell_type_counter, cell_cycle_phase_counter);
    }

    this->WriteCellResultsToFiles(cell_type_counter, cell_cycle_phase_counter);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class VertexBasedTissue<1>;
template class VertexBasedTissue<2>;
template class VertexBasedTissue<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedTissue)
