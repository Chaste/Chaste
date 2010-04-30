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
#include "VertexMeshWriter.hpp"

template<unsigned DIM>
VertexBasedTissue<DIM>::VertexBasedTissue(MutableVertexMesh<DIM, DIM>& rMesh,
                                          std::vector<TissueCell>& rCells,
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
        boost::shared_ptr<AbstractCellMutationState> p_state = this->rGetCellUsingLocationIndex(*iter).GetMutationState();
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
double VertexBasedTissue<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB)
{
    double adhesion_parameter;

    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = TissueConfig::Instance()->GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
    }
    else
    {
        adhesion_parameter = TissueConfig::Instance()->GetNagaiHondaCellCellAdhesionEnergyParameter();
    }
    return adhesion_parameter;
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
c_vector<double, DIM> VertexBasedTissue<DIM>::GetLocationOfCellCentre(TissueCell& rCell)
{
    return mrMesh.GetCentroidOfElement(this->mCellLocationMap[&rCell]);
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
VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::GetElementCorrespondingToCell(TissueCell& rCell)
{
    return mrMesh.GetElement(this->mCellLocationMap[&rCell]);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::GetNumElements()
{
    return mrMesh.GetNumElements();
}


template<unsigned DIM>
TissueCell* VertexBasedTissue<DIM>::AddCell(TissueCell& rNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCell* pParentCell)
{
    // Get the element associated with this cell
    VertexElement<DIM, DIM>* p_element = GetElementCorrespondingToCell(*pParentCell);

    // Divide the element
    unsigned new_element_index;
    if ( norm_2(rCellDivisionVector) < DBL_EPSILON )
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
    this->mCells.push_back(rNewCell);

    // Update location cell map
    TissueCell* p_created_cell = &(this->mCells.back());
    this->mLocationCellMap[new_element_index] = p_created_cell;
    this->mCellLocationMap[p_created_cell] = new_element_index;
    return p_created_cell;
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if (it->IsDead())
        {
            // Remove the element from the mesh
            num_removed++;
            mrMesh.DeleteElementPriorToReMesh(this->mCellLocationMap[&(*it)]);
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

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + displacement;

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}


template<unsigned DIM>
bool VertexBasedTissue<DIM>::IsCellAssociatedWithADeletedLocation(TissueCell& rCell)
{
    return GetElementCorrespondingToCell(rCell)->IsDeleted();;
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    VertexElementMap element_map(mrMesh.GetNumAllElements());

    mrMesh.ReMesh(element_map);

    if (!element_map.IsIdentityMap())
    {
        // Fix up the mappings between TissueCells and VertexElements
        std::map<TissueCell*, unsigned> old_map = this->mCellLocationMap;

        this->mCellLocationMap.clear();
        this->mLocationCellMap.clear();

        for (std::list<TissueCell>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            // This shouldn't ever happen, as the cell vector only contains living cells
            unsigned old_elem_index = old_map[&(*cell_iter)];
            assert(!element_map.IsDeleted(old_elem_index));
            unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);

            this->mLocationCellMap[new_elem_index] = &(*cell_iter);
            this->mCellLocationMap[&(*cell_iter)] = new_elem_index;
        }
    }
    // Check that each VertexElement has a TissueCell associated with it in the updated tissue
    Validate();
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::Validate()
{
    std::vector<bool> validated_element = std::vector<bool>(this->GetNumElements(), false);

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = GetLocationIndexUsingCell(*cell_iter);
        validated_element[elem_index] = true;
    }

    for (unsigned i=0; i<validated_element.size(); i++)
    {
        if (!validated_element[i])
        {
            std::stringstream ss;
            ss << "Element " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}


template<unsigned DIM>
double VertexBasedTissue<DIM>::GetTargetAreaOfCell(const TissueCell& rCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = TissueConfig::Instance()->GetMatureCellTargetArea();

    double cell_age = rCell.GetAge();
    double g1_duration = rCell.GetCellCycleModel()->GetG1Duration();

    // If the cell is differentiated then its G1 duration is infinite
    if (g1_duration == DBL_MAX) // don't use magic number, compare to DBL_MAX
    {
        // This is just for fixed cell cycle models, need to work out how to find the g1 duration
        g1_duration = TissueConfig::Instance()->GetTransitCellG1Duration();
    }

    if (rCell.GetCellProliferativeType() == APOPTOTIC)
    {
        // Age of cell when apoptosis begins
        if (rCell.GetStartOfApoptosisTime() - rCell.GetBirthTime() < g1_duration)
        {
            cell_target_area *= 0.5*(1 + (rCell.GetStartOfApoptosisTime() - rCell.GetBirthTime())/g1_duration);
        }

        // The target area of an apoptotic cell decreases linearly to zero (and past it negative)
        cell_target_area = cell_target_area - cell_target_area/(TissueConfig::Instance()->GetApoptosisTime())*(SimulationTime::Instance()->GetTime()-rCell.GetStartOfApoptosisTime());

        // Don't allow a negative target area
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else
    {
        // In the case of a proliferating cell, the target area increases
        // linearly from A/2 to A over the course of the G1 phase
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
    for (std::list<TissueCell>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);

        // First write the number of Nodes belonging to this VertexElement
        *mpVizElementsFile <<  mrMesh.GetElement(elem_index)->GetNumNodes() << " ";

        // Then write the global index of each Node in this element
        unsigned num_nodes_in_this_element = mrMesh.GetElement(elem_index)->GetNumNodes();
        for (unsigned i=0; i<num_nodes_in_this_element; i++)
        {
            *mpVizElementsFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i) << " ";
        }
    }
    *mpVizElementsFile << "\n";

#ifdef CHASTE_VTK
    VertexMeshWriter<DIM, DIM> mesh_writer(mDirPath, "results", false);
    std::stringstream time;
    time << p_time->GetTimeStepsElapsed();

    std::vector<double> cell_types;
    std::vector<double> cell_ancestors;
    std::vector<double> cell_mutation_states;
    std::vector<double> cell_ages;
    std::vector<double> cell_cycle_phases;
    std::vector<double> cell_areas;

    // Loop over Voronoi elements
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mrMesh.GetElementIteratorBegin();
         elem_iter != mrMesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in the Voronoi tessellation mesh
        unsigned elem_index = elem_iter->GetIndex();

        // Get the cell corresponding to this element
        TissueCell* p_cell = this->mLocationCellMap[elem_index];

        if (TissueConfig::Instance()->GetOutputCellAncestors())
        {
            double ancestor_index = (p_cell->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)p_cell->GetAncestor();
            cell_ancestors.push_back(ancestor_index);
        }
        if (TissueConfig::Instance()->GetOutputCellProliferativeTypes())
        {
            double cell_type = p_cell->GetCellProliferativeType();
            cell_types.push_back(cell_type);
        }
        if (TissueConfig::Instance()->GetOutputCellMutationStates())
        {
            double mutation_state = p_cell->GetMutationState()->GetColour();
            cell_mutation_states.push_back(mutation_state);
        }
        if (TissueConfig::Instance()->GetOutputCellAges())
        {
            double age = p_cell->GetAge();
            cell_ages.push_back(age); 
        }
        if (TissueConfig::Instance()->GetOutputCellCyclePhases())
        {
            double cycle_phase = p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase();
            cell_cycle_phases.push_back(cycle_phase);
        }
        if (TissueConfig::Instance()->GetOutputCellVolumes())
        {
            double cell_area;
            if (DIM==2)
            {
                cell_area = mrMesh.GetAreaOfElement(elem_index);
            }
            else // DIM==3
            {
                cell_area = mrMesh.GetVolumeOfElement(elem_index);
            }
            cell_areas.push_back(cell_area);
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
        mesh_writer.AddCellData("Cell areas", cell_areas);
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
