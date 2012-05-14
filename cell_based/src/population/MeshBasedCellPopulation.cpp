/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "MeshBasedCellPopulation.hpp"
#include "TrianglesMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "CellBasedEventHandler.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Cylindrical2dMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "Exception.hpp"

template<unsigned DIM>
MeshBasedCellPopulation<DIM>::MeshBasedCellPopulation(MutableMesh<DIM, DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh,
                                      bool validate)
    : AbstractCentreBasedCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mpVoronoiTessellation(NULL),
      mDeleteMesh(deleteMesh),
      mUseAreaBasedDampingConstant(false),
      mAreaBasedDampingConstantParameter(0.1),
      mOutputVoronoiData(false),
      mOutputCellPopulationVolumes(false),
      mWriteVtkAsPoints(false)
{
	mpMutableMesh = static_cast<MutableMesh<DIM, DIM>* >(&(this->mrMesh));
    // This must always be true
    assert(this->mCells.size() <= this->mrMesh.GetNumNodes());

    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
MeshBasedCellPopulation<DIM>::MeshBasedCellPopulation(MutableMesh<DIM, DIM>& rMesh)
    : AbstractCentreBasedCellPopulation<DIM>(rMesh)
{
	mpMutableMesh = static_cast<MutableMesh<DIM, DIM>* >(&(this->mrMesh));
    mpVoronoiTessellation = NULL;
    mDeleteMesh = true;
}

template<unsigned DIM>
MeshBasedCellPopulation<DIM>::~MeshBasedCellPopulation()
{
    delete mpVoronoiTessellation;

    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
bool MeshBasedCellPopulation<DIM>::UseAreaBasedDampingConstant()
{
    return mUseAreaBasedDampingConstant;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::SetAreaBasedDampingConstant(bool useAreaBasedDampingConstant)
{
    assert(DIM==2);
    mUseAreaBasedDampingConstant = useAreaBasedDampingConstant;
}

template<unsigned DIM>
unsigned MeshBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mpMutableMesh->AddNode(pNewNode);
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
	static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).SetNode(nodeIndex, rNewLocation, false);
}

template<unsigned DIM>
double MeshBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    double damping_multiplier = AbstractCentreBasedCellPopulation<DIM>::GetDampingConstant(nodeIndex);

    if (mUseAreaBasedDampingConstant)
    {
        /**
         * We use a linear dependence of the form
         *
         * new_damping_const = old_damping_const * (d0+d1*A)
         *
         * where d0, d1 are parameters, A is the cell's area, and old_damping_const
         * is the damping constant if not using mUseAreaBasedDampingConstant
         */

        assert(DIM==2);

        double rest_length = 1.0;
        double d0 = mAreaBasedDampingConstantParameter;

        /**
         * Compute the parameter d1 such that d0+A*d1=1, where A is the equilibrium area
         * of a cell (this is equal to sqrt(3)/4, which is a third of the area of a regular
         * hexagon of edge length 1)
         */
        double d1 = 2.0*(1.0 - d0)/(sqrt(3)*rest_length*rest_length);

        double area_cell = GetVolumeOfVoronoiElement(nodeIndex);

        /**
         * The cell area should not be too large - the next assertion is to avoid
         * getting an infinite cell area, which may occur if area-based viscosity
         * is chosen in the absence of ghost nodes.
         */
        assert(area_cell < 1000);

        damping_multiplier = d0 + area_cell*d1;
    }

    return damping_multiplier;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::Validate()
{
    std::vector<bool> validated_node = std::vector<bool>(this->GetNumNodes(), false);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = GetLocationIndexUsingCell(*cell_iter);
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            EXCEPTION("Node " << i << " does not appear to have a cell associated with it");
        }
    }
}

template<unsigned DIM>
MutableMesh<DIM, DIM>& MeshBasedCellPopulation<DIM>::rGetMesh()
{
    return *mpMutableMesh;
}

template<unsigned DIM>
const MutableMesh<DIM, DIM>& MeshBasedCellPopulation<DIM>::rGetMesh() const
{
    return *mpMutableMesh;
}

template<unsigned DIM>
unsigned MeshBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<CellPtr>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if ((*it)->IsDead())
        {
            // Check if this cell is in a marked spring
            std::vector<const std::pair<CellPtr,CellPtr>*> pairs_to_remove; // Pairs that must be purged
            for (std::set<std::pair<CellPtr,CellPtr> >::iterator it1 = mMarkedSprings.begin();
                 it1 != mMarkedSprings.end();
                 ++it1)
            {
                const std::pair<CellPtr,CellPtr>& r_pair = *it1;

                for (unsigned i=0; i<2; i++)
                {
                    CellPtr p_cell = (i==0 ? r_pair.first : r_pair.second);

                    if (p_cell == *it)
                    {
                        // Remember to purge this spring
                        pairs_to_remove.push_back(&r_pair);
                        break;
                    }
                }
            }

            // Purge any marked springs that contained this cell
            for (std::vector<const std::pair<CellPtr,CellPtr>* >::iterator pair_it = pairs_to_remove.begin();
                 pair_it != pairs_to_remove.end();
                 ++pair_it)
            {
                mMarkedSprings.erase(**pair_it);
            }

            // Remove the node from the mesh
            num_removed++;
            static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).DeleteNodePriorToReMesh(this->GetLocationIndexUsingCell((*it)));

            // Update mappings between cells and location indices
            unsigned location_index_of_removed_node = this->GetLocationIndexUsingCell((*it));
            this->RemoveCellUsingLocationIndex(location_index_of_removed_node, (*it));

            // Update vector of cells
            it = this->mCells.erase(it);
            --it;
        }
    }

    return num_removed;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    NodeMap map(static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).GetNumAllNodes());
    static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).ReMesh(map);

    if (!map.IsIdentityMap())
    {
        UpdateGhostNodesAfterReMesh(map);

        // Update the mappings between cells and location indices
        std::map<Cell*, unsigned> old_map = this->mCellLocationMap;

        // Remove any dead pointers from the maps (needed to avoid archiving errors)
        this->mLocationCellMap.clear();
        this->mCellLocationMap.clear();

        for (std::list<CellPtr>::iterator it = this->mCells.begin();
             it != this->mCells.end();
             ++it)
        {
            unsigned old_node_index = old_map[(*it).get()];

            // This shouldn't ever happen, as the cell vector only contains living cells
            assert(!map.IsDeleted(old_node_index));

            unsigned new_node_index = map.GetNewIndex(old_node_index);
            this->SetCellUsingLocationIndex(new_node_index,*it);
        }

        this->Validate();
    }

    // Purge any marked springs that are no longer springs
    std::vector<const std::pair<CellPtr,CellPtr>*> springs_to_remove;
    for (std::set<std::pair<CellPtr,CellPtr> >::iterator spring_it = mMarkedSprings.begin();
         spring_it != mMarkedSprings.end();
         ++spring_it)
    {
        CellPtr p_cell_1 = spring_it->first;
        CellPtr p_cell_2 = spring_it->second;
        Node<DIM>* p_node_1 = this->GetNodeCorrespondingToCell(p_cell_1);
        Node<DIM>* p_node_2 = this->GetNodeCorrespondingToCell(p_cell_2);

        bool joined = false;

        // For each element containing node1, if it also contains node2 then the cells are joined
        std::set<unsigned> node2_elements = p_node_2->rGetContainingElementIndices();
        for (typename Node<DIM>::ContainingElementIterator elem_iter = p_node_1->ContainingElementsBegin();
             elem_iter != p_node_1->ContainingElementsEnd();
             ++elem_iter)
        {
            if (node2_elements.find(*elem_iter) != node2_elements.end())
            {
                joined = true;
                break;
            }
        }

        // If no longer joined, remove this spring from the set
        if (!joined)
        {
            springs_to_remove.push_back(&(*spring_it));
        }
    }

    // Remove any springs necessary
    for (std::vector<const std::pair<CellPtr,CellPtr>* >::iterator spring_it = springs_to_remove.begin();
         spring_it != springs_to_remove.end();
         ++spring_it)
    {
        mMarkedSprings.erase(**spring_it);
    }

    // Tessellate if needed
    TessellateIfNeeded();

    static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).SetMeshHasChangedSinceLoading();
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::TessellateIfNeeded()
{
    if (DIM==2 || DIM==3)
    {
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::TESSELLATION);
        if (mUseAreaBasedDampingConstant || mOutputVoronoiData || mOutputCellPopulationVolumes || this->mOutputCellVolumes)
        {
            CreateVoronoiTessellation();
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::TESSELLATION);
    }
}

template<unsigned DIM>
Node<DIM>* MeshBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned MeshBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumAllNodes();
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{
}

template<unsigned DIM>
CellPtr MeshBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    assert(pNewCell);
    assert(pParentCell);

    // Add new cell to cell population
    CellPtr p_created_cell = AbstractCentreBasedCellPopulation<DIM>::AddCell(pNewCell, rCellDivisionVector, pParentCell);
    assert(p_created_cell == pNewCell);

    // Mark spring between parent cell and new cell
    std::pair<CellPtr,CellPtr> cell_pair = CreateCellPair(pParentCell, p_created_cell);
    MarkSpring(cell_pair);

    // Return pointer to new cell
    return p_created_cell;
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractCentreBasedCellPopulation<DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizElementsFile = output_file_handler.OpenOutputFile("results.vizelements");

    if (mOutputVoronoiData)
    {
        mpVoronoiFile = output_file_handler.OpenOutputFile("voronoi.dat");
    }
    if (mOutputCellPopulationVolumes)
    {
        mpCellPopulationVolumesFile = output_file_handler.OpenOutputFile("cellpopulationareas.dat");
    }
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::CloseOutputFiles()
{
    AbstractCentreBasedCellPopulation<DIM>::CloseOutputFiles();

    mpVizElementsFile->close();

    if (mOutputVoronoiData)
    {
        mpVoronoiFile->close();
    }
    if (mOutputCellPopulationVolumes)
    {
        mpCellPopulationVolumesFile->close();
    }
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::WriteResultsToFiles()
{
    if (SimulationTime::Instance()->GetTimeStepsElapsed() == 0 && this->mpVoronoiTessellation == NULL)
    {
        TessellateIfNeeded();//Update isn't run on time-step zero
    }

    AbstractCentreBasedCellPopulation<DIM>::WriteResultsToFiles();

    // Write element data to file
    *mpVizElementsFile << SimulationTime::Instance()->GetTime() << "\t";

    for (typename MutableMesh<DIM,DIM>::ElementIterator elem_iter = static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).GetElementIteratorBegin();
         elem_iter != static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).GetElementIteratorEnd();
         ++elem_iter)
    {
        bool element_contains_dead_cells_or_deleted_nodes = false;

        // Hack that covers the case where the element contains a node that is associated with a cell that has just been killed (#1129)
        for (unsigned i=0; i<DIM+1; i++)
        {
            unsigned node_index = elem_iter->GetNodeGlobalIndex(i);

            if (this->GetNode(node_index)->IsDeleted())
            {
                element_contains_dead_cells_or_deleted_nodes = true;
                break;
            }
            else if (this->IsCellAttachedToLocationIndex(node_index))
            {
                if (this->GetCellUsingLocationIndex(node_index)->IsDead())
                {
                    element_contains_dead_cells_or_deleted_nodes = true;
                    break;
                }
            }
        }
        if (!element_contains_dead_cells_or_deleted_nodes)
        {
            for (unsigned i=0; i<DIM+1; i++)
            {
                *mpVizElementsFile << elem_iter->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
    *mpVizElementsFile << "\n";

    // Write data to file.
    if (mOutputVoronoiData)
    {
        WriteVoronoiResultsToFile();
    }
    if (mOutputCellPopulationVolumes)
    {
        WriteCellPopulationVolumeResultsToFile();
    }

}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    // Write time to file
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();

    unsigned num_points = GetNumNodes();
    if (!mWriteVtkAsPoints && (mpVoronoiTessellation != NULL))
    {
        num_points = mpVoronoiTessellation->GetNumElements();
    }

    std::vector<double> cell_types(num_points);
    std::vector<double> cell_ancestors(num_points);
    std::vector<double> cell_mutation_states(num_points);
    std::vector<double> cell_ages(num_points);
    std::vector<double> cell_cycle_phases(num_points);
    std::vector<std::vector<double> > cellwise_data;

    unsigned num_cell_data_items = 0;
    if (this->Begin()->HasCellData())
    {
        //We assume that the first cell is representative of all cells
        num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    }

    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cellwise_data_var(num_points);
        cellwise_data.push_back(cellwise_data_var);
    }

    if (mWriteVtkAsPoints)
    {
        VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);

        // Loop over cells
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
             cell_iter != this->End();
             ++cell_iter)
        {
            // Get the node index corresponding to this cell
            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

            // Get this cell-cycle model
            AbstractCellCycleModel* p_model = cell_iter->GetCellCycleModel();

            if (this->mOutputCellAncestors)
            {
                double ancestor_index = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_iter->GetAncestor();
                cell_ancestors[node_index] = ancestor_index;
            }
            if (this->mOutputCellProliferativeTypes)
            {
                cell_types[node_index] = p_model->GetCellProliferativeType();
            }
            if (this->mOutputCellMutationStates)
            {
                cell_mutation_states[node_index] = cell_iter->GetMutationState()->GetColour();
            }
            if (this->mOutputCellAges)
            {
                cell_ages[node_index] = cell_iter->GetAge();
            }
            if (this->mOutputCellCyclePhases)
            {
                cell_cycle_phases[node_index] = p_model->GetCurrentCellCyclePhase();
            }
            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cellwise_data[var][node_index] = cell_iter->GetCellData()->GetItem(var);
            }
        }

        if (this->mOutputCellProliferativeTypes)
        {
            mesh_writer.AddPointData("Cell types", cell_types);
        }
        if (this->mOutputCellAncestors)
        {
            mesh_writer.AddPointData("Ancestors", cell_ancestors);
        }
        if (this->mOutputCellMutationStates)
        {
            mesh_writer.AddPointData("Mutation states", cell_mutation_states);
        }
        if (this->mOutputCellAges)
        {
            mesh_writer.AddPointData("Ages", cell_ages);
        }
        if (this->mOutputCellCyclePhases)
        {
            mesh_writer.AddPointData("Cycle phases", cell_cycle_phases);
        }
        if (num_cell_data_items > 0)
        {
            for (unsigned var=0; var<cellwise_data.size(); var++)
            {
                std::stringstream data_name;
                data_name << "Cellwise data " << var;
                std::vector<double> cellwise_data_var = cellwise_data[var];
                mesh_writer.AddPointData(data_name.str(), cellwise_data_var);
            }
        }

        {
            // Make a copy of the nodes in a disposable mesh for writing
            std::vector<Node<DIM>* > nodes;
            for (unsigned index=0; index<this->mrMesh.GetNumNodes(); index++)
            {
                Node<DIM>* p_node = this->mrMesh.GetNode(index);
                nodes.push_back(p_node);
            }

            NodesOnlyMesh<DIM> mesh;
            mesh.ConstructNodesWithoutMesh(nodes);
            mesh_writer.WriteFilesUsingMesh(mesh);
        }

        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
    else if (mpVoronoiTessellation != NULL)
    {
        VertexMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results", false);
        std::vector<double> cell_volumes(num_points);

        // Loop over elements of mpVoronoiTessellation
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
             elem_iter != mpVoronoiTessellation->GetElementIteratorEnd();
             ++elem_iter)
        {
            // Get index of this element in mpVoronoiTessellation
            unsigned elem_index = elem_iter->GetIndex();

            // Get the index of the corresponding node in mrMesh
            unsigned node_index = mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

            // There should be no ghost nodes
            assert(!this->IsGhostNode(node_index));

            // Get the cell corresponding to this element
            CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

            // Get this cell-cycle model
            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();

            if (this->mOutputCellAncestors)
            {
                double ancestor_index = (p_cell->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)p_cell->GetAncestor();
                cell_ancestors[elem_index] = ancestor_index;
            }
            if (this->mOutputCellProliferativeTypes)
            {
                cell_types[elem_index] = p_model->GetCellProliferativeType();
            }
            if (this->mOutputCellMutationStates)
            {
                cell_mutation_states[elem_index] = p_cell->GetMutationState()->GetColour();
            }
            if (this->mOutputCellAges)
            {
                cell_ages[elem_index] = p_cell->GetAge();
            }
            if (this->mOutputCellCyclePhases)
            {
                cell_cycle_phases[elem_index] = p_model->GetCurrentCellCyclePhase();
            }
            if (this->mOutputCellVolumes)
            {
                cell_volumes[elem_index] = mpVoronoiTessellation->GetVolumeOfElement(elem_index);
            }
            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cellwise_data[var][elem_index] = p_cell->GetCellData()->GetItem(var);
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
        if (num_cell_data_items > 0)
        {
            for (unsigned var=0; var<cellwise_data.size(); var++)
            {
                std::stringstream data_name;
                data_name << "Cellwise data " << var;
                std::vector<double> cellwise_data_var = cellwise_data[var];
                mesh_writer.AddCellData(data_name.str(), cellwise_data_var);
            }
        }

        mesh_writer.WriteVtkUsingMesh(*mpVoronoiTessellation, time.str());
        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
#endif //CHASTE_VTK
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::WriteVoronoiResultsToFile()
{
    assert(DIM==2 || DIM==3);
    assert(mpVoronoiTessellation != NULL);
    // Write time to file
    *mpVoronoiFile << SimulationTime::Instance()->GetTime() << " ";

    // Loop over elements of mpVoronoiTessellation
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
         elem_iter != mpVoronoiTessellation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in mpVoronoiTessellation
        unsigned elem_index = elem_iter->GetIndex();

        // Get the index of the corresponding node in mrMesh
        unsigned node_index = mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        // Write node index and location to file
        *mpVoronoiFile << node_index << " ";
        c_vector<double, DIM> node_location = this->GetNode(node_index)->rGetLocation();
        for (unsigned i=0; i<DIM; i++)
        {
            *mpVoronoiFile << node_location[i] << " ";
        }

        double cell_volume = mpVoronoiTessellation->GetVolumeOfElement(elem_index);
        double cell_surface_area = mpVoronoiTessellation->GetSurfaceAreaOfElement(elem_index);
        *mpVoronoiFile << cell_volume << " " << cell_surface_area << " ";
    }
    *mpVoronoiFile << "\n";
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::WriteCellPopulationVolumeResultsToFile()
{
    assert(DIM==2 || DIM==3);

    // Write time to file
    *mpCellPopulationVolumesFile << SimulationTime::Instance()->GetTime() << " ";
    assert (this->mpVoronoiTessellation != NULL);

    // Don't use the Voronoi tessellation to calculate the total area of the mesh because it gives huge areas for boundary cells
    double total_area = static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).GetVolume();
    double apoptotic_area = 0.0;

    // Loop over elements of mpVoronoiTessellation
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
         elem_iter != mpVoronoiTessellation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in mpVoronoiTessellation
        unsigned elem_index = elem_iter->GetIndex();

        // Get the index of the corresponding node in mrMesh
        unsigned node_index = mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        // Discount ghost nodes
        if (!this->IsGhostNode(node_index))
        {
            // Get the cell corresponding to this node
            CellPtr p_cell =  this->GetCellUsingLocationIndex(node_index);

            // Only bother calculating the area/volume of apoptotic cells
            bool cell_is_apoptotic = p_cell->HasCellProperty<ApoptoticCellProperty>();
            if (cell_is_apoptotic)
            {
                double cell_volume = mpVoronoiTessellation->GetVolumeOfElement(elem_index);
                apoptotic_area += cell_volume;
            }
        }
    }
    *mpCellPopulationVolumesFile << total_area << " " << apoptotic_area << "\n";
}

template<unsigned DIM>
double MeshBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Ensure that the Voronoi tessellation exists
    if (mpVoronoiTessellation == NULL)
    {
        CreateVoronoiTessellation();
    }

    // Get the node index corresponding to this cell
    unsigned node_index = this->GetLocationIndexUsingCell(pCell);

    // Get the element index of the Voronoi tessellation corresponding to this node index
    unsigned element_index = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(node_index);
    
    // Get the cell's volume from the Voronoi tessellation
    double cell_volume = mpVoronoiTessellation->GetVolumeOfElement(element_index);

    return cell_volume;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::WriteCellVolumeResultsToFile()
{
    assert (mpVoronoiTessellation != NULL);
    assert(DIM==2 || DIM==3);

    // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    ///\todo It would be simpler to merely iterate over the cell population here

    // Loop over elements of mpVoronoiTessellation
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
         elem_iter != mpVoronoiTessellation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in mpVoronoiTessellation
        unsigned elem_index = elem_iter->GetIndex();

        // Get the index of the corresponding node in mrMesh
        unsigned node_index = mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        // Discount ghost nodes
        if (!this->IsGhostNode(node_index))
        {
            // Write node index to file
            *(this->mpCellVolumesFile) << node_index << " ";

            // Get the cell corresponding to this node
            CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

            // Write cell ID to file
            unsigned cell_index = p_cell->GetCellId();
            *(this->mpCellVolumesFile) << cell_index << " ";

            // Write node location to file
            c_vector<double, DIM> node_location = this->GetNode(node_index)->rGetLocation();
            for (unsigned i=0; i<DIM; i++)
            {
                *(this->mpCellVolumesFile) << node_location[i] << " ";
            }

            // Write cell volume (in 3D) or area (in 2D) to file
            double cell_volume = this->GetVolumeOfCell(p_cell);
            *(this->mpCellVolumesFile) << cell_volume << " ";
        }
    }
    *(this->mpCellVolumesFile) << "\n";

}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::SetWriteVtkAsPoints(bool writeVtkAsPoints)
{
    mWriteVtkAsPoints = writeVtkAsPoints;
}

template<unsigned DIM>
bool MeshBasedCellPopulation<DIM>::GetWriteVtkAsPoints()
{
    return mWriteVtkAsPoints;
}

//////////////////////////////////////////////////////////////////////////////
//                          Spring iterator class                           //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
Node<DIM>* MeshBasedCellPopulation<DIM>::SpringIterator::GetNodeA()
{
    return mEdgeIter.GetNodeA();
}

template<unsigned DIM>
Node<DIM>* MeshBasedCellPopulation<DIM>::SpringIterator::GetNodeB()
{
    return mEdgeIter.GetNodeB();
}

template<unsigned DIM>
CellPtr MeshBasedCellPopulation<DIM>::SpringIterator::GetCellA()
{
    assert((*this) != mrCellPopulation.SpringsEnd());
    return mrCellPopulation.GetCellUsingLocationIndex(mEdgeIter.GetNodeA()->GetIndex());
}

template<unsigned DIM>
CellPtr MeshBasedCellPopulation<DIM>::SpringIterator::GetCellB()
{
    assert((*this) != mrCellPopulation.SpringsEnd());
    return mrCellPopulation.GetCellUsingLocationIndex(mEdgeIter.GetNodeB()->GetIndex());
}

template<unsigned DIM>
bool MeshBasedCellPopulation<DIM>::SpringIterator::operator!=(const MeshBasedCellPopulation<DIM>::SpringIterator& rOther)
{
    return (mEdgeIter != rOther.mEdgeIter);
}

template<unsigned DIM>
typename MeshBasedCellPopulation<DIM>::SpringIterator& MeshBasedCellPopulation<DIM>::SpringIterator::operator++()
{
    bool edge_is_ghost = false;

    do
    {
        ++mEdgeIter;
        if (*this != mrCellPopulation.SpringsEnd())
        {
            bool a_is_ghost = mrCellPopulation.IsGhostNode(mEdgeIter.GetNodeA()->GetIndex());
            bool b_is_ghost = mrCellPopulation.IsGhostNode(mEdgeIter.GetNodeB()->GetIndex());

            edge_is_ghost = (a_is_ghost || b_is_ghost);
        }
    }
    while (*this!=mrCellPopulation.SpringsEnd() && edge_is_ghost);

    return (*this);
}

template<unsigned DIM>
MeshBasedCellPopulation<DIM>::SpringIterator::SpringIterator(
            MeshBasedCellPopulation<DIM>& rCellPopulation,
            typename MutableMesh<DIM,DIM>::EdgeIterator edgeIter)
    : mrCellPopulation(rCellPopulation),
      mEdgeIter(edgeIter)
{
    if (mEdgeIter!=static_cast<MutableMesh<DIM, DIM>*>(&(this->mrCellPopulation.mrMesh))->EdgesEnd())
    {
        bool a_is_ghost = mrCellPopulation.IsGhostNode(mEdgeIter.GetNodeA()->GetIndex());
        bool b_is_ghost = mrCellPopulation.IsGhostNode(mEdgeIter.GetNodeB()->GetIndex());

        if (a_is_ghost || b_is_ghost)
        {
            ++(*this);
        }
    }
}

template<unsigned DIM>
typename MeshBasedCellPopulation<DIM>::SpringIterator MeshBasedCellPopulation<DIM>::SpringsBegin()
{
    return SpringIterator(*this, static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).EdgesBegin());
}

template<unsigned DIM>
typename MeshBasedCellPopulation<DIM>::SpringIterator MeshBasedCellPopulation<DIM>::SpringsEnd()
{
    return SpringIterator(*this, static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).EdgesEnd());
}

/**
 *
 */
template<>
void MeshBasedCellPopulation<2>::CreateVoronoiTessellation()
{
    delete mpVoronoiTessellation;

    // Check if the mesh associated with this cell population is periodic
    bool is_mesh_periodic = false;
    if (dynamic_cast<Cylindrical2dMesh*>(&mrMesh))
    {
        is_mesh_periodic = true;
    }

    mpVoronoiTessellation = new VertexMesh<2, 2>(static_cast<MutableMesh<2, 2> &>((this->mrMesh)), is_mesh_periodic);
}

/**
 * The cylindrical mesh is only defined in 2D, hence there is
 * a separate definition for this method in 3D, which doesn't have the capability
 * of dealing with periodic boundaries in 3D. This is /todo #1374.
 */
template<>
void MeshBasedCellPopulation<3>::CreateVoronoiTessellation()
{
    delete mpVoronoiTessellation;
    mpVoronoiTessellation = new VertexMesh<3, 3>(static_cast<MutableMesh<3, 3> &>((this->mrMesh)));
}

/**
 * The VoronoiTessellation class is only defined in 2D or 3D, hence there
 * are two definitions to this method (one templated and one not).
 */
template<>
void MeshBasedCellPopulation<1>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}

template<unsigned DIM>
VertexMesh<DIM, DIM>* MeshBasedCellPopulation<DIM>::GetVoronoiTessellation()
{
    assert(mpVoronoiTessellation!=NULL);
    return mpVoronoiTessellation;
}

template<unsigned DIM>
double MeshBasedCellPopulation<DIM>::GetVolumeOfVoronoiElement(unsigned index)
{
    unsigned element_index = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
    double volume = mpVoronoiTessellation->GetVolumeOfElement(element_index);
    return volume;
}

template<unsigned DIM>
double MeshBasedCellPopulation<DIM>::GetSurfaceAreaOfVoronoiElement(unsigned index)
{
    unsigned element_index = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
    double surface_area = mpVoronoiTessellation->GetSurfaceAreaOfElement(element_index);
    return surface_area;
}

template<unsigned DIM>
double MeshBasedCellPopulation<DIM>::GetVoronoiEdgeLength(unsigned index1, unsigned index2)
{
    unsigned element_index1 = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index1);
    unsigned element_index2 = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index2);
    try
    {
        double edge_length = mpVoronoiTessellation->GetEdgeLength(element_index1, element_index2);
        return edge_length;
    }
    catch (Exception& e)
    {
        // The edge was between two (potentially infinite) cells on the boundary of the mesh
        EXCEPTION("Spring iterator tried to calculate interaction between degenerate cells on the boundary of the mesh.  Have you set ghost layers correctly?");
    }
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::CheckCellPointers()
{
    bool res = true;
    for (std::list<CellPtr>::iterator it=this->mCells.begin();
         it!=this->mCells.end();
         ++it)
    {
        CellPtr p_cell = *it;
        assert(p_cell);
        AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
        assert(p_model);

        // Check cell exists in cell population
        unsigned node_index = this->GetLocationIndexUsingCell(p_cell);
        std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
        CellPtr p_cell_in_cell_population = this->GetCellUsingLocationIndex(node_index);
#define COVERAGE_IGNORE //Debugging code.  Shouldn't fail under normal conditions
        if (p_cell_in_cell_population != p_cell)
        {
            std::cout << "  Mismatch with cell population" << std::endl << std::flush;
            res = false;
        }

        // Check model links back to cell
        if (p_model->GetCell() != p_cell)
        {
            std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
            res = false;
        }
    }
    assert(res);
#undef COVERAGE_IGNORE

    res = true;
    for (std::set<std::pair<CellPtr,CellPtr> >::iterator it1 = mMarkedSprings.begin();
         it1 != mMarkedSprings.end();
         ++it1)
    {
        const std::pair<CellPtr,CellPtr>& r_pair = *it1;

        for (unsigned i=0; i<2; i++)
        {
            CellPtr p_cell = (i==0 ? r_pair.first : r_pair.second);

            assert(p_cell);
            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();
            assert(p_model);
            unsigned node_index = this->GetLocationIndexUsingCell(p_cell);
            std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;

#define COVERAGE_IGNORE //Debugging code.  Shouldn't fail under normal conditions
            // Check cell is alive
            if (p_cell->IsDead())
            {
                std::cout << "  Cell is dead" << std::endl << std::flush;
                res = false;
            }

            // Check cell exists in cell population
            CellPtr p_cell_in_cell_population = this->GetCellUsingLocationIndex(node_index);
            if (p_cell_in_cell_population != p_cell)
            {
                std::cout << "  Mismatch with cell population" << std::endl << std::flush;
                res = false;
            }

            // Check model links back to cell
            if (p_model->GetCell() != p_cell)
            {
                std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
                res = false;
            }
        }
#undef COVERAGE_IGNORE
    }
    assert(res);
}

template<unsigned DIM>
std::pair<CellPtr,CellPtr> MeshBasedCellPopulation<DIM>::CreateCellPair(CellPtr pCell1, CellPtr pCell2)
{
    assert(pCell1);
    assert(pCell2);

    std::pair<CellPtr,CellPtr> cell_pair;

    if (pCell1->GetCellId() < pCell2->GetCellId())
    {
        cell_pair.first = pCell1;
        cell_pair.second = pCell2;
    }
    else
    {
        cell_pair.first = pCell2;
        cell_pair.second = pCell1;
    }
    return cell_pair;
}

template<unsigned DIM>
bool MeshBasedCellPopulation<DIM>::IsMarkedSpring(const std::pair<CellPtr,CellPtr>& rCellPair)
{
    // the pair should be ordered like this (CreateCellPair will ensure this)
    assert(rCellPair.first->GetCellId() < rCellPair.second->GetCellId());

    return mMarkedSprings.find(rCellPair) != mMarkedSprings.end();
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::MarkSpring(std::pair<CellPtr,CellPtr>& rCellPair)
{
    // the pair should be ordered like this (CreateCellPair will ensure this)
    assert(rCellPair.first->GetCellId() < rCellPair.second->GetCellId());

    mMarkedSprings.insert(rCellPair);
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::UnmarkSpring(std::pair<CellPtr,CellPtr>& rCellPair)
{
    // the pair should be ordered like this (CreateCellPair will ensure this)
    assert(rCellPair.first->GetCellId() < rCellPair.second->GetCellId());

    mMarkedSprings.erase(rCellPair);
}

template<unsigned DIM>
double MeshBasedCellPopulation<DIM>::GetAreaBasedDampingConstantParameter()
{
    return mAreaBasedDampingConstantParameter;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::SetAreaBasedDampingConstantParameter(double areaBasedDampingConstantParameter)
{
    assert(areaBasedDampingConstantParameter >= 0.0);
    mAreaBasedDampingConstantParameter = areaBasedDampingConstantParameter;
}

template<unsigned DIM>
bool MeshBasedCellPopulation<DIM>::GetOutputVoronoiData()
{
    return mOutputVoronoiData;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::SetOutputVoronoiData(bool outputVoronoiData)
{
    mOutputVoronoiData = outputVoronoiData;
}

template<unsigned DIM>
bool MeshBasedCellPopulation<DIM>::GetOutputCellPopulationVolumes()
{
    return mOutputCellPopulationVolumes;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::SetOutputCellPopulationVolumes(bool outputCellPopulationVolumes)
{
    mOutputCellPopulationVolumes = outputCellPopulationVolumes;
}

template<unsigned DIM>
void MeshBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<UseAreaBasedDampingConstant>" << mUseAreaBasedDampingConstant << "</UseAreaBasedDampingConstant>\n";
    *rParamsFile << "\t\t<AreaBasedDampingConstantParameter>" <<  mAreaBasedDampingConstantParameter << "</AreaBasedDampingConstantParameter>\n";
    *rParamsFile << "\t\t<OutputVoronoiData>" <<  mOutputVoronoiData << "</OutputVoronoiData>\n";
    *rParamsFile << "\t\t<OutputCellPopulationVolumes>" << mOutputCellPopulationVolumes << "</OutputCellPopulationVolumes>\n";

    // Call method on direct parent class
    AbstractCentreBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
double MeshBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = this->mrMesh.GetWidth(rDimension);
    return width;
}

template<unsigned DIM>
std::set<unsigned> MeshBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    // Get pointer to this node
    Node<DIM>* p_node = this->mrMesh.GetNode(index);

    // Loop over containing elements
    std::set<unsigned> neighbouring_node_indices;
    for (typename Node<DIM>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
         elem_iter != p_node->ContainingElementsEnd();
         ++elem_iter)
    {
        // Get pointer to this containing element
        Element<DIM,DIM>* p_element = static_cast<MutableMesh<DIM, DIM>&>((this->mrMesh)).GetElement(*elem_iter);

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

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MeshBasedCellPopulation<1>;
template class MeshBasedCellPopulation<2>;
template class MeshBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedCellPopulation)
