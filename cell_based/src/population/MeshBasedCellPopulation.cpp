/*

Copyright (c) 2005-2013, University of Oxford.
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


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::MeshBasedCellPopulation(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh,
                                      bool validate)
    : AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>(rMesh, rCells, locationIndices),
      mpVoronoiTessellation(NULL),
      mDeleteMesh(deleteMesh),
      mUseAreaBasedDampingConstant(false),
      mAreaBasedDampingConstantParameter(0.1),
      mOutputVoronoiData(false),
      mOutputCellPopulationVolumes(false),
      mWriteVtkAsPoints(false),
      mOutputMeshInVtk(false),
      mHasVariableRestLength(false)
{
    mpMutableMesh = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>* >(&(this->mrMesh));
    // This must always be true
    assert(this->mCells.size() <= this->mrMesh.GetNumNodes());

    if (validate)
    {
        Validate();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::MeshBasedCellPopulation(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    : AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>(rMesh)
{
    mpMutableMesh = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>* >(&(this->mrMesh));
    mpVoronoiTessellation = NULL;
    mDeleteMesh = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::~MeshBasedCellPopulation()
{
    delete mpVoronoiTessellation;

    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::UseAreaBasedDampingConstant()
{
    return mUseAreaBasedDampingConstant;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetAreaBasedDampingConstant(bool useAreaBasedDampingConstant)
{
    assert(SPACE_DIM==2);
    mUseAreaBasedDampingConstant = useAreaBasedDampingConstant;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    return mpMutableMesh->AddNode(pNewNode);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM>& rNewLocation)
{
    static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).SetNode(nodeIndex, rNewLocation, false);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetDampingConstant(unsigned nodeIndex)
{
    double damping_multiplier = AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetDampingConstant(nodeIndex);

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

        assert(SPACE_DIM==2);

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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::Validate()
{
    std::vector<bool> validated_node = std::vector<bool>(this->GetNumNodes(), false);

    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMesh<ELEMENT_DIM,SPACE_DIM>& MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::rGetMesh()
{
    return *mpMutableMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const MutableMesh<ELEMENT_DIM,SPACE_DIM>& MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::rGetMesh() const
{
    return *mpMutableMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::RemoveDeadCells()
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
            for (std::set<std::pair<CellPtr,CellPtr> >::iterator it1 = this->mMarkedSprings.begin();
                 it1 != this->mMarkedSprings.end();
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
                this->mMarkedSprings.erase(**pair_it);
            }

            // Remove the node from the mesh
            num_removed++;
            static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).DeleteNodePriorToReMesh(this->GetLocationIndexUsingCell((*it)));

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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::Update(bool hasHadBirthsOrDeaths)
{
    NodeMap map(static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).GetNumAllNodes());
    static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).ReMesh(map);

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
    for (std::set<std::pair<CellPtr,CellPtr> >::iterator spring_it = this->mMarkedSprings.begin();
         spring_it != this->mMarkedSprings.end();
         ++spring_it)
    {
        CellPtr p_cell_1 = spring_it->first;
        CellPtr p_cell_2 = spring_it->second;
        Node<SPACE_DIM>* p_node_1 = this->GetNodeCorrespondingToCell(p_cell_1);
        Node<SPACE_DIM>* p_node_2 = this->GetNodeCorrespondingToCell(p_cell_2);

        bool joined = false;

        // For each element containing node1, if it also contains node2 then the cells are joined
        std::set<unsigned> node2_elements = p_node_2->rGetContainingElementIndices();
        for (typename Node<SPACE_DIM>::ContainingElementIterator elem_iter = p_node_1->ContainingElementsBegin();
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
        this->mMarkedSprings.erase(**spring_it);
    }

    // Tessellate if needed
    TessellateIfNeeded();

    static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).SetMeshHasChangedSinceLoading();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::TessellateIfNeeded()
{
     if ((SPACE_DIM==2 || SPACE_DIM==3)&&(ELEMENT_DIM==SPACE_DIM))
    {
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::TESSELLATION);
        if (mUseAreaBasedDampingConstant || mOutputVoronoiData || mOutputCellPopulationVolumes || this->mOutputCellVolumes)
        {
            CreateVoronoiTessellation();
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::TESSELLATION);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::DivideLongSprings(double springDivisionThreshold)
{
    // Only implemented for 2D elements
    assert(ELEMENT_DIM==2);

    std::vector<c_vector<unsigned, 5> > new_nodes;
    new_nodes = rGetMesh().SplitLongEdges(springDivisionThreshold);

    // Add new cells onto new nodes
    for (unsigned index=0; index<new_nodes.size(); index++)
    {
        // Copy the cell attached to one of the neighbouring nodes onto the new node
        unsigned new_node_index = new_nodes[index][0];
        unsigned node_a_index = new_nodes[index][1];
        unsigned node_b_index = new_nodes[index][2];

         CellPtr p_neighbour_cell = this->GetCellUsingLocationIndex(node_a_index);

        // Create copy of cell property collection to modify for daughter cell
        CellPropertyCollection daughter_property_collection = p_neighbour_cell->rGetCellPropertyCollection();

        // Remove the CellId from the daughter cell a new one will be assigned in the constructor
        daughter_property_collection.RemoveProperty<CellId>();

        CellPtr p_new_cell(new Cell(p_neighbour_cell->GetMutationState(),
                                      p_neighbour_cell->GetCellCycleModel()->CreateCellCycleModel(),
                                      false,
                                      daughter_property_collection));

        // Add new cell to cell population
        this->mCells.push_back(p_new_cell);
        this->AddCellUsingLocationIndex(new_node_index,p_new_cell);

        // Update rest lengths

        // remove old node pair // note node_a_index < node_b_index
        std::pair<unsigned,unsigned> node_pair = CreateOrderedPair(node_a_index, node_b_index);
        double old_rest_length  = mSpringRestLengths[node_pair];

        std::map<std::pair<unsigned,unsigned>, double>::iterator  iter = mSpringRestLengths.find(node_pair);
        mSpringRestLengths.erase(iter);

        //Add new pairs
        node_pair = CreateOrderedPair(node_a_index, new_node_index);
        mSpringRestLengths[node_pair]= old_rest_length/2;

        node_pair = CreateOrderedPair(node_b_index, new_node_index);
        mSpringRestLengths[node_pair]= old_rest_length/2;
        // If necessary add other new spring rest lengths
        for (unsigned pair_index = 3; pair_index < 5; pair_index++)
        {
            unsigned other_node_index = new_nodes[index][pair_index];

            if (other_node_index != UNSIGNED_UNSET)
            {
                node_pair = CreateOrderedPair(other_node_index, new_node_index);
                double new_rest_length = norm_2(rGetMesh().GetVectorFromAtoB(rGetMesh().GetNode(new_node_index)->rGetLocation(),
                                                    rGetMesh().GetNode(other_node_index)->rGetLocation()));

                mSpringRestLengths[node_pair] = new_rest_length;
            }

        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumAllNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPtr MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AddCell(CellPtr pNewCell, const c_vector<double,SPACE_DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    assert(pNewCell);
    assert(pParentCell);

    // Add new cell to cell population
    CellPtr p_created_cell = AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AddCell(pNewCell, rCellDivisionVector, pParentCell);
    assert(p_created_cell == pNewCell);

    // Mark spring between parent cell and new cell
    std::pair<CellPtr,CellPtr> cell_pair = this->CreateCellPair(pParentCell, p_created_cell);
    this->MarkSpring(cell_pair);

    // Return pointer to new cell
    return p_created_cell;
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::CloseOutputFiles()
{
    AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::CloseOutputFiles();

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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteResultsToFiles()
{
    if (SimulationTime::Instance()->GetTimeStepsElapsed() == 0 && this->mpVoronoiTessellation == NULL)
    {
        TessellateIfNeeded();//Update isn't run on time-step zero
    }

    AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteResultsToFiles();

    // Write element data to file
    *mpVizElementsFile << SimulationTime::Instance()->GetTime() << "\t";

    for (typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).GetElementIteratorBegin();
         elem_iter != static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).GetElementIteratorEnd();
         ++elem_iter)
    {
        bool element_contains_dead_cells_or_deleted_nodes = false;

        // Hack that covers the case where the element contains a node that is associated with a cell that has just been killed (#1129)
        for (unsigned i=0; i<ELEMENT_DIM+1; i++)
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
            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteVtkResultsToFile()
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
    std::vector<std::string> cell_data_names;
    //We assume that the first cell is representative of all cells
    num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    cell_data_names = this->Begin()->GetCellData()->GetKeys();

    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cellwise_data_var(num_points);
        cellwise_data.push_back(cellwise_data_var);
    }

    if (mOutputMeshInVtk)
    {
        VtkMeshWriter<ELEMENT_DIM,SPACE_DIM> mesh_writer(this->mDirPath, "mesh_"+time.str(), false);
        mesh_writer.WriteFilesUsingMesh(rGetMesh());
    }

    if (mWriteVtkAsPoints)
    {
        VtkMeshWriter<SPACE_DIM,SPACE_DIM> cells_writer(this->mDirPath, "results_"+time.str(), false);

        // Loop over cells
        for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->Begin();
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
                cell_types[node_index] = cell_iter->GetCellProliferativeType()->GetColour();
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
                cellwise_data[var][node_index] = cell_iter->GetCellData()->GetItem(cell_data_names[var]);
            }
        }

        if (this->mOutputCellProliferativeTypes)
        {
            cells_writer.AddPointData("Cell types", cell_types);
        }
        if (this->mOutputCellAncestors)
        {
            cells_writer.AddPointData("Ancestors", cell_ancestors);
        }
        if (this->mOutputCellMutationStates)
        {
            cells_writer.AddPointData("Mutation states", cell_mutation_states);
        }
        if (this->mOutputCellAges)
        {
            cells_writer.AddPointData("Ages", cell_ages);
        }
        if (this->mOutputCellCyclePhases)
        {
            cells_writer.AddPointData("Cycle phases", cell_cycle_phases);
        }
        if (num_cell_data_items > 0)
        {
            for (unsigned var=0; var<cellwise_data.size(); var++)
            {
                cells_writer.AddPointData(cell_data_names[var], cellwise_data[var]);
            }
        }

        {
            // Make a copy of the nodes in a disposable mesh for writing
            std::vector<Node<SPACE_DIM>* > nodes;
            for (unsigned index=0; index<this->mrMesh.GetNumNodes(); index++)
            {
                Node<SPACE_DIM>* p_node = this->mrMesh.GetNode(index);
                nodes.push_back(p_node);
            }

            NodesOnlyMesh<SPACE_DIM> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);	// Arbitrary cut off as connectivity not used.
            cells_writer.WriteFilesUsingMesh(mesh);
        }

        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
    else if (mpVoronoiTessellation != NULL)
    {
        VertexMeshWriter<ELEMENT_DIM,SPACE_DIM> mesh_writer(this->mDirPath, "results", false);
        std::vector<double> cell_volumes(num_points);

        // Loop over elements of mpVoronoiTessellation
        for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
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
                cell_types[elem_index] = p_cell->GetCellProliferativeType()->GetColour();
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
                cellwise_data[var][elem_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
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
                mesh_writer.AddCellData(cell_data_names[var], cellwise_data[var]);
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteVoronoiResultsToFile()
{
    assert(SPACE_DIM==2 || SPACE_DIM==3);
    assert(mpVoronoiTessellation != NULL);
    // Write time to file
    *mpVoronoiFile << SimulationTime::Instance()->GetTime() << " ";

    // Loop over elements of mpVoronoiTessellation
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
         elem_iter != mpVoronoiTessellation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in mpVoronoiTessellation
        unsigned elem_index = elem_iter->GetIndex();

        // Get the index of the corresponding node in mrMesh
        unsigned node_index = mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        // Write node index and location to file
        *mpVoronoiFile << node_index << " ";
        c_vector<double, SPACE_DIM> node_location = this->GetNode(node_index)->rGetLocation();
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *mpVoronoiFile << node_location[i] << " ";
        }

        double cell_volume = mpVoronoiTessellation->GetVolumeOfElement(elem_index);
        double cell_surface_area = mpVoronoiTessellation->GetSurfaceAreaOfElement(elem_index);
        *mpVoronoiFile << cell_volume << " " << cell_surface_area << " ";
    }
    *mpVoronoiFile << "\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteCellPopulationVolumeResultsToFile()
{
    assert(SPACE_DIM==2 || SPACE_DIM==3);

    // Write time to file
    *mpCellPopulationVolumesFile << SimulationTime::Instance()->GetTime() << " ";
    assert (this->mpVoronoiTessellation != NULL);

    // Don't use the Voronoi tessellation to calculate the total area of the mesh because it gives huge areas for boundary cells
    double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).GetVolume();
    double apoptotic_area = 0.0;

    // Loop over elements of mpVoronoiTessellation
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetVolumeOfCell(CellPtr pCell)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteCellVolumeResultsToFile()
{
    assert (mpVoronoiTessellation != NULL);
    assert(SPACE_DIM==2 || SPACE_DIM==3);

    // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    ///\todo It would be simpler to merely iterate over the cell population here

    // Loop over elements of mpVoronoiTessellation
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
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
            c_vector<double, SPACE_DIM> node_location = this->GetNode(node_index)->rGetLocation();
            for (unsigned i=0; i<SPACE_DIM; i++)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetWriteVtkAsPoints(bool writeVtkAsPoints)
{
    mWriteVtkAsPoints = writeVtkAsPoints;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetWriteVtkAsPoints()
{
    return mWriteVtkAsPoints;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetOutputMeshInVtk(bool outputMeshInVtk)
{
    mOutputMeshInVtk = outputMeshInVtk;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetOutputMeshInVtk()
{
    return mOutputMeshInVtk;
}

//////////////////////////////////////////////////////////////////////////////
//                          Spring iterator class                           //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::GetNodeA()
{
    return mEdgeIter.GetNodeA();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::GetNodeB()
{
    return mEdgeIter.GetNodeB();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPtr MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::GetCellA()
{
    assert((*this) != mrCellPopulation.SpringsEnd());
    return mrCellPopulation.GetCellUsingLocationIndex(mEdgeIter.GetNodeA()->GetIndex());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPtr MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::GetCellB()
{
    assert((*this) != mrCellPopulation.SpringsEnd());
    return mrCellPopulation.GetCellUsingLocationIndex(mEdgeIter.GetNodeB()->GetIndex());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::operator!=(const MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator& rOther)
{
    return (mEdgeIter != rOther.mEdgeIter);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator& MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::operator++()
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::SpringIterator(
            MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
            typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::EdgeIterator edgeIter)
    : mrCellPopulation(rCellPopulation),
      mEdgeIter(edgeIter)
{
    if (mEdgeIter!=static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation.mrMesh))->EdgesEnd())
    {
        bool a_is_ghost = mrCellPopulation.IsGhostNode(mEdgeIter.GetNodeA()->GetIndex());
        bool b_is_ghost = mrCellPopulation.IsGhostNode(mEdgeIter.GetNodeB()->GetIndex());

        if (a_is_ghost || b_is_ghost)
        {
            ++(*this);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringsBegin()
{
    return SpringIterator(*this, static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).EdgesBegin());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringsEnd()
{
    return SpringIterator(*this, static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).EdgesEnd());
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
 * Can't tessellate 2d meshes in 3d space yet.
 */
template<>
void MeshBasedCellPopulation<2,3>::CreateVoronoiTessellation()
{
    // We don't allow tessellation yet.
    NEVER_REACHED;
}

/**
 * The cylindrical mesh is only defined in 2D, hence there is
 * a separate definition for this method in 3D, which doesn't have the capability
 * of dealing with periodic boundaries in 3D. This is \todo #1374.
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
void MeshBasedCellPopulation<1, 1>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}

/**
 * The VoronoiTessellation class is only defined in 2D or 3D, hence there
 * are two definitions to this method (one templated and one not).
 */
template<>
void MeshBasedCellPopulation<1, 2>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}

/**
 * The VoronoiTessellation class is only defined in 2D or 3D, hence there
 * are two definitions to this method (one templated and one not).
 */
template<>
void MeshBasedCellPopulation<1, 3>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM,SPACE_DIM>* MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetVoronoiTessellation()
{
    assert(mpVoronoiTessellation!=NULL);
    return mpVoronoiTessellation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetVolumeOfVoronoiElement(unsigned index)
{
    unsigned element_index = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
    double volume = mpVoronoiTessellation->GetVolumeOfElement(element_index);
    return volume;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetSurfaceAreaOfVoronoiElement(unsigned index)
{
    unsigned element_index = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(index);
    double surface_area = mpVoronoiTessellation->GetSurfaceAreaOfElement(element_index);
    return surface_area;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetVoronoiEdgeLength(unsigned index1, unsigned index2)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::CheckCellPointers()
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
    for (std::set<std::pair<CellPtr,CellPtr> >::iterator it1 = this->mMarkedSprings.begin();
         it1 != this->mMarkedSprings.end();
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetAreaBasedDampingConstantParameter()
{
    return mAreaBasedDampingConstantParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetAreaBasedDampingConstantParameter(double areaBasedDampingConstantParameter)
{
    assert(areaBasedDampingConstantParameter >= 0.0);
    mAreaBasedDampingConstantParameter = areaBasedDampingConstantParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetOutputVoronoiData()
{
    return mOutputVoronoiData;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetOutputVoronoiData(bool outputVoronoiData)
{
    mOutputVoronoiData = outputVoronoiData;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetOutputCellPopulationVolumes()
{
    return mOutputCellPopulationVolumes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetOutputCellPopulationVolumes(bool outputCellPopulationVolumes)
{
    mOutputCellPopulationVolumes = outputCellPopulationVolumes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >& MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::rGetNodePairs()
{
    //mNodePairs.Clear();
    NEVER_REACHED;
    return mNodePairs;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<UseAreaBasedDampingConstant>" << mUseAreaBasedDampingConstant << "</UseAreaBasedDampingConstant>\n";
    *rParamsFile << "\t\t<AreaBasedDampingConstantParameter>" <<  mAreaBasedDampingConstantParameter << "</AreaBasedDampingConstantParameter>\n";
    *rParamsFile << "\t\t<OutputVoronoiData>" <<  mOutputVoronoiData << "</OutputVoronoiData>\n";
    *rParamsFile << "\t\t<OutputCellPopulationVolumes>" << mOutputCellPopulationVolumes << "</OutputCellPopulationVolumes>\n";
    *rParamsFile << "\t\t<WriteVtkAsPoints>" << mWriteVtkAsPoints << "</WriteVtkAsPoints>\n";
    *rParamsFile << "\t\t<OutputMeshInVtk>" << mOutputMeshInVtk << "</OutputMeshInVtk>\n";
    *rParamsFile << "\t\t<HasVariableRestLength>" << mHasVariableRestLength << "</HasVariableRestLength>\n";

    // Call method on direct parent class
    AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = this->mrMesh.GetWidth(rDimension);
    return width;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    // Get pointer to this node
    Node<SPACE_DIM>* p_node = this->mrMesh.GetNode(index);

    // Loop over containing elements
    std::set<unsigned> neighbouring_node_indices;
    for (typename Node<SPACE_DIM>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
         elem_iter != p_node->ContainingElementsEnd();
         ++elem_iter)
    {
        // Get pointer to this containing element
        Element<ELEMENT_DIM,SPACE_DIM>* p_element = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).GetElement(*elem_iter);

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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::CalculateRestLengths()
{
    mSpringRestLengths.clear();

    // Iterate over all springs and add calculate separation of adjacent  node pairs
    for (SpringIterator spring_iterator = SpringsBegin();
         spring_iterator != SpringsEnd();
         ++spring_iterator)
    {
        // Note that nodeA_global_index is always less than nodeB_global_index
        Node<SPACE_DIM>* p_nodeA = spring_iterator.GetNodeA();
        Node<SPACE_DIM>* p_nodeB = spring_iterator.GetNodeB();

        unsigned nodeA_global_index = p_nodeA->GetIndex();
        unsigned nodeB_global_index = p_nodeB->GetIndex();

        // Calculate the distance between nodes
        c_vector<double, SPACE_DIM> node_a_location = p_nodeA->rGetLocation();
        c_vector<double, SPACE_DIM> node_b_location = p_nodeB->rGetLocation();

       double separation = norm_2(rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location));

       std::pair<unsigned,unsigned> node_pair;
       if (nodeA_global_index<nodeB_global_index)
       {
            node_pair.first = nodeA_global_index;
            node_pair.second  = nodeB_global_index;
       }
       else
       {
            node_pair.first = nodeB_global_index;
            node_pair.second  = nodeA_global_index;
       }

       mSpringRestLengths[node_pair]= separation;
    }
    mHasVariableRestLength = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetRestLength(unsigned indexA, unsigned indexB)
{
    if (mHasVariableRestLength)
    {
        if(indexA>indexB)
        {
            unsigned temp = indexA;
            indexA = indexB;
            indexB = temp;
        }

        std::pair<unsigned,unsigned> node_pair (indexA, indexB);

        std::map<std::pair<unsigned,unsigned>, double>::const_iterator  iter = mSpringRestLengths.find(node_pair);


        if (iter != mSpringRestLengths.end() )
        {
            // Return the stored rest length.
            return iter->second;
        }
        else
        {
            EXCEPTION("Tried to get a rest length of an edge that doesn't exist. You can only use variable rest lengths if SetUpdateCellPopulationRule is set on the simulation.");
        }
    }
    else
    {
        return 1.0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<unsigned,unsigned>  MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::CreateOrderedPair(unsigned index1, unsigned index2)
{
    assert(index1 != index2);

    std::pair<unsigned,unsigned> ordered_pair;

    if (index1<index2)
    {
        ordered_pair.first = index1;
        ordered_pair.second  = index2;
    }
    else
    {
        ordered_pair.first = index2;
        ordered_pair.second  = index1;
    }
    return ordered_pair;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MeshBasedCellPopulation<1,1>;
template class MeshBasedCellPopulation<1,2>;
template class MeshBasedCellPopulation<2,2>;
template class MeshBasedCellPopulation<1,3>;
template class MeshBasedCellPopulation<2,3>;
template class MeshBasedCellPopulation<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MeshBasedCellPopulation)
