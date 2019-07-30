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

#include "MeshBasedCellPopulation.hpp"
#include "VtkMeshWriter.hpp"
#include "CellBasedEventHandler.hpp"
#include "Cylindrical2dMesh.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellId.hpp"
#include "CellVolumesWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellPopulationAreaWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::MeshBasedCellPopulation(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh,
                                      bool validate)
    : AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>(rMesh, rCells, locationIndices),
      mpVoronoiTessellation(nullptr),
      mDeleteMesh(deleteMesh),
      mUseAreaBasedDampingConstant(false),
      mAreaBasedDampingConstantParameter(0.1),
      mWriteVtkAsPoints(false),
      mOutputMeshInVtk(false),
      mHasVariableRestLength(false)
{
    mpMutableMesh = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>* >(&(this->mrMesh));

    assert(this->mCells.size() <= this->mrMesh.GetNumNodes());

    if (validate)
    {
        Validate();
    }

    // Initialise the applied force at each node to zero
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::MeshBasedCellPopulation(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    : AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>(rMesh)
{
    mpMutableMesh = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>* >(&(this->mrMesh));
    mpVoronoiTessellation = nullptr;
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
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE
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

    /*
     * The next code block computes the area-dependent damping constant as given by equation
     * (5) in the following reference: van Leeuwen et al. 2009. An integrative computational model
     * for intestinal tissue renewal. Cell Prolif. 42(5):617-636. doi:10.1111/j.1365-2184.2009.00627.x
     */
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
        assert(SPACE_DIM == 2); // LCOV_EXCL_LINE

        double rest_length = 1.0;
        double d0 = mAreaBasedDampingConstantParameter;

        /**
         * Compute the parameter d1 such that d0+A*d1=1, where A is the equilibrium area
         * of a cell (this is equal to sqrt(3.0)/4, which is a third of the area of a regular
         * hexagon of edge length 1)
         */
        double d1 = 2.0*(1.0 - d0)/(sqrt(3.0)*rest_length*rest_length);

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
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Node " << i << " does not appear to have a cell associated with it");
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
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetTetrahedralMeshForPdeModifier()
{
    return mpMutableMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<CellPtr>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         )
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
        }
        else
        {
            ++it;
        }
    }

    return num_removed;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::Update(bool hasHadBirthsOrDeaths)
{
    ///\todo check if there is a more efficient way of keeping track of node velocity information (#2404)
    bool output_node_velocities = (this-> template HasWriter<NodeVelocityWriter>());

    /**
     * If node radii are set, then we must keep a record of these, since they will be cleared during
     * the remeshing process. We then restore these attributes to the nodes after calling ReMesh().
     *
     * At present, we check whether node radii are set by interrogating the radius of the first node
     * in the mesh and asking if it is strictly greater than zero (the default value, as set in the
     * NodeAttributes constructor). Hence, we assume that either ALL node radii are set, or NONE are.
     *
     * \todo There may be a better way of checking if node radii are set (#2694)
     */
    std::map<unsigned, double> old_node_radius_map;
    old_node_radius_map.clear();
    if (this->mrMesh.GetNodeIteratorBegin()->HasNodeAttributes())
    {
        if (this->mrMesh.GetNodeIteratorBegin()->GetRadius() > 0.0)
        {
            for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
                 node_iter != this->mrMesh.GetNodeIteratorEnd();
                 ++node_iter)
            {
                unsigned node_index = node_iter->GetIndex();
                old_node_radius_map[node_index] = node_iter->GetRadius();
            }
        }
    }

    std::map<unsigned, c_vector<double, SPACE_DIM> > old_node_applied_force_map;
    old_node_applied_force_map.clear();
    if (output_node_velocities)
    {
        /*
         * If outputting node velocities, we must keep a record of the applied force at each
         * node, since this will be cleared during the remeshing process. We then restore
         * these attributes to the nodes after calling ReMesh().
         */
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
             node_iter != this->mrMesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            old_node_applied_force_map[node_index] = node_iter->rGetAppliedForce();
        }
    }

    NodeMap node_map(this->mrMesh.GetNumAllNodes());

    // We must use a static_cast to call ReMesh() as this method is not defined in parent mesh classes
    static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((this->mrMesh)).ReMesh(node_map);

    if (!node_map.IsIdentityMap())
    {
        UpdateGhostNodesAfterReMesh(node_map);

        // Update the mappings between cells and location indices
        std::map<Cell*, unsigned> old_cell_location_map = this->mCellLocationMap;

        // Remove any dead pointers from the maps (needed to avoid archiving errors)
        this->mLocationCellMap.clear();
        this->mCellLocationMap.clear();

        for (std::list<CellPtr>::iterator it = this->mCells.begin(); it != this->mCells.end(); ++it)
        {
            unsigned old_node_index = old_cell_location_map[(*it).get()];

            // This shouldn't ever happen, as the cell vector only contains living cells
            assert(!node_map.IsDeleted(old_node_index));

            unsigned new_node_index = node_map.GetNewIndex(old_node_index);
            this->SetCellUsingLocationIndex(new_node_index,*it);

            if (old_node_radius_map[old_node_index] > 0.0)
            {
                this->GetNode(new_node_index)->SetRadius(old_node_radius_map[old_node_index]);
            }
            if (output_node_velocities)
            {
                this->GetNode(new_node_index)->AddAppliedForceContribution(old_node_applied_force_map[old_node_index]);
            }
        }

        this->Validate();
    }
    else
    {
        if (old_node_radius_map[this->mCellLocationMap[(*(this->mCells.begin())).get()]] > 0.0)
        {
            for (std::list<CellPtr>::iterator it = this->mCells.begin(); it != this->mCells.end(); ++it)
            {
                unsigned node_index = this->mCellLocationMap[(*it).get()];
                this->GetNode(node_index)->SetRadius(old_node_radius_map[node_index]);
            }
        }
        if (output_node_velocities)
        {
            for (std::list<CellPtr>::iterator it = this->mCells.begin(); it != this->mCells.end(); ++it)
            {
                unsigned node_index = this->mCellLocationMap[(*it).get()];
                this->GetNode(node_index)->AddAppliedForceContribution(old_node_applied_force_map[node_index]);
            }
        }
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
        if (mUseAreaBasedDampingConstant ||
            this-> template HasWriter<VoronoiDataWriter>() ||
            this-> template HasWriter<CellPopulationAreaWriter>() ||
            this-> template HasWriter<CellVolumesWriter>())
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
    assert(ELEMENT_DIM == 2); // LCOV_EXCL_LINE

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
                                    p_neighbour_cell->GetSrnModel()->CreateSrnModel(),
                                    false,
                                    daughter_property_collection));

        // Add new cell to cell population
        this->mCells.push_back(p_new_cell);
        this->AddCellUsingLocationIndex(new_node_index,p_new_cell);

        // Update rest lengths

        // Remove old node pair // note node_a_index < node_b_index
        std::pair<unsigned,unsigned> node_pair = this->CreateOrderedPair(node_a_index, node_b_index);
        double old_rest_length  = mSpringRestLengths[node_pair];

        std::map<std::pair<unsigned,unsigned>, double>::iterator  iter = mSpringRestLengths.find(node_pair);
        mSpringRestLengths.erase(iter);

        // Add new pairs
        node_pair = this->CreateOrderedPair(node_a_index, new_node_index);
        mSpringRestLengths[node_pair] = 0.5*old_rest_length;

        node_pair = this->CreateOrderedPair(node_b_index, new_node_index);
        mSpringRestLengths[node_pair] = 0.5*old_rest_length;

        // If necessary add other new spring rest lengths
        for (unsigned pair_index=3; pair_index<5; pair_index++)
        {
            unsigned other_node_index = new_nodes[index][pair_index];

            if (other_node_index != UNSIGNED_UNSET)
            {
                node_pair = this->CreateOrderedPair(other_node_index, new_node_index);
                double new_rest_length = rGetMesh().GetDistanceBetweenNodes(new_node_index, other_node_index);
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
CellPtr MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    assert(pNewCell);
    assert(pParentCell);

    // Add new cell to population
    CellPtr p_created_cell = AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AddCell(pNewCell, pParentCell);
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
void MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellPopulationElementWriter>())
        {
            this-> template AddPopulationWriter<CellPopulationElementWriter>();
        }
    }

    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::OpenWritersFiles(rOutputFileHandler);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteResultsToFiles(const std::string& rDirectory)
{
    if (SimulationTime::Instance()->GetTimeStepsElapsed() == 0 && this->mpVoronoiTessellation == nullptr)
    {
        TessellateIfNeeded(); // Update isn't run on time-step zero
    }

    AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteResultsToFiles(rDirectory);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM> > pPopulationCountWriter)
{
    pPopulationCountWriter->Visit(this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    // Store the present time as a string
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    // Store the number of cells for which to output data to VTK
    unsigned num_cells_from_mesh = GetNumNodes();
    if (!mWriteVtkAsPoints && (mpVoronoiTessellation != nullptr))
    {
        num_cells_from_mesh = mpVoronoiTessellation->GetNumElements();
    }

    // When outputting any CellData, we assume that the first cell is representative of all cells
    unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

    std::vector<std::vector<double> > cell_data;
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cell_data_var(num_cells_from_mesh);
        cell_data.push_back(cell_data_var);
    }

    if (mOutputMeshInVtk)
    {
        // Create mesh writer for VTK output
        VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(rDirectory, "mesh_"+time.str(), false);
        mesh_writer.WriteFilesUsingMesh(rGetMesh());
    }

    if (mWriteVtkAsPoints)
    {
        // Create mesh writer for VTK output
        VtkMeshWriter<SPACE_DIM, SPACE_DIM> cells_writer(rDirectory, "results_"+time.str(), false);

        // Iterate over any cell writers that are present
        unsigned num_cells = this->GetNumAllCells();
        for (typename std::vector<boost::shared_ptr<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
             cell_writer_iter != this->mCellWriters.end();
             ++cell_writer_iter)
        {
            // Create vector to store VTK cell data
            std::vector<double> vtk_cell_data(num_cells);

            // Loop over cells
            for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->Begin();
                 cell_iter != this->End();
                 ++cell_iter)
            {
                // Get the node index corresponding to this cell
                unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

                // Populate the vector of VTK cell data
                vtk_cell_data[node_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(*cell_iter, this);
            }

            cells_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
        }

        // Loop over cells
        for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->Begin();
             cell_iter != this->End();
             ++cell_iter)
        {
            // Get the node index corresponding to this cell
            unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cell_data[var][node_index] = cell_iter->GetCellData()->GetItem(cell_data_names[var]);
            }
        }
        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cells_writer.AddPointData(cell_data_names[var], cell_data[var]);
        }

        // Make a copy of the nodes in a disposable mesh for writing
        {
            std::vector<Node<SPACE_DIM>* > nodes;
            for (unsigned index=0; index<this->mrMesh.GetNumNodes(); index++)
            {
                Node<SPACE_DIM>* p_node = this->mrMesh.GetNode(index);
                nodes.push_back(p_node);
            }

            NodesOnlyMesh<SPACE_DIM> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5); // Arbitrary cut off as connectivity not used.
            cells_writer.WriteFilesUsingMesh(mesh);
        }

        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
    else if (mpVoronoiTessellation != nullptr)
    {
        // Create mesh writer for VTK output
        VertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(rDirectory, "results", false);
        std::vector<double> cell_volumes(num_cells_from_mesh);

        // Iterate over any cell writers that are present
        unsigned num_cells = this->GetNumAllCells();
        for (typename std::vector<boost::shared_ptr<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
             cell_writer_iter != this->mCellWriters.end();
             ++cell_writer_iter)
        {
            // Create vector to store VTK cell data
            std::vector<double> vtk_cell_data(num_cells);

            // Loop over elements of mpVoronoiTessellation
            for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
                 elem_iter != mpVoronoiTessellation->GetElementIteratorEnd();
                 ++elem_iter)
            {
                // Get index of this element in mpVoronoiTessellation
                unsigned elem_index = elem_iter->GetIndex();

                // Get the cell corresponding to this element, via the index of the corresponding node in mrMesh
                unsigned node_index = mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);
                CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

                // Populate the vector of VTK cell data
                vtk_cell_data[elem_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
            }

            mesh_writer.AddCellData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
        }

        // Loop over elements of mpVoronoiTessellation
        for (typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_iter = mpVoronoiTessellation->GetElementIteratorBegin();
             elem_iter != mpVoronoiTessellation->GetElementIteratorEnd();
             ++elem_iter)
        {
            // Get index of this element in mpVoronoiTessellation
            unsigned elem_index = elem_iter->GetIndex();

            // Get the cell corresponding to this element, via the index of the corresponding node in mrMesh
            unsigned node_index = mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);
            CellPtr p_cell = this->GetCellUsingLocationIndex(node_index);

            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cell_data[var][elem_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
            }
        }

        for (unsigned var=0; var<cell_data.size(); var++)
        {
            mesh_writer.AddCellData(cell_data_names[var], cell_data[var]);
        }

        mesh_writer.WriteVtkUsingMesh(*mpVoronoiTessellation, time.str());
        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetVolumeOfCell(CellPtr pCell)
{
    double cell_volume = 0;

    if (ELEMENT_DIM == SPACE_DIM)
    {
        // Ensure that the Voronoi tessellation exists
        if (mpVoronoiTessellation == nullptr)
        {
            CreateVoronoiTessellation();
        }

        // Get the node index corresponding to this cell
        unsigned node_index = this->GetLocationIndexUsingCell(pCell);

        // Try to get the element index of the Voronoi tessellation corresponding to this node index
        try
        {
            unsigned element_index = mpVoronoiTessellation->GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(node_index);

            // Get the cell's volume from the Voronoi tessellation
            cell_volume = mpVoronoiTessellation->GetVolumeOfElement(element_index);
        }
        catch (Exception&)
        {
            // If it doesn't exist this must be a boundary cell, so return infinite volume
            cell_volume = DBL_MAX;
        }
    }
    else if (SPACE_DIM==3 && ELEMENT_DIM==2)
    {
        unsigned node_index = this->GetLocationIndexUsingCell(pCell);

        Node<SPACE_DIM>* p_node = rGetMesh().GetNode(node_index);

        assert(!(p_node->rGetContainingElementIndices().empty()));

        for (typename Node<SPACE_DIM>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
             elem_iter != p_node->ContainingElementsEnd();
             ++elem_iter)
        {
            Element<ELEMENT_DIM,SPACE_DIM>* p_element = rGetMesh().GetElement(*elem_iter);

            c_matrix<double, SPACE_DIM, ELEMENT_DIM> jacob;
            double det;

            p_element->CalculateJacobian(jacob, det);

            cell_volume += fabs(p_element->GetVolume(det));
        }

        // This calculation adds a third of each element to the total area
        cell_volume /= 3.0;
    }
    else
    {
        // Not implemented for other dimensions
        NEVER_REACHED;
    }

    return cell_volume;
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
    if (bool(dynamic_cast<Cylindrical2dMesh*>(&(this->mrMesh))))
    {
        *pVizSetupFile << "MeshWidth\t" << this->GetWidth(0) << "\n";
    }
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
bool MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator::operator!=(const typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator& rOther)
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
    if (bool(dynamic_cast<Cylindrical2dMesh*>(&mrMesh)))
    {
        is_mesh_periodic = true;
        mpVoronoiTessellation = new Cylindrical2dVertexMesh(static_cast<Cylindrical2dMesh &>(this->mrMesh));
    }
    else
    {
        mpVoronoiTessellation = new VertexMesh<2, 2>(static_cast<MutableMesh<2, 2> &>((this->mrMesh)), is_mesh_periodic);
    }
}

/**
 * Can't tessellate 2d meshes in 3d space yet.
 */
// LCOV_EXCL_START
template<>
void MeshBasedCellPopulation<2,3>::CreateVoronoiTessellation()
{
    // We don't allow tessellation yet.
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

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
// LCOV_EXCL_START
template<>
void MeshBasedCellPopulation<1, 1>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}
// LCOV_EXCL_STOP


/**
 * The VoronoiTessellation class is only defined in 2D or 3D, hence there
 * are two definitions to this method (one templated and one not).
 */
// LCOV_EXCL_START
template<>
void MeshBasedCellPopulation<1, 2>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

/**
 * The VoronoiTessellation class is only defined in 2D or 3D, hence there
 * are two definitions to this method (one templated and one not).
 */
// LCOV_EXCL_START
template<>
void MeshBasedCellPopulation<1, 3>::CreateVoronoiTessellation()
{
    // No 1D Voronoi tessellation
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM,SPACE_DIM>* MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetVoronoiTessellation()
{
    assert(mpVoronoiTessellation!=nullptr);
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
    catch (Exception&)
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
// LCOV_EXCL_START //Debugging code.  Shouldn't fail under normal conditions
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
    UNUSED_OPT(res);
    assert(res);
// LCOV_EXCL_STOP

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

// LCOV_EXCL_START //Debugging code.  Shouldn't fail under normal conditions
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
// LCOV_EXCL_STOP
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

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >& MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::rGetNodePairs()
{
    //mNodePairs.Clear();
    NEVER_REACHED;
    return mNodePairs;
}
// LCOV_EXCL_STOP

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<UseAreaBasedDampingConstant>" << mUseAreaBasedDampingConstant << "</UseAreaBasedDampingConstant>\n";
    *rParamsFile << "\t\t<AreaBasedDampingConstantParameter>" <<  mAreaBasedDampingConstantParameter << "</AreaBasedDampingConstantParameter>\n";
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
        double separation = rGetMesh().GetDistanceBetweenNodes(nodeA_global_index, nodeB_global_index);

        // Order node indices
        std::pair<unsigned,unsigned> node_pair = this->CreateOrderedPair(nodeA_global_index, nodeB_global_index);

        mSpringRestLengths[node_pair] = separation;
    }
    mHasVariableRestLength = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::GetRestLength(unsigned indexA, unsigned indexB)
{
    if (mHasVariableRestLength)
    {
        std::pair<unsigned,unsigned> node_pair = this->CreateOrderedPair(indexA, indexB);
        std::map<std::pair<unsigned,unsigned>, double>::const_iterator iter = mSpringRestLengths.find(node_pair);

        if (iter != mSpringRestLengths.end())
        {
            // Return the stored rest length
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
void MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SetRestLength(unsigned indexA, unsigned indexB, double restLength)
{
    if (mHasVariableRestLength)
    {
        std::pair<unsigned,unsigned> node_pair = this->CreateOrderedPair(indexA, indexB);
        std::map<std::pair<unsigned,unsigned>, double>::iterator iter = mSpringRestLengths.find(node_pair);

        if (iter != mSpringRestLengths.end())
        {
            // modify the stored rest length
            iter->second = restLength;
        }
        else
        {
            EXCEPTION("Tried to set the rest length of an edge not in the mesh.");
        }
    }
    else
    {
        EXCEPTION("Tried to set a rest length in a simulation with fixed rest length. You can only use variable rest lengths if SetUpdateCellPopulationRule is set on the simulation.");
    }
}

// Explicit instantiation
template class MeshBasedCellPopulation<1,1>;
template class MeshBasedCellPopulation<1,2>;
template class MeshBasedCellPopulation<2,2>;
template class MeshBasedCellPopulation<1,3>;
template class MeshBasedCellPopulation<2,3>;
template class MeshBasedCellPopulation<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MeshBasedCellPopulation)
