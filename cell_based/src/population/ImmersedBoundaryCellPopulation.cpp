/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "ImmersedBoundaryCellPopulation.hpp"

#include "CellPopulationElementWriter.hpp"
#include "ImmersedBoundaryMeshWriter.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "Warnings.hpp"

template <unsigned DIM>
ImmersedBoundaryCellPopulation<DIM>::ImmersedBoundaryCellPopulation(ImmersedBoundaryMesh<DIM, DIM>& rMesh,
                                                                    std::vector<CellPtr>& rCells,
                                                                    bool deleteMesh,
                                                                    bool validate,
                                                                    const std::vector<unsigned> locationIndices)
        : AbstractOffLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices),
          mDeleteMesh(deleteMesh),
          mIntrinsicSpacing(0.01),
          mPopulationHasActiveSources(false),
          mReMeshFrequency(UINT_MAX)
{
    mpImmersedBoundaryMesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&(this->mrMesh));
    mpVertexBasedDivisionRule.reset(new ShortAxisVertexBasedDivisionRule<DIM>());

    // If no location indices are specified, associate with elements from the mesh (assumed to be sequentially ordered).
    std::list<CellPtr>::iterator it = this->mCells.begin();
    for (unsigned i = 0; it != this->mCells.end(); ++it, ++i)
    {
        unsigned index = locationIndices.empty() ? i : locationIndices[i]; // assume that the ordering matches
        AbstractCellPopulation<DIM, DIM>::AddCellUsingLocationIndex(index, *it);
    }

    // Check each element has only one cell attached
    if (validate)
    {
        Validate();
    }

    mInteractionDistance = 0.05 * CalculateIntrinsicCellSize();

    // Set the mesh division spacing distance
    rMesh.SetElementDivisionSpacing(0.25 * mInteractionDistance);
}

template <unsigned DIM>
ImmersedBoundaryCellPopulation<DIM>::ImmersedBoundaryCellPopulation(ImmersedBoundaryMesh<DIM, DIM>& rMesh)
        : AbstractOffLatticeCellPopulation<DIM>(rMesh),
          mDeleteMesh(true)
{
    mpImmersedBoundaryMesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&(this->mrMesh));
}

template <unsigned DIM>
ImmersedBoundaryCellPopulation<DIM>::~ImmersedBoundaryCellPopulation()
{
    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::CalculateIntrinsicCellSize()
{
    double average_intrinsic_size = 0.0;

    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_iter = mpImmersedBoundaryMesh->GetElementIteratorBegin();
         elem_iter != mpImmersedBoundaryMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        average_intrinsic_size += mpImmersedBoundaryMesh->GetVolumeOfElement(elem_iter->GetIndex());
    }

    return sqrt(average_intrinsic_size / mpImmersedBoundaryMesh->GetNumElements());
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    return 0.0;
}

template <unsigned DIM>
ImmersedBoundaryMesh<DIM, DIM>& ImmersedBoundaryCellPopulation<DIM>::rGetMesh()
{
    return *mpImmersedBoundaryMesh;
}

template <unsigned DIM>
const ImmersedBoundaryMesh<DIM, DIM>& ImmersedBoundaryCellPopulation<DIM>::rGetMesh() const
{
    return *mpImmersedBoundaryMesh;
}

template <unsigned DIM>
ImmersedBoundaryElement<DIM, DIM>* ImmersedBoundaryCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mpImmersedBoundaryMesh->GetElement(elementIndex);
}

template <unsigned DIM>
ImmersedBoundaryElement<DIM - 1, DIM>* ImmersedBoundaryCellPopulation<DIM>::GetLamina(unsigned laminaIndex)
{
    return mpImmersedBoundaryMesh->GetLamina(laminaIndex);
}

template <unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template <unsigned DIM>
c_vector<double, DIM> ImmersedBoundaryCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mpImmersedBoundaryMesh->GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

template <unsigned DIM>
Node<DIM>* ImmersedBoundaryCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetInteractionDistance(double new_distance)
{
    assert(new_distance >= 0.0);
    mInteractionDistance = new_distance;
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetInteractionDistance()
{
    return mInteractionDistance;
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetReMeshFrequency(unsigned newFrequency)
{
    mReMeshFrequency = newFrequency;
}

template <unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::GetReMeshFrequency()
{
    return mReMeshFrequency;
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetIntrinsicSpacing()
{
    return mIntrinsicSpacing;
}

//\todo: implement this method. Decide what "neighbouring" should be for IB cells
template <unsigned DIM>
std::set<unsigned> ImmersedBoundaryCellPopulation<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    // The set to return
    std::set<unsigned> neighbouring_indices;
    return neighbouring_indices;
}

template <unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    std::cout << "Trying to add new node from within cell population\n";
    return 0; //mpImmersedBoundaryMesh->AddNode(pNewNode);
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mpImmersedBoundaryMesh->SetNode(nodeIndex, rNewLocation);
}

template <unsigned DIM>
ImmersedBoundaryElement<DIM, DIM>* ImmersedBoundaryCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mpImmersedBoundaryMesh->GetElement(this->GetLocationIndexUsingCell(pCell));
}

template <unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::GetNumElements()
{
    return mpImmersedBoundaryMesh->GetNumElements();
}

template <unsigned DIM>
CellPtr ImmersedBoundaryCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    // Get the element associated with this cell
    ImmersedBoundaryElement<DIM, DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Get the orientation of division
    c_vector<double, DIM> division_vector = unit_vector<double>(2, 0);

    // Divide the element
    unsigned new_elem_idx = mpImmersedBoundaryMesh->DivideElementAlongGivenAxis(p_element, division_vector, true);
    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(new_elem_idx, p_created_cell);
    this->mCellLocationMap[p_created_cell.get()] = new_elem_idx;

    return p_created_cell;
}

template <unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    //\todo: may need to implement this

    /*for (std::list<CellPtr>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         )
    {
        if ((*it)->IsDead())
        {
            // Count the cell as dead
            num_removed++;

            // Remove the element from the mesh if it is not deleted yet
            ///\todo (#2489) this should cause an error - we should fix this!
            if (!(this->GetElement(this->GetLocationIndexUsingCell((*it)))->IsDeleted()))
            {
                // This warning relies on the fact that there is only one other possibility for
                // vertex elements to be marked as deleted: a T2 swap
                WARN_ONCE_ONLY("A Cell is removed without performing a T2 swap. This could leave a void in the mesh.");
                mpImmersedBoundaryMesh->DeleteElementPriorToReMesh(this->GetLocationIndexUsingCell((*it)));
            }

            // Delete the cell
            it = this->mCells.erase(it);
        }
        else
        {
            ++it;
        }
    }*/
    return num_removed;
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::UpdateNodeLocations(double dt)
{
    // Helper variables, pre-declared for efficiency
    unsigned num_grid_pts_x = this->rGetMesh().GetNumGridPtsX();
    unsigned num_grid_pts_y = this->rGetMesh().GetNumGridPtsY();

    double characteristic_spacing = this->rGetMesh().GetCharacteristicNodeSpacing();
    double grid_spacing_x = 1.0 / (double)num_grid_pts_x;
    double grid_spacing_y = 1.0 / (double)num_grid_pts_y;

    unsigned first_idx_x;
    unsigned first_idx_y;

    std::vector<unsigned> x_indices(4);
    std::vector<unsigned> y_indices(4);

    std::vector<double> x_deltas(4);
    std::vector<double> y_deltas(4);

    double delta;

    c_vector<double, DIM> node_location;
    c_vector<double, DIM> displacement;

    // Get references to the fluid velocity grid
    const multi_array<double, 3>& vel_grids = this->rGetMesh().rGet2dVelocityGrids();

    // Iterate over all nodes
    for (typename ImmersedBoundaryMesh<DIM, DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin(false);
         node_iter != this->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // Get location of current node
        node_location = node_iter->rGetLocation();

        // Get first grid index in each dimension, taking account of possible wrap-around
        first_idx_x = unsigned(floor(node_location[0] / grid_spacing_x)) + num_grid_pts_x - 1;
        first_idx_y = unsigned(floor(node_location[1] / grid_spacing_y)) + num_grid_pts_y - 1;

        // Calculate all four indices and deltas in each dimension
        for (unsigned i = 0; i < 4; i++)
        {
            x_indices[i] = (first_idx_x + i) % num_grid_pts_x;
            y_indices[i] = (first_idx_y + i) % num_grid_pts_y;

            x_deltas[i] = Delta1D(fabs(x_indices[i] * grid_spacing_x - node_location[0]), grid_spacing_x);
            y_deltas[i] = Delta1D(fabs(y_indices[i] * grid_spacing_x - node_location[1]), grid_spacing_y);
        }

        // Loop over the 4x4 grid which will influence the displacement of the current node
        for (unsigned x_idx = 0; x_idx < 4; x_idx++)
        {
            for (unsigned y_idx = 0; y_idx < 4; y_idx++)
            {
                // The applied velocity is weighted by the delta function
                delta = x_deltas[x_idx] * y_deltas[y_idx];
                displacement[0] += vel_grids[0][x_indices[x_idx]][y_indices[y_idx]] * delta;
                displacement[1] += vel_grids[1][x_indices[x_idx]][y_indices[y_idx]] * delta;
            }
        }

        // Normalise by timestep
        displacement *= dt;

        //If the displacement is too big, warn the user once and scale it back
        if (norm_2(displacement) > characteristic_spacing)
        {
            WARN_ONCE_ONLY("Nodes are moving more than half the CharacteristicNodeSpacing. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
            displacement *= characteristic_spacing / norm_2(displacement);
        }

        // Get new node location
        node_location += displacement;

        // Account for periodic boundary
        for (unsigned i = 0; i < DIM; i++)
        {
            node_location[i] = fmod(node_location[i] + 1.0, 1.0);
        }

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(node_location);

        // Move the node
        this->SetNode(node_iter->GetIndex(), new_point);
    }

    // If active sources, we need to update those location as well
    if (this->DoesPopulationHaveActiveSources())
    {
        std::vector<FluidSource<DIM>*>& r_element_sources = this->rGetMesh().rGetElementFluidSources();
        std::vector<FluidSource<DIM>*>& r_balance_sources = this->rGetMesh().rGetBalancingFluidSources();

        // Construct a vector of all sources combined
        std::vector<FluidSource<DIM>*> combined_sources;
        combined_sources.insert(combined_sources.end(), r_element_sources.begin(), r_element_sources.end());
        combined_sources.insert(combined_sources.end(), r_balance_sources.begin(), r_balance_sources.end());

        c_vector<double, DIM> source_location;

        // Iterate over all sources and update their locations
        for (unsigned source_idx = 0; source_idx < combined_sources.size(); source_idx++)
        {
            // Get location of current node
            source_location = combined_sources[source_idx]->rGetLocation();

            // Get first grid index in each dimension, taking account of possible wrap-around
            first_idx_x = unsigned(floor(source_location[0] / grid_spacing_x)) + num_grid_pts_x - 1;
            first_idx_y = unsigned(floor(source_location[1] / grid_spacing_y)) + num_grid_pts_y - 1;

            // Calculate all four indices and deltas in each dimension
            for (unsigned i = 0; i < 4; i++)
            {
                x_indices[i] = (first_idx_x + i) % num_grid_pts_x;
                y_indices[i] = (first_idx_y + i) % num_grid_pts_y;

                x_deltas[i] = Delta1D(fabs(x_indices[i] * grid_spacing_x - source_location[0]), grid_spacing_x);
                y_deltas[i] = Delta1D(fabs(y_indices[i] * grid_spacing_x - source_location[1]), grid_spacing_y);
            }

            // Loop over the 4x4 grid which will influence the displacement of the current node
            for (unsigned x_idx = 0; x_idx < 4; x_idx++)
            {
                for (unsigned y_idx = 0; y_idx < 4; y_idx++)
                {
                    // The applied velocity is weighted by the delta function
                    delta = x_deltas[x_idx] * y_deltas[y_idx];
                    displacement[0] += vel_grids[0][x_indices[x_idx]][y_indices[y_idx]] * delta;
                    displacement[1] += vel_grids[1][x_indices[x_idx]][y_indices[y_idx]] * delta;
                }
            }

            // Normalise by timestep
            displacement *= dt;

            //If the displacement is too big, warn the user once and scale it back
            if (norm_2(displacement) > characteristic_spacing)
            {
                WARN_ONCE_ONLY("Sources are moving more than half the CharacteristicNodeSpacing. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
                displacement *= characteristic_spacing / norm_2(displacement);
            }

            // Get new node location
            source_location += displacement;

            // Account for periodic boundary
            for (unsigned i = 0; i < DIM; i++)
            {
                source_location[i] = fmod(source_location[i] + 1.0, 1.0);
            }

            // Move the node
            combined_sources[source_idx]->rGetModifiableLocation() = source_location;
        }
    }

    // Finally, call ReMesh if required
    if (SimulationTime::Instance()->GetTimeStepsElapsed() % mReMeshFrequency == 0)
    {
        mpImmersedBoundaryMesh->ReMesh();
    }
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::Delta1D(double dist, double spacing)
{
    return (0.25 * (1.0 + cos(M_PI * dist / (2 * spacing))));
}

template <unsigned DIM>
bool ImmersedBoundaryCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();
    ;
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    // I don't think this is needed for IB
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::Validate()
{
    // Check each element has only one cell attached
    std::vector<unsigned> validated_element = std::vector<unsigned>(this->GetNumElements(), 0);
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);
        validated_element[elem_index]++;
    }

    for (unsigned i = 0; i < validated_element.size(); i++)
    {
        if (validated_element[i] == 0)
        {
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Element " << i << " does not appear to have a cell associated with it");
        }

        if (validated_element[i] > 1)
        {
            // This should never be reached as you can only set one cell per element index
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Element " << i << " appears to have " << validated_element[i] << " cells associated with it");
        }
    }
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::CheckForStepSizeException(unsigned nodeIndex, c_vector<double, DIM>& displacement, double dt)
{
    //\todo: deal with this method...
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    //pPopulationWriter->Visit(this);
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter)
{
    //pPopulationCountWriter->Visit(this);
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Get the vertex element index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get the cell's volume from the vertex mesh
    double cell_volume = mpImmersedBoundaryMesh->GetVolumeOfElement(elem_index);

    return cell_volume;
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    // Create mesh writer for VTK output
    ImmersedBoundaryMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results", false);

    // Calculated the cell overlap information, and get the number of cell parts needed for each element
    mesh_writer.CalculateCellOverlaps(*mpImmersedBoundaryMesh);
    std::vector<unsigned> num_cell_parts = mesh_writer.rGetNumCellParts();

    // Iterate over any cell writers that are present
    for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
         cell_writer_iter != this->mCellWriters.end();
         ++cell_writer_iter)
    {
        // Create vector to store VTK cell data
        std::vector<double> vtk_cell_data;

        // Iterate over immersed boundary elements
        for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_iter = mpImmersedBoundaryMesh->GetElementIteratorBegin();
             elem_iter != mpImmersedBoundaryMesh->GetElementIteratorEnd();
             ++elem_iter)
        {
            // Get index of this element in the vertex mesh
            unsigned elem_index = elem_iter->GetIndex();

            // Get the cell corresponding to this element
            CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
            assert(p_cell);

            // Populate the vector of VTK cell data.  We loop over the number of output cells as this takes into
            // account that some elements will be broken into pieces for visualisation
            for (unsigned cell_part = 0; cell_part < num_cell_parts[elem_index]; cell_part++)
            {
                vtk_cell_data.push_back((*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this));
            }
        }

        // Iterate over immersed boundary laminas (no associated cell) to ensure vtk_cell_data is the correct size
        for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryLaminaIterator lam_iter = mpImmersedBoundaryMesh->GetLaminaIteratorBegin();
             lam_iter != mpImmersedBoundaryMesh->GetLaminaIteratorEnd();
             ++lam_iter)
        {
            vtk_cell_data.push_back(-1.0);
        }

        mesh_writer.AddCellData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    }

    // When outputting any CellData, we assume that the first cell is representative of all cells
    unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

    std::vector<std::vector<double> > cell_data;
    for (unsigned var = 0; var < num_cell_data_items; var++)
    {
        std::vector<double> cell_data_var;
        cell_data.push_back(cell_data_var);
    }

    // Iterate over immersed boundary elements
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_iter = mpImmersedBoundaryMesh->GetElementIteratorBegin();
         elem_iter != mpImmersedBoundaryMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in the vertex mesh
        unsigned elem_index = elem_iter->GetIndex();

        // Get the cell corresponding to this element
        CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
        assert(p_cell);

        for (unsigned var = 0; var < num_cell_data_items; var++)
        {
            // Populate the vector of VTK cell data.  We loop over the number of output cells as this takes into
            // account that some elements will be broken into pieces for visualisation
            for (unsigned cell_part = 0; cell_part < num_cell_parts[elem_index]; cell_part++)
            {
                cell_data[var].push_back(p_cell->GetCellData()->GetItem(cell_data_names[var]));
            }
        }
    }

    // Iterate over immersed boundary laminas (no associated cell) to ensure cell_data is the correct size
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryLaminaIterator lam_iter = mpImmersedBoundaryMesh->GetLaminaIteratorBegin();
         lam_iter != mpImmersedBoundaryMesh->GetLaminaIteratorEnd();
         ++lam_iter)
    {
        for (unsigned var = 0; var < num_cell_data_items; var++)
        {
            cell_data[var].push_back(DOUBLE_UNSET);
        }
    }

    for (unsigned var = 0; var < num_cell_data_items; var++)
    {
        mesh_writer.AddCellData(cell_data_names[var], cell_data[var]);
    }

    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    mesh_writer.WriteVtkUsingMesh(*mpImmersedBoundaryMesh, time.str());

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this->template HasWriter<CellPopulationElementWriter>())
        {
            this->template AddPopulationWriter<CellPopulationElementWriter>();
        }
    }

    AbstractCellPopulation<DIM>::OpenWritersFiles(rOutputFileHandler);
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Add the division rule parameters
    *rParamsFile << "\t\t<VertexBasedDivisionRule>\n";
    mpVertexBasedDivisionRule->OutputCellVertexBasedDivisionRuleInfo(rParamsFile);
    *rParamsFile << "\t\t</VertexBasedDivisionRule>\n";

    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    double width = this->mrMesh.GetWidth(rDimension);
    return width;
}

///\todo: implement this.  May need to put method back in to mesh class
template <unsigned DIM>
std::set<unsigned> ImmersedBoundaryCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    std::set<unsigned> neighbouring_indices;
    return neighbouring_indices;
}

template <unsigned DIM>
TetrahedralMesh<DIM, DIM>* ImmersedBoundaryCellPopulation<DIM>::GetTetrahedralMeshForPdeModifier()
{
    //    // This method only works in 2D sequential
    //    assert(DIM == 2);
    //    assert(PetscTools::IsSequential());
    //
    //    unsigned num_vertex_nodes = mpMutableVertexMesh->GetNumNodes();
    //    unsigned num_vertex_elements = mpMutableVertexMesh->GetNumElements();
    //
    //    std::string mesh_file_name = "mesh";
    //
    //    // Get a unique temporary foldername
    //    std::stringstream pid;
    //    pid << getpid();
    //    OutputFileHandler output_file_handler("2D_temporary_tetrahedral_mesh_" + pid.str());
    //    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();
    //
    //    // Compute the number of nodes in the TetrahedralMesh
    //    unsigned num_tetrahedral_nodes = num_vertex_nodes + num_vertex_elements;
    //
    //    // Write node file
    //    out_stream p_node_file = output_file_handler.OpenOutputFile(mesh_file_name+".node");
    //    (*p_node_file) << std::scientific;
    //    (*p_node_file) << std::setprecision(20);
    //    (*p_node_file) << num_tetrahedral_nodes << "\t2\t0\t1" << std::endl;
    //
    //    // Begin by writing each node in the VertexMesh
    //    for (unsigned node_index=0; node_index<num_vertex_nodes; node_index++)
    //    {
    //        Node<DIM>* p_node = mpMutableVertexMesh->GetNode(node_index);
    //
    //        ///\todo will the nodes in mpMutableVertexMesh always have indices 0,1,2,...? (#2221)
    //        unsigned index = p_node->GetIndex();
    //        const c_vector<double, DIM>& r_location = p_node->rGetLocation();
    //        unsigned is_boundary_node = p_node->IsBoundaryNode() ? 1 : 0;
    //
    //        (*p_node_file) << index << "\t" << r_location[0] << "\t" << r_location[1] << "\t" << is_boundary_node << std::endl;
    //    }
    //
    //    // Now write an additional node at each VertexElement's centroid
    //    unsigned num_tetrahedral_elements = 0;
    //    for (unsigned vertex_elem_index=0; vertex_elem_index<num_vertex_elements; vertex_elem_index++)
    //    {
    //        unsigned index = num_vertex_nodes + vertex_elem_index;
    //
    //        c_vector<double, DIM> location = mpMutableVertexMesh->GetCentroidOfElement(vertex_elem_index);
    //
    //        // Any node located at a VertexElement's centroid will not be a boundary node
    //        unsigned is_boundary_node = 0;
    //        (*p_node_file) << index << "\t" << location[0] << "\t" << location[1] << "\t" << is_boundary_node << std::endl;
    //
    //        // Also keep track of how many tetrahedral elements there will be
    //        num_tetrahedral_elements += mpMutableVertexMesh->GetElement(vertex_elem_index)->GetNumNodes();
    //    }
    //    p_node_file->close();
    //
    //    // Write element file
    //    out_stream p_elem_file = output_file_handler.OpenOutputFile(mesh_file_name+".ele");
    //    (*p_elem_file) << std::scientific;
    //    (*p_elem_file) << num_tetrahedral_elements << "\t3\t0" << std::endl;
    //
    //    std::set<std::pair<unsigned, unsigned> > tetrahedral_edges;
    //
    //    unsigned tetrahedral_elem_index = 0;
    //    for (unsigned vertex_elem_index=0; vertex_elem_index<num_vertex_elements; vertex_elem_index++)
    //    {
    //        VertexElement<DIM, DIM>* p_vertex_element = mpMutableVertexMesh->GetElement(vertex_elem_index);
    //
    //        // Iterate over nodes owned by this VertexElement
    //        unsigned num_nodes_in_vertex_element = p_vertex_element->GetNumNodes();
    //        for (unsigned local_index=0; local_index<num_nodes_in_vertex_element; local_index++)
    //        {
    //            unsigned node_0_index = p_vertex_element->GetNodeGlobalIndex(local_index);
    //            unsigned node_1_index = p_vertex_element->GetNodeGlobalIndex((local_index+1)%num_nodes_in_vertex_element);
    //            unsigned node_2_index = num_vertex_nodes + vertex_elem_index;
    //
    //            (*p_elem_file) << tetrahedral_elem_index++ << "\t" << node_0_index << "\t" << node_1_index << "\t" << node_2_index << std::endl;
    //
    //            // Add edges to the set if they are not already present
    //            std::pair<unsigned, unsigned> edge_0 = this->CreateOrderedPair(node_0_index, node_1_index);
    //            std::pair<unsigned, unsigned> edge_1 = this->CreateOrderedPair(node_1_index, node_2_index);
    //            std::pair<unsigned, unsigned> edge_2 = this->CreateOrderedPair(node_2_index, node_0_index);
    //
    //            tetrahedral_edges.insert(edge_0);
    //            tetrahedral_edges.insert(edge_1);
    //            tetrahedral_edges.insert(edge_2);
    //        }
    //    }
    //    p_elem_file->close();
    //
    //    // Write edge file
    //    out_stream p_edge_file = output_file_handler.OpenOutputFile(mesh_file_name+".edge");
    //    (*p_edge_file) << std::scientific;
    //    (*p_edge_file) << tetrahedral_edges.size() << "\t1" << std::endl;
    //
    //    unsigned edge_index = 0;
    //    for (std::set<std::pair<unsigned, unsigned> >::iterator edge_iter = tetrahedral_edges.begin();
    //         edge_iter != tetrahedral_edges.end();
    //         ++edge_iter)
    //    {
    //        std::pair<unsigned, unsigned> this_edge = *edge_iter;
    //
    //        // To be a boundary edge both nodes need to be boundary nodes.
    //        bool is_boundary_edge = false;
    //        if (this_edge.first  < mpMutableVertexMesh->GetNumNodes() &&
    //            this_edge.second  < mpMutableVertexMesh->GetNumNodes())
    //        {
    //            is_boundary_edge = (mpMutableVertexMesh->GetNode(this_edge.first)->IsBoundaryNode() &&
    //                                mpMutableVertexMesh->GetNode(this_edge.second)->IsBoundaryNode() );
    //        }
    //        unsigned is_boundary_edge_unsigned = is_boundary_edge ? 1 : 0;
    //
    //        (*p_edge_file) << edge_index++ << "\t" << this_edge.first << "\t" << this_edge.second << "\t" << is_boundary_edge_unsigned << std::endl;
    //    }
    //    p_edge_file->close();
    //
    //    // Having written the mesh to file, now construct it using TrianglesMeshReader
    //    TetrahedralMesh<DIM, DIM>* p_mesh = new TetrahedralMesh<DIM, DIM>;
    //
    //    // Nested scope so reader is destroyed before we remove the temporary files
    //    {
    //        TrianglesMeshReader<DIM, DIM> mesh_reader(output_dir + mesh_file_name);
    //        p_mesh->ConstructFromMeshReader(mesh_reader);
    //    }
    //
    //    // Delete the temporary files
    //    output_file_handler.FindFile("").Remove();
    //
    //    // The original files have been deleted, it is better if the mesh object forgets about them
    //    p_mesh->SetMeshHasChangedSinceLoading();
    //
    //    return p_mesh;

    return NULL;
}

template <unsigned DIM>
bool ImmersedBoundaryCellPopulation<DIM>::IsPdeNodeAssociatedWithNonApoptoticCell(unsigned pdeNodeIndex)
{
    //    bool non_apoptotic_cell_present = true;
    //
    //    if (pdeNodeIndex < this->GetNumNodes())
    //    {
    //        std::set<unsigned> containing_element_indices = this->GetNode(pdeNodeIndex)->rGetContainingElementIndices();
    //
    //        for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
    //             iter != containing_element_indices.end();
    //             iter++)
    //        {
    //            if (this->GetCellUsingLocationIndex(*iter)->template HasCellProperty<ApoptoticCellProperty>() )
    //            {
    //                non_apoptotic_cell_present = false;
    //                break;
    //            }
    //        }
    //    }
    //    else
    //    {
    //        /*
    //         * This node of the tetrahedral finite element mesh is in the centre of the element of the
    //         * vertex-based cell population, so we can use an offset to compute which cell to interrogate.
    //         */
    //        non_apoptotic_cell_present = !(this->GetCellUsingLocationIndex(pdeNodeIndex - this->GetNumNodes())->template HasCellProperty<ApoptoticCellProperty>());
    //    }
    //
    //    return non_apoptotic_cell_present;

    return true;
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetCellDataItemAtPdeNode(unsigned pdeNodeIndex,
                                                                     std::string& rVariableName,
                                                                     bool dirichletBoundaryConditionApplies,
                                                                     double dirichletBoundaryValue)
{
    //    unsigned num_nodes = this->GetNumNodes();
    //    double value = 0.0;
    //
    //    // Cells correspond to nodes in the centre of the vertex element; nodes on vertices have averaged values from containing cells
    //
    //    if (pdeNodeIndex >= num_nodes)
    //    {
    //        // Offset to relate elements in vertex mesh to nodes in tetrahedral mesh
    //        assert(pdeNodeIndex-num_nodes < num_nodes);
    //
    //        CellPtr p_cell = this->GetCellUsingLocationIndex(pdeNodeIndex - num_nodes);
    //        value = p_cell->GetCellData()->GetItem(rVariableName);
    //    }
    //    else
    //    {
    //        ///\todo Work out a better way to do the nodes not associated with cells
    //        if (dirichletBoundaryConditionApplies)
    //        {
    //            // We need to impose the Dirichlet boundaries again here as not represented in cell data
    //            value = dirichletBoundaryValue;
    //        }
    //        else
    //        {
    //            assert(pdeNodeIndex < num_nodes);
    //            Node<DIM>* p_node = this->GetNode(pdeNodeIndex);
    //
    //            // Average over data from containing elements (cells)
    //            std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
    //            for (std::set<unsigned>::iterator index_iter = containing_elements.begin();
    //                 index_iter != containing_elements.end();
    //                 ++index_iter)
    //            {
    //                assert(*index_iter < num_nodes);
    //                CellPtr p_cell = this->GetCellUsingLocationIndex(*index_iter);
    //                value += p_cell->GetCellData()->GetItem(rVariableName);
    //            }
    //            value /= containing_elements.size();
    //        }
    //    }
    //
    //    return value;

    return 0.0;
}

template <unsigned DIM>
boost::shared_ptr<AbstractVertexBasedDivisionRule<DIM> > ImmersedBoundaryCellPopulation<DIM>::GetVertexBasedDivisionRule()
{
    return mpVertexBasedDivisionRule;
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetVertexBasedDivisionRule(boost::shared_ptr<AbstractVertexBasedDivisionRule<DIM> > pVertexBasedDivisionRule)
{
    mpVertexBasedDivisionRule = pVertexBasedDivisionRule;
}

template <unsigned DIM>
bool ImmersedBoundaryCellPopulation<DIM>::DoesPopulationHaveActiveSources()
{
    return mPopulationHasActiveSources;
}

template <unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetIfPopulationHasActiveSources(bool hasActiveSources)
{
    mPopulationHasActiveSources = hasActiveSources;
}

template <unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetDefaultTimeStep()
{
    return 0.002;
}

// Explicit instantiation
template class ImmersedBoundaryCellPopulation<1>;
template class ImmersedBoundaryCellPopulation<2>;
template class ImmersedBoundaryCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryCellPopulation)
