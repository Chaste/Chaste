/*

Copyright (c) 2005-2015, University of Oxford.
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
#include <boost/foreach.hpp>
#include "ImmersedBoundaryMeshWriter.hpp"
#include "Warnings.hpp"
#include "ChasteSyscalls.hpp"
#include "IsNan.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include <boost/multi_array.hpp>
#include "Debug.hpp"


// Cell population writers
#include "CellPopulationElementWriter.hpp"

#include "RandomNumberGenerator.hpp"


template<unsigned DIM>
ImmersedBoundaryCellPopulation<DIM>::ImmersedBoundaryCellPopulation(ImmersedBoundaryMesh<DIM, DIM>& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDeleteMesh(deleteMesh)
{
    mpImmersedBoundaryMesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>* >(&(this->mrMesh));
    mpVertexBasedDivisionRule.reset(new ShortAxisVertexBasedDivisionRule<DIM>());

    mInteractionDistance = 0.05 * sqrt(mpImmersedBoundaryMesh->GetVolumeOfElement(mpImmersedBoundaryMesh->GetNumElements()-1));

    // If no location indices are specified, associate with elements from the mesh (assumed to be sequentially ordered).
    std::list<CellPtr>::iterator it = this->mCells.begin();
    for (unsigned i=0; it != this->mCells.end(); ++it, ++i)
    {
        unsigned index = locationIndices.empty() ? i : locationIndices[i]; // assume that the ordering matches
        AbstractCellPopulation<DIM, DIM>::AddCellUsingLocationIndex(index,*it);
    }

    // Check each element has only one cell attached
    if (validate)
    {
        Validate();
    }

    // Default active sources to true (safer - simulation will always work with sources set to zero, just a bit slower)
    mPopulationHasActiveSources = true;

    // Set the intrinsic spacing to a default 0.01
    //\todo should this be static?
    mIntrinsicSpacing = 0.01;
}

template<unsigned DIM>
ImmersedBoundaryCellPopulation<DIM>::ImmersedBoundaryCellPopulation(ImmersedBoundaryMesh<DIM, DIM>& rMesh)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh),
      mDeleteMesh(true)
{
    mpImmersedBoundaryMesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>* >(&(this->mrMesh));
}

template<unsigned DIM>
ImmersedBoundaryCellPopulation<DIM>::~ImmersedBoundaryCellPopulation()
{
    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    return 0.0;
}

template<unsigned DIM>
ImmersedBoundaryMesh<DIM, DIM>& ImmersedBoundaryCellPopulation<DIM>::rGetMesh()
{
    return *mpImmersedBoundaryMesh;
}

template<unsigned DIM>
const ImmersedBoundaryMesh<DIM, DIM>& ImmersedBoundaryCellPopulation<DIM>::rGetMesh() const
{
    return *mpImmersedBoundaryMesh;
}

template<unsigned DIM>
ImmersedBoundaryElement<DIM, DIM>* ImmersedBoundaryCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mpImmersedBoundaryMesh->GetElement(elementIndex);
}

template<unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
c_vector<double, DIM> ImmersedBoundaryCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mpImmersedBoundaryMesh->GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

template<unsigned DIM>
Node<DIM>* ImmersedBoundaryCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetInteractionDistance(double new_distance)
{
    assert(new_distance >= 0.0);
    mInteractionDistance = new_distance;
}

template<unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetInteractionDistance()
{
    return mInteractionDistance;
}

template<unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetIntrinsicSpacing()
{
    return mIntrinsicSpacing;
}


//\todo: implement this method. Decide what "neighbouring" should be for IB cells
template<unsigned DIM>
std::set<unsigned> ImmersedBoundaryCellPopulation<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    // The set to return
    std::set<unsigned> neighbouring_indices;
    return neighbouring_indices;
}

template<unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    std::cout << "Trying to add new node from within cell population\n";
    return 0;//mpImmersedBoundaryMesh->AddNode(pNewNode);
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mpImmersedBoundaryMesh->SetNode(nodeIndex, rNewLocation);
}

template<unsigned DIM>
ImmersedBoundaryElement<DIM, DIM>* ImmersedBoundaryCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mpImmersedBoundaryMesh->GetElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
unsigned ImmersedBoundaryCellPopulation<DIM>::GetNumElements()
{
    return mpImmersedBoundaryMesh->GetNumElements();
}

template<unsigned DIM>
CellPtr ImmersedBoundaryCellPopulation<DIM>::AddCell(CellPtr pNewCell,
                                                const c_vector<double,DIM>& rCellDivisionVector,
                                                CellPtr pParentCell)
{
    /*
    // Get the element associated with this cell
    ImmersedBoundaryElement<DIM, DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element
    unsigned new_element_index = mpImmersedBoundaryMesh->DivideElementAlongGivenAxis(p_element,
                                                                                  rCellDivisionVector,
                                                                                  true);
    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(new_element_index,p_created_cell);
    this->mCellLocationMap[p_created_cell.get()] = new_element_index;
    return p_created_cell;*/

    //\todo: may need to implement this when cells start dividing
    CellPtr p_created_cell = this->mCells.back();
    return p_created_cell;
}

template<unsigned DIM>
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

template<unsigned DIM>
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
        for (unsigned i = 0 ; i < 4 ; i ++)
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
        ChastePoint <DIM> new_point(node_location);

        // Move the node
        this->SetNode(node_iter->GetIndex(), new_point);
    }
    // If active sources, we need to update those location as well
    if (this->DoesPopulationHaveActiveSources())
    {
        std::vector<FluidSource<DIM> *> &r_element_sources = this->rGetMesh().rGetElementFluidSources();
        std::vector<FluidSource<DIM> *> &r_balance_sources = this->rGetMesh().rGetBalancingFluidSources();

        // Construct a vector of all sources combined
        std::vector<FluidSource<DIM> *> combined_sources;
        combined_sources.insert(combined_sources.end(), r_element_sources.begin(), r_element_sources.end());
        combined_sources.insert(combined_sources.end(), r_balance_sources.begin(), r_balance_sources.end());

        c_vector<double, DIM> source_location;

        // Iterate over all sources and update their locations
        for (unsigned source_idx = 0 ; source_idx < combined_sources.size() ; source_idx++)
        {
            // Get location of current node
            source_location = combined_sources[source_idx]->rGetLocation();

            // Get first grid index in each dimension, taking account of possible wrap-around
            first_idx_x = unsigned(floor(source_location[0] / grid_spacing_x)) + num_grid_pts_x - 1;
            first_idx_y = unsigned(floor(source_location[1] / grid_spacing_y)) + num_grid_pts_y - 1;

            // Calculate all four indices and deltas in each dimension
            for (unsigned i = 0 ; i < 4 ; i ++)
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
}

template<unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::Delta1D(double dist, double spacing)
{
    return (0.25 * (1.0 + cos(M_PI * dist / (2 * spacing))));
}

template<unsigned DIM>
bool ImmersedBoundaryCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();;
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    // I don't think this is needed for IB
}

template<unsigned DIM>
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

    for (unsigned i=0; i<validated_element.size(); i++)
    {
        if (validated_element[i] == 0)
        {
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() <<", Element " << i << " does not appear to have a cell associated with it");
        }

        if (validated_element[i] > 1)
        {
            // This should never be reached as you can only set one cell per element index
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() <<", Element " << i << " appears to have " << validated_element[i] << " cells associated with it");
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    //pPopulationWriter->Visit(this);
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter)
{
    //pPopulationCountWriter->Visit(this);
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    //pCellWriter->VisitCell(pCell, this);
}

template<unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Get the vertex element index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get the cell's volume from the vertex mesh
    double cell_volume = mpImmersedBoundaryMesh->GetVolumeOfElement(elem_index);

    return cell_volume;
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
//#ifdef CHASTE_VTK
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

        // Iterate over vertex elements ///\todo #2512 - replace with loop over cells
        for (typename ImmersedBoundaryMesh<DIM,DIM>::ImmersedBoundaryElementIterator elem_iter = mpImmersedBoundaryMesh->GetElementIteratorBegin();
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
            for (unsigned cell_part = 0 ; cell_part < num_cell_parts[elem_index] ; cell_part++)
            {
                vtk_cell_data.push_back((*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this));
            }
        }

        mesh_writer.AddCellData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    }

    // When outputting any CellData, we assume that the first cell is representative of all cells
    unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

    std::vector<std::vector<double> > cell_data;
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cell_data_var;
        cell_data.push_back(cell_data_var);
    }

    // Loop over vertex elements ///\todo #2512 - replace with loop over cells
    for (typename ImmersedBoundaryMesh<DIM,DIM>::ImmersedBoundaryElementIterator elem_iter = mpImmersedBoundaryMesh->GetElementIteratorBegin();
         elem_iter != mpImmersedBoundaryMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in the vertex mesh
        unsigned elem_index = elem_iter->GetIndex();

        // Get the cell corresponding to this element
        CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
        assert(p_cell);

        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            // Populate the vector of VTK cell data.  We loop over the number of output cells as this takes into
            // account that some elements will be broken into pieces for visualisation
            for (unsigned cell_part = 0 ; cell_part < num_cell_parts[elem_index] ; cell_part++)
            {
                cell_data[var].push_back(p_cell->GetCellData()->GetItem(cell_data_names[var]));
            }
        }
    }
    for (unsigned var=0; var<num_cell_data_items; var++)
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
//#endif //CHASTE_VTK
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellPopulationElementWriter>())
        {
            this-> template AddPopulationWriter<CellPopulationElementWriter>();
        }
    }

    AbstractCellPopulation<DIM>::OpenWritersFiles(rOutputFileHandler);
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Add the division rule parameters
    *rParamsFile << "\t\t<VertexBasedDivisionRule>\n";
    mpVertexBasedDivisionRule->OutputCellVertexBasedDivisionRuleInfo(rParamsFile);
    *rParamsFile << "\t\t</VertexBasedDivisionRule>\n";

    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
double ImmersedBoundaryCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = this->mrMesh.GetWidth(rDimension);

    return width;
}


//\todo: implement this.  May need to put method back in to mesh class
template<unsigned DIM>
std::set<unsigned> ImmersedBoundaryCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    std::set<unsigned> neighbouring_indices;
    return neighbouring_indices;
}

template<unsigned DIM>
boost::shared_ptr<AbstractVertexBasedDivisionRule<DIM> > ImmersedBoundaryCellPopulation<DIM>::GetVertexBasedDivisionRule()
{
    return mpVertexBasedDivisionRule;
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetVertexBasedDivisionRule(boost::shared_ptr<AbstractVertexBasedDivisionRule<DIM> > pVertexBasedDivisionRule)
{
    mpVertexBasedDivisionRule = pVertexBasedDivisionRule;
}

template<unsigned DIM>
bool ImmersedBoundaryCellPopulation<DIM>::DoesPopulationHaveActiveSources()
{
    return mPopulationHasActiveSources;
}

template<unsigned DIM>
void ImmersedBoundaryCellPopulation<DIM>::SetIfPopulationHasActiveSources(bool hasActiveSources)
{
    mPopulationHasActiveSources = hasActiveSources;
}

// Explicit instantiation
template class ImmersedBoundaryCellPopulation<1>;
template class ImmersedBoundaryCellPopulation<2>;
template class ImmersedBoundaryCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryCellPopulation)
