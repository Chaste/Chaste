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

    mInteractionDistance = 2 * mpImmersedBoundaryMesh->GetCharacteristicNodeSpacing();

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
    // Get references to the fluid velocity grid
    const std::vector<std::vector<double> >& velocity_grid_x = this->rGetMesh().rGetFluidVelocityGridX();
    const std::vector<std::vector<double> >& velocity_grid_y = this->rGetMesh().rGetFluidVelocityGridY();

    /**
    * Set up all necessary variables before loop for efficiency
    */

    // Vectors for displacement and location of node
    c_vector<double, DIM> displacement;
    c_vector<double, DIM> node_location;

    // Get the size of the mesh and set up force grids
    unsigned num_grid_pts_x = this->rGetMesh().GetNumGridPtsX();
    unsigned num_grid_pts_y = this->rGetMesh().GetNumGridPtsY();

    // Helper variables
    double characteristic_spacing = mpImmersedBoundaryMesh->GetCharacteristicNodeSpacing();
    double step_size_x = 1.0 / (double)num_grid_pts_x;
    double step_size_y = 1.0 / (double)num_grid_pts_y;
    double dist_x;
    double dist_y;
    int first_idx_x;
    int first_idx_y;

    // Iterate over all nodes
    for (typename ImmersedBoundaryMesh<DIM, DIM>::NodeIterator node_iter = this->rGetMesh().GetNodeIteratorBegin(false);
            node_iter != this->rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
    {
        // Get location of current node
        node_location = node_iter->rGetLocation();

        // Get first index in fluid velocity grid grid (x and y)
        first_idx_x = (int) floor(node_location[0] / step_size_x) - 1;
        first_idx_y = (int) floor(node_location[1] / step_size_y) - 1;

        // Loop over the 4x4 grid which will influence the displacement of the current node
        for (unsigned x_idx = 0; x_idx < 4; x_idx++)
        {
            // Calculate distance between current x index and node, then account for possible wrap-around
            dist_x = fabs((double) (first_idx_x + x_idx) * step_size_x - node_location[0]);
            if (first_idx_x == -1)
            {
                first_idx_x += num_grid_pts_x;
            }

            for (unsigned y_idx = 0; y_idx < 4; y_idx++)
            {
                // Calculate distance between current y index and node, then account for possible wrap-around
                dist_y = fabs((double) (first_idx_y + y_idx) * step_size_y - node_location[1]);
                if (first_idx_y == -1)
                {
                    first_idx_y += num_grid_pts_y;
                }

                // The applied velocity is weighted by the delta function
                displacement[0] += velocity_grid_x[(first_idx_y + y_idx) % num_grid_pts_y][(first_idx_x + x_idx) % num_grid_pts_x] * Delta1D(dist_x, step_size_x) * Delta1D(dist_y, step_size_y);
                displacement[1] += velocity_grid_y[(first_idx_y + y_idx) % num_grid_pts_y][(first_idx_x + x_idx) % num_grid_pts_x] * Delta1D(dist_x, step_size_x) * Delta1D(dist_y, step_size_y);
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

// Explicit instantiation
template class ImmersedBoundaryCellPopulation<1>;
template class ImmersedBoundaryCellPopulation<2>;
template class ImmersedBoundaryCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryCellPopulation)
