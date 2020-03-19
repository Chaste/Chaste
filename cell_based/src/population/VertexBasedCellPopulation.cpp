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

#include "VertexBasedCellPopulation.hpp"
#include "Warnings.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"
#include "StepSizeException.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "SmartPointers.hpp"
#include "T2SwapCellKiller.hpp"
#include "ApoptoticCellProperty.hpp"
#include "SrnCellModel.hpp"
#include "CellPopulationElementWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT2SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"
#include "AbstractCellBasedSimulation.hpp"

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDeleteMesh(deleteMesh),
      mOutputCellRearrangementLocations(true),
      mRestrictVertexMovement(true)
{
    mpMutableVertexMesh = static_cast<MutableVertexMesh<DIM, DIM>* >(&(this->mrMesh));
    mpVertexBasedDivisionRule.reset(new ShortAxisVertexBasedDivisionRule<DIM>());

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

    // If cells contain an SRN model, then we need to track mesh operations
    // and update SRNs accordingly. Here we assume the first cell is a representative of other cells
    if ((*this->mCells.begin())->HasSrnModel())
    {
        mPopulationSrn.SetVertexCellPopulation(this);
        mpMutableVertexMesh->SetMeshOperationTracking(true);
    }
}

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh),
      mDeleteMesh(true),
      mOutputCellRearrangementLocations(true),
      mRestrictVertexMovement(true)
{
    mpMutableVertexMesh = static_cast<MutableVertexMesh<DIM, DIM>* >(&(this->mrMesh));
}

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::~VertexBasedCellPopulation()
{
    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    // Take the average of the cells containing this vertex
    double average_damping_constant = 0.0;

    std::set<unsigned> containing_elements = GetNode(nodeIndex)->rGetContainingElementIndices();

    unsigned num_containing_elements = containing_elements.size();
    if (num_containing_elements == 0)
    {
        EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Node " << nodeIndex << " is not contained in any elements, so GetDampingConstant() returns zero");
    }

    double temp = 1.0/((double) num_containing_elements);
    for (std::set<unsigned>::iterator iter = containing_elements.begin();
         iter != containing_elements.end();
         ++iter)
    {
        CellPtr p_cell = this->GetCellUsingLocationIndex(*iter);
        bool cell_is_wild_type = p_cell->GetMutationState()->IsType<WildTypeCellMutationState>();

        if (cell_is_wild_type)
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
    return *mpMutableVertexMesh;
}

template<unsigned DIM>
const MutableVertexMesh<DIM, DIM>& VertexBasedCellPopulation<DIM>::rGetMesh() const
{
    return *mpMutableVertexMesh;
}

template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mpMutableVertexMesh->GetElement(elementIndex);
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
c_vector<double, DIM> VertexBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mpMutableVertexMesh->GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

template<unsigned DIM>
Node<DIM>* VertexBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
std::set<unsigned> VertexBasedCellPopulation<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);
    return this->rGetMesh().GetNeighbouringElementIndices(elem_index);
}

template<unsigned int DIM>
std::set<std::pair<unsigned int, unsigned int>>
VertexBasedCellPopulation<DIM>::GetNeighbouringEdgeIndices(CellPtr pCell, unsigned EdgeLocalIndex)
{
    std::set<std::pair<unsigned, unsigned>> neighbours;
    auto cellLocationIndex = this->GetLocationIndexUsingCell(pCell);
    auto element = this->GetElement(cellLocationIndex);
    auto globalEdgeIndex = element->GetEdgeGlobalIndex(EdgeLocalIndex);
    auto neighbourElementIndices = element->GetNeighbouringElementAtEdgeIndex(EdgeLocalIndex);

    // Normally there is only one neighbouring element
    for (auto neighbourElementIndex : neighbourElementIndices)
    {
        auto neighbourElement = this->GetElement(neighbourElementIndex);
        // Iterate over neighbouring element indices
        for (unsigned eIndex = 0; eIndex < neighbourElement->GetNumEdges(); eIndex++)
        {
            // If the neighbours edge matches EdgeLocalIndex
            if (neighbourElement->GetEdge(eIndex)->GetIndex() == globalEdgeIndex)
            {
                neighbours.insert(std::pair<unsigned, unsigned>(neighbourElementIndex, eIndex));
            }
        }
    }

    return neighbours;
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mpMutableVertexMesh->AddNode(pNewNode);
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mpMutableVertexMesh->SetNode(nodeIndex, rNewLocation);
}

template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mpMutableVertexMesh->GetElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::GetNumElements()
{
    return mpMutableVertexMesh->GetNumElements();
}

template<unsigned DIM>
CellPtr VertexBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    // Get the element associated with this cell
    VertexElement<DIM, DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Get the orientation of division
    c_vector<double, DIM> division_vector = mpVertexBasedDivisionRule->CalculateCellDivisionVector(pParentCell, *this);

    // Divide the element
    unsigned new_element_index = mpMutableVertexMesh->DivideElementAlongGivenAxis(p_element, division_vector, true);
    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(new_element_index,p_created_cell);
    this->mCellLocationMap[p_created_cell.get()] = new_element_index;

    return p_created_cell;
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator it = this->mCells.begin();
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
                mpMutableVertexMesh->DeleteElementPriorToReMesh(this->GetLocationIndexUsingCell((*it)));
            }

            // Delete the cell
            it = this->mCells.erase(it);
        }
        else
        {
            ++it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::CheckForStepSizeException(unsigned nodeIndex, c_vector<double,DIM>& rDisplacement, double dt)
{
    double length = norm_2(rDisplacement);

    if (mRestrictVertexMovement)
    {
        if (length > 0.5*mpMutableVertexMesh->GetCellRearrangementThreshold())
        {
            rDisplacement *= 0.5*mpMutableVertexMesh->GetCellRearrangementThreshold()/length;

            std::ostringstream message;
            message << "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted ";
            message << "so the motion has been restricted. Use a smaller timestep to avoid these warnings.";

            double suggested_step = 0.95*dt*((0.5*mpMutableVertexMesh->GetCellRearrangementThreshold())/length);

            // The first time we see this behaviour, throw a StepSizeException, but not more than once
            if (mThrowStepSizeException)
            {
                mThrowStepSizeException = false;
                throw StepSizeException(suggested_step, message.str(), false);
            }
        }
    }
}

template<unsigned DIM>
bool VertexBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    VertexElementMap element_map(mpMutableVertexMesh->GetNumAllElements());
    mpMutableVertexMesh->ReMesh(element_map);

    if (!element_map.IsIdentityMap())
    {
        // Fix up the mappings between CellPtrs and VertexElements
        ///\todo We want to make these maps private, so we need a better way of doing the code below.
        std::map<Cell*, unsigned> old_map = this->mCellLocationMap;

        this->mCellLocationMap.clear();
        this->mLocationCellMap.clear();

        for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            // The cell vector should only ever contain living cells
            unsigned old_elem_index = old_map[(*cell_iter).get()];
            assert(!element_map.IsDeleted(old_elem_index));

            unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);
            this->SetCellUsingLocationIndex(new_elem_index, *cell_iter);
        }

        // Check that each VertexElement has only one CellPtr associated with it in the updated cell population
        Validate();
    }
    //First cell is representative of other cells
    bool EdgeModelOrNot = (*this->mCells.begin())->GetSrnModel()->HasEdgeModel();
    if (EdgeModelOrNot)
    {
        // Note that SRN update after divisions is handled through Cell::Divide() method
        mPopulationSrn.UpdateSrnAfterBirthOrDeath(element_map);
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
void VertexBasedCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter)
{
    pPopulationCountWriter->Visit(this);
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::GetRosetteRankOfCell(CellPtr pCell)
{
    // Get the vertex element index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get the element rosette rank from the vertex mesh
    unsigned rosette_rank = mpMutableVertexMesh->GetRosetteRankOfElement(elem_index);

    return rosette_rank;
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Get the vertex element index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get the cell's volume from the vertex mesh
    double cell_volume = mpMutableVertexMesh->GetVolumeOfElement(elem_index);

    return cell_volume;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{

    auto cells = this->rGetCells();
    boost::shared_ptr<CellEdgeData> p_cell_edge_data = (*cells.begin())->GetCellEdgeData();
    //If edge SRNs are specified, then write vtk results into a mesh where quantities
    //associated with each edge are taken into account. We assume that the first cell is
    //representative of all cells.
    //Edge VTKs are also written if cells contain CellEdgeData
    if (cells.size() > 0)
    {
        //If cells contain edge data
        if (p_cell_edge_data != nullptr)
        {
            this->WriteCellEdgeVtkResultsToFile(rDirectory);
            return;
        }
    }

    this->WriteCellVtkResultsToFile(rDirectory);

}

template<unsigned int DIM>
void VertexBasedCellPopulation<DIM>::WriteCellVtkResultsToFile(const std::string &rDirectory)
{
#ifdef CHASTE_VTK

    // Create mesh writer for VTK output
    VertexMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results", false);

    // Iterate over any cell writers that are present
    unsigned num_cells = this->GetNumAllCells();
    for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
            cell_writer_iter != this->mCellWriters.end();
            ++cell_writer_iter)
    {
        // Create vector to store VTK cell data
        std::vector<double> vtk_cell_data(num_cells);

        // Iterate over vertex elements ///\todo #2512 - replace with loop over cells
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
                elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
                ++elem_iter)
        {
            // Get index of this element in the vertex mesh
            unsigned elem_index = elem_iter->GetIndex();

            // Get the cell corresponding to this element
            CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
            assert(p_cell);

            // Populate the vector of VTK cell data
            vtk_cell_data[elem_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
        }

        mesh_writer.AddCellData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
    }

    // When outputting any CellData, we assume that the first cell is representative of all cells
    unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

    std::vector<std::vector<double> > cell_data;
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cell_data_var(num_cells);
        cell_data.push_back(cell_data_var);
    }

    // Loop over vertex elements ///\todo #2512 - replace with loop over cells
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
            elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
            ++elem_iter)
    {
        // Get index of this element in the vertex mesh
        unsigned elem_index = elem_iter->GetIndex();

        // Get the cell corresponding to this element
        CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
        assert(p_cell);

        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cell_data[var][elem_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
        }
    }
    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        mesh_writer.AddCellData(cell_data_names[var], cell_data[var]);
    }

    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    mesh_writer.WriteVtkUsingMesh(*mpMutableVertexMesh, time.str());

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";

#endif //CHASTE_VTK
}

template<unsigned int DIM>
void VertexBasedCellPopulation<DIM>::WriteCellEdgeVtkResultsToFile(const std::string &rDirectory)
{
#ifdef CHASTE_VTK
    //Writes cell only data
    {
        // Create mesh writer for VTK output
        VertexMeshWriter<DIM, DIM> mesh_writer(rDirectory, "cell_results", false);

        // Iterate over any cell writers that are present
        unsigned num_cells = this->GetNumAllCells();
        for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = this->mCellWriters.begin();
                cell_writer_iter != this->mCellWriters.end();
                ++cell_writer_iter)
        {
            // Create vector to store VTK cell data
            std::vector<double> vtk_cell_data(num_cells);

            // Iterate over vertex elements ///\todo #2512 - replace with loop over cells
            for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
                    elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
                    ++elem_iter)
            {
                // Get index of this element in the vertex mesh
                unsigned elem_index = elem_iter->GetIndex();

                // Get the cell corresponding to this element
                CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
                assert(p_cell);

                // Populate the vector of VTK cell data
                vtk_cell_data[elem_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this);
            }

            mesh_writer.AddCellData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
        }

        // When outputting any CellData, we assume that the first cell is representative of all cells
        unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
        std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

        std::vector<std::vector<double> > cell_data;
        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            std::vector<double> cell_data_var(num_cells);
            cell_data.push_back(cell_data_var);
        }

        // Loop over vertex elements ///\todo #2512 - replace with loop over cells
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
                elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
                ++elem_iter)
        {
            // Get index of this element in the vertex mesh
            unsigned elem_index = elem_iter->GetIndex();

            // Get the cell corresponding to this element
            CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
            assert(p_cell);

            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cell_data[var][elem_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
            }
        }
        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            mesh_writer.AddCellData(cell_data_names[var], cell_data[var]);
        }

        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
        std::stringstream time;
        time << num_timesteps;

        mesh_writer.WriteVtkUsingMesh(*mpMutableVertexMesh, time.str());

        *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"cell_results_";
        *(this->mpVtkMetaFile) << num_timesteps;
        *(this->mpVtkMetaFile) << ".vtu\"/>\n";
    }
    // Create mesh writer for VTK output
    TrapezoidEdgeVertexMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results", false);
    unsigned num_edges = 0;
    //here elements are synonymous with cells
    const unsigned n_cells = this->GetNumElements();
    //Similarly as in TrapEdgeVerteMeshWriter, but instead of nodes
    //we fill edge arrays
    //The first value MUST be zero
    std::vector<unsigned> cell_offset_dist(n_cells);
    //The order of stored data is illustrated below:
    // [_____|_][____|_]
    //  ^^^^  ^
    // edge   cell interior
    for (unsigned i=1; i<n_cells; ++i)
    {
        cell_offset_dist[i] = cell_offset_dist[i-1]+this->GetElement(i-1)->GetNumEdges()+1;
        num_edges += this->GetElement(i)->GetNumEdges();
    }
    //Total number of edges
    num_edges+=this->GetElement(0)->GetNumEdges();

    // Iterate over any cell writers that are present
    // The data that is written below is associated with the entire cell
    // E.g. the cell age. Thus, edges also share the same characteristic
    for (auto cell_writer: this->mCellWriters)
    {
        // Create vector to store VTK cell data
        std::vector<double> vtk_cell_data(num_edges+n_cells);
        // Iterate over vertex elements ///\todo #2512 - replace with loop over cells
        for (auto elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
             elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
             ++elem_iter)
        {
            // Get index of this element in the vertex mesh
            unsigned elem_index = elem_iter->GetIndex();
            // Get the cell corresponding to this element
            CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
            assert(p_cell);
            // Write the same property in all the edges
            auto element = this->GetElement(elem_index);
            unsigned elem_num_edges = element->GetNumEdges();
            //Edge data
            for (unsigned e = 0; e < elem_num_edges; ++e)
            {
                // Populate the vector of VTK cell data
                vtk_cell_data[cell_offset_dist[elem_index]+e]
                              = cell_writer->GetCellDataForVtkOutput(p_cell, this);
            }
            //Internal cell data
            vtk_cell_data[cell_offset_dist[elem_index]+elem_num_edges]
                          =cell_writer->GetCellDataForVtkOutput(p_cell, this);
        }
        mesh_writer.AddCellData(cell_writer->GetVtkCellDataName(), vtk_cell_data);
    }
    //When outputting CellData and CellEdgeData, we assume that the first cell
    //and its edges are representative of the population.
    //Get cell/edge data names
    const unsigned int n_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

    const unsigned int n_edge_data_items = this->Begin()->GetCellEdgeData()->GetNumItems();
    std::vector<std::string> edge_data_names = this->Begin()->GetCellEdgeData()->GetKeys();

    //Total number of data items. Each data item (edge+interior) has its own value
    const unsigned int n_data_items = n_edge_data_items + n_cell_data_items;
    std::vector<std::string> data_names(n_data_items);
    for (unsigned int var=0; var<n_edge_data_items; ++var)
    {
        data_names[var] = edge_data_names[var];
    }
    for (unsigned int var=n_edge_data_items; var<n_data_items; ++var)
    {
        data_names[var] = cell_data_names[var-n_edge_data_items];
    }
    std::vector<std::vector<double> > data_values(n_data_items,
                                                  std::vector<double>(num_edges + n_cells));

    // Writing CellEdgeData values to the edges of the cells
    // Loop over vertex elements ///\todo #2512 - replace with loop over cells
    for (auto elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
            elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
            ++elem_iter)
    {
        const unsigned int elem_index = elem_iter ->GetIndex();
        CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
        assert(p_cell);
        //Edge data for interior data is set to zero.
        auto element = this->GetElement(elem_index);
        unsigned elem_num_edges = element->GetNumEdges();
        for (unsigned int var = 0; var < n_edge_data_items; ++var)
        {
            // Write the same property in all the edges
            for (unsigned int e = 0; e < elem_num_edges; ++e)
                data_values[var][cell_offset_dist[elem_index]+e]
                                 = p_cell->GetCellEdgeData()->GetItem(data_names[var])[e];
            //Cell interior is set to zero
            data_values[var][cell_offset_dist[elem_index]+elem_num_edges] = 0.0;
        }
    }

    // Writing CellData values to the interior of the cells
    // Loop over vertex elements ///\todo #2512 - replace with loop over cells
    for (auto elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
         elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in the vertex mesh
        unsigned elem_index = elem_iter->GetIndex();
        unsigned elem_num_edges = elem_iter->GetNumEdges();

        // Get the cell corresponding to this element
        CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
        assert(p_cell);

        for (unsigned int var = n_edge_data_items; var<n_data_items; ++var)
        {
            //Data in the edges are set to 0
            for (unsigned int e = 0; e < elem_num_edges; ++e)
            {
                data_values[var][cell_offset_dist[elem_index]+e] = 0.0;
            }
            //Filling cell interior data
            data_values[var][cell_offset_dist[elem_index]+elem_num_edges]
                             = p_cell->GetCellData()->GetItem(data_names[var]);

        }
    }

    for (unsigned var=0; var<n_data_items; var++)
    {
        mesh_writer.AddCellData(data_names[var], data_values[var]);
    }

    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    mesh_writer.WriteVtkUsingMesh(*mpMutableVertexMesh, time.str());

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::OpenWritersFiles(OutputFileHandler& rOutputFileHandler)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellPopulationElementWriter>())
        {
            this-> template AddPopulationWriter<CellPopulationElementWriter>();
        }
    }

    if (mOutputCellRearrangementLocations)
    {
        if (!this-> template HasWriter<VertexT1SwapLocationsWriter>())
        {
            this-> template AddPopulationWriter<VertexT1SwapLocationsWriter>();
        }
        if (!this-> template HasWriter<VertexT2SwapLocationsWriter>())
        {
            this-> template AddPopulationWriter<VertexT2SwapLocationsWriter>();
        }
        if (!this-> template HasWriter<VertexT3SwapLocationsWriter>())
        {
            this-> template AddPopulationWriter<VertexT3SwapLocationsWriter>();
        }
    }


    AbstractCellPopulation<DIM>::OpenWritersFiles(rOutputFileHandler);
}


template<unsigned DIM>
bool VertexBasedCellPopulation<DIM>::GetOutputCellRearrangementLocations()
{
    return mOutputCellRearrangementLocations;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetOutputCellRearrangementLocations(bool outputCellRearrangementLocations)
{
    mOutputCellRearrangementLocations = outputCellRearrangementLocations;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<CellRearrangementThreshold>" << mpMutableVertexMesh->GetCellRearrangementThreshold() << "</CellRearrangementThreshold>\n";
    *rParamsFile << "\t\t<T2Threshold>" <<  mpMutableVertexMesh->GetT2Threshold() << "</T2Threshold>\n";
    *rParamsFile << "\t\t<CellRearrangementRatio>" << mpMutableVertexMesh->GetCellRearrangementRatio() << "</CellRearrangementRatio>\n";
    *rParamsFile << "\t\t<OutputCellRearrangementLocations>" << mOutputCellRearrangementLocations << "</OutputCellRearrangementLocations>\n";

    // Add the division rule parameters
    *rParamsFile << "\t\t<VertexBasedDivisionRule>\n";
    mpVertexBasedDivisionRule->OutputCellVertexBasedDivisionRuleInfo(rParamsFile);
    *rParamsFile << "\t\t</VertexBasedDivisionRule>\n";

    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = this->mrMesh.GetWidth(rDimension);

    return width;
}

template<unsigned DIM>
std::set<unsigned> VertexBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    return mpMutableVertexMesh->GetNeighbouringNodeIndices(index);
}

template<unsigned DIM>
boost::shared_ptr<AbstractVertexBasedDivisionRule<DIM> > VertexBasedCellPopulation<DIM>::GetVertexBasedDivisionRule()
{
    return mpVertexBasedDivisionRule;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetVertexBasedDivisionRule(boost::shared_ptr<AbstractVertexBasedDivisionRule<DIM> > pVertexBasedDivisionRule)
{
    mpVertexBasedDivisionRule = pVertexBasedDivisionRule;
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetTetrahedralMeshForPdeModifier()
{
    // This method only works in 2D sequential
    assert(DIM == 2);                        // LCOV_EXCL_LINE - disappears at compile time.
    assert(PetscTools::IsSequential());

    unsigned num_vertex_nodes = mpMutableVertexMesh->GetNumNodes();
    unsigned num_vertex_elements = mpMutableVertexMesh->GetNumElements();

    std::string mesh_file_name = "mesh";

    // Get a unique temporary foldername
    std::stringstream pid;
    pid << getpid();
    OutputFileHandler output_file_handler("2D_temporary_tetrahedral_mesh_" + pid.str());
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

    // Compute the number of nodes in the TetrahedralMesh
    unsigned num_tetrahedral_nodes = num_vertex_nodes + num_vertex_elements;

    // Write node file
    out_stream p_node_file = output_file_handler.OpenOutputFile(mesh_file_name+".node");
    (*p_node_file) << std::scientific;
    (*p_node_file) << std::setprecision(20);
    (*p_node_file) << num_tetrahedral_nodes << "\t2\t0\t1" << std::endl;

    // Begin by writing each node in the VertexMesh
    for (unsigned node_index=0; node_index<num_vertex_nodes; node_index++)
    {
        Node<DIM>* p_node = mpMutableVertexMesh->GetNode(node_index);

        ///\todo will the nodes in mpMutableVertexMesh always have indices 0,1,2,...? (#2221)
        unsigned index = p_node->GetIndex();
        const c_vector<double, DIM>& r_location = p_node->rGetLocation();
        unsigned is_boundary_node = p_node->IsBoundaryNode() ? 1 : 0;

        (*p_node_file) << index << "\t" << r_location[0] << "\t" << r_location[1] << "\t" << is_boundary_node << std::endl;
    }

    // Now write an additional node at each VertexElement's centroid
    unsigned num_tetrahedral_elements = 0;
    for (unsigned vertex_elem_index=0; vertex_elem_index<num_vertex_elements; vertex_elem_index++)
    {
        unsigned index = num_vertex_nodes + vertex_elem_index;

        c_vector<double, DIM> location = mpMutableVertexMesh->GetCentroidOfElement(vertex_elem_index);

        // Any node located at a VertexElement's centroid will not be a boundary node
        unsigned is_boundary_node = 0;
        (*p_node_file) << index << "\t" << location[0] << "\t" << location[1] << "\t" << is_boundary_node << std::endl;

        // Also keep track of how many tetrahedral elements there will be
        num_tetrahedral_elements += mpMutableVertexMesh->GetElement(vertex_elem_index)->GetNumNodes();
    }
    p_node_file->close();

    // Write element file
    out_stream p_elem_file = output_file_handler.OpenOutputFile(mesh_file_name+".ele");
    (*p_elem_file) << std::scientific;
    (*p_elem_file) << num_tetrahedral_elements << "\t3\t0" << std::endl;

    std::set<std::pair<unsigned, unsigned> > tetrahedral_edges;

    unsigned tetrahedral_elem_index = 0;
    for (unsigned vertex_elem_index=0; vertex_elem_index<num_vertex_elements; vertex_elem_index++)
    {
        VertexElement<DIM, DIM>* p_vertex_element = mpMutableVertexMesh->GetElement(vertex_elem_index);

        // Iterate over nodes owned by this VertexElement
        unsigned num_nodes_in_vertex_element = p_vertex_element->GetNumNodes();
        for (unsigned local_index=0; local_index<num_nodes_in_vertex_element; local_index++)
        {
            unsigned node_0_index = p_vertex_element->GetNodeGlobalIndex(local_index);
            unsigned node_1_index = p_vertex_element->GetNodeGlobalIndex((local_index+1)%num_nodes_in_vertex_element);
            unsigned node_2_index = num_vertex_nodes + vertex_elem_index;

            (*p_elem_file) << tetrahedral_elem_index++ << "\t" << node_0_index << "\t" << node_1_index << "\t" << node_2_index << std::endl;

            // Add edges to the set if they are not already present
            std::pair<unsigned, unsigned> edge_0 = this->CreateOrderedPair(node_0_index, node_1_index);
            std::pair<unsigned, unsigned> edge_1 = this->CreateOrderedPair(node_1_index, node_2_index);
            std::pair<unsigned, unsigned> edge_2 = this->CreateOrderedPair(node_2_index, node_0_index);

            tetrahedral_edges.insert(edge_0);
            tetrahedral_edges.insert(edge_1);
            tetrahedral_edges.insert(edge_2);
        }
    }
    p_elem_file->close();

    // Write edge file
    out_stream p_edge_file = output_file_handler.OpenOutputFile(mesh_file_name+".edge");
    (*p_edge_file) << std::scientific;
    (*p_edge_file) << tetrahedral_edges.size() << "\t1" << std::endl;

    unsigned edge_index = 0;
    for (std::set<std::pair<unsigned, unsigned> >::iterator edge_iter = tetrahedral_edges.begin();
         edge_iter != tetrahedral_edges.end();
         ++edge_iter)
    {
        std::pair<unsigned, unsigned> this_edge = *edge_iter;

        // To be a boundary edge both nodes need to be boundary nodes.
        bool is_boundary_edge = false;
        if (this_edge.first  < mpMutableVertexMesh->GetNumNodes() &&
            this_edge.second  < mpMutableVertexMesh->GetNumNodes())
        {
            is_boundary_edge = (mpMutableVertexMesh->GetNode(this_edge.first)->IsBoundaryNode() &&
                                mpMutableVertexMesh->GetNode(this_edge.second)->IsBoundaryNode() );
        }
        unsigned is_boundary_edge_unsigned = is_boundary_edge ? 1 : 0;

        (*p_edge_file) << edge_index++ << "\t" << this_edge.first << "\t" << this_edge.second << "\t" << is_boundary_edge_unsigned << std::endl;
    }
    p_edge_file->close();

    // Having written the mesh to file, now construct it using TrianglesMeshReader
    TetrahedralMesh<DIM, DIM>* p_mesh = new TetrahedralMesh<DIM, DIM>;

    // Nested scope so reader is destroyed before we remove the temporary files
    {
        TrianglesMeshReader<DIM, DIM> mesh_reader(output_dir + mesh_file_name);
        p_mesh->ConstructFromMeshReader(mesh_reader);
    }

    // Delete the temporary files
    output_file_handler.FindFile("").Remove();

    // The original files have been deleted, it is better if the mesh object forgets about them
    p_mesh->SetMeshHasChangedSinceLoading();

    return p_mesh;
}

template<unsigned DIM>
bool VertexBasedCellPopulation<DIM>::IsPdeNodeAssociatedWithNonApoptoticCell(unsigned pdeNodeIndex)
{
    bool non_apoptotic_cell_present = true;

    if (pdeNodeIndex < this->GetNumNodes())
    {
        std::set<unsigned> containing_element_indices = this->GetNode(pdeNodeIndex)->rGetContainingElementIndices();

        for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
             iter != containing_element_indices.end();
             iter++)
        {
            if (this->GetCellUsingLocationIndex(*iter)->template HasCellProperty<ApoptoticCellProperty>() )
            {
                non_apoptotic_cell_present = false;
                break;
            }
        }
    }
    else
    {
        /*
         * This node of the tetrahedral finite element mesh is in the centre of the element of the
         * vertex-based cell population, so we can use an offset to compute which cell to interrogate.
         */
        non_apoptotic_cell_present = !(this->GetCellUsingLocationIndex(pdeNodeIndex - this->GetNumNodes())->template HasCellProperty<ApoptoticCellProperty>());
    }

    return non_apoptotic_cell_present;
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetCellDataItemAtPdeNode(
        unsigned pdeNodeIndex,
        std::string& rVariableName,
        bool dirichletBoundaryConditionApplies,
        double dirichletBoundaryValue)
{
    unsigned num_nodes = this->GetNumNodes();
    double value = 0.0;

    // Cells correspond to nodes in the centre of the vertex element; nodes on vertices have averaged values from containing cells

    if (pdeNodeIndex >= num_nodes)
    {
        // Offset to relate elements in vertex mesh to nodes in tetrahedral mesh
        assert(pdeNodeIndex-num_nodes < num_nodes);

        CellPtr p_cell = this->GetCellUsingLocationIndex(pdeNodeIndex - num_nodes);
        value = p_cell->GetCellData()->GetItem(rVariableName);
    }
    else
    {
        ///\todo Work out a better way to do the nodes not associated with cells
        if (dirichletBoundaryConditionApplies)
        {
            // We need to impose the Dirichlet boundaries again here as not represented in cell data
            value = dirichletBoundaryValue;
        }
        else
        {
            assert(pdeNodeIndex < num_nodes);
            Node<DIM>* p_node = this->GetNode(pdeNodeIndex);

            // Average over data from containing elements (cells)
            std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            for (std::set<unsigned>::iterator index_iter = containing_elements.begin();
                 index_iter != containing_elements.end();
                 ++index_iter)
            {
                assert(*index_iter < num_nodes);
                CellPtr p_cell = this->GetCellUsingLocationIndex(*index_iter);
                value += p_cell->GetCellData()->GetItem(rVariableName);
            }
            value /= containing_elements.size();
        }
    }

    return value;
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetDefaultTimeStep()
{
    return 0.002;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
    if (bool(dynamic_cast<Cylindrical2dVertexMesh*>(&(this->mrMesh))))
    {
        *pVizSetupFile << "MeshWidth\t" << this->GetWidth(0) << "\n";
    }
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SimulationSetupHook(AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    MAKE_PTR_ARGS(T2SwapCellKiller<DIM>, p_t2_swap_cell_killer, (this));
    pSimulation->AddCellKiller(p_t2_swap_cell_killer);
}


template<unsigned DIM>
bool VertexBasedCellPopulation<DIM>::GetRestrictVertexMovementBoolean()
{
    return mRestrictVertexMovement;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetRestrictVertexMovementBoolean(bool restrictMovement)
{
    mRestrictVertexMovement = restrictMovement;
}

template<unsigned DIM>
VertexBasedPopulationSrn<DIM>& VertexBasedCellPopulation<DIM>::rGetVertexBasedPopulationSrn()
{
    return mPopulationSrn;
}

// Explicit instantiation
template class VertexBasedCellPopulation<1>;
template class VertexBasedCellPopulation<2>;
template class VertexBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedCellPopulation)
