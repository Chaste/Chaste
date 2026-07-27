/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "SemBasedCellPopulation.hpp"
#include "StepSizeException.hpp"
#include "VtkMeshWriter.hpp"
#include "WildTypeCellMutationState.hpp"

#include <list>
#include <set>
#include <sstream>

template<unsigned DIM>
SemBasedCellPopulation<DIM>::SemBasedCellPopulation(SemMesh<DIM>& rMesh,
                           std::vector<CellPtr>& rCells,
                           bool deleteMesh,
                           bool validate,
                           const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDeleteMesh(deleteMesh)
{
    mpSemMesh = static_cast<SemMesh<DIM>*>(&(this->mrMesh));

    if (this->mCells.size() != mpSemMesh->GetNumElements())
    {
        EXCEPTION("There must be precisely one CellPtr for each SemElement");
    }

    std::set<unsigned> used_location_indices;
    std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
    for (unsigned cell_index=0; cell_iter != this->mCells.end(); ++cell_iter, ++cell_index)
    {
        unsigned location_index = locationIndices.empty() ? cell_index : locationIndices[cell_index];
        if (location_index >= mpSemMesh->GetNumElements())
        {
            EXCEPTION("A supplied location index does not correspond to a SemElement");
        }
        if (!used_location_indices.insert(location_index).second)
        {
            EXCEPTION("A SemElement location index is assigned to more than one cell");
        }
        this->AddCellUsingLocationIndex(location_index, *cell_iter);
    }

    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
SemBasedCellPopulation<DIM>::SemBasedCellPopulation(SemMesh<DIM>& rMesh)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh),
      mDeleteMesh(true)
{
    mpSemMesh = static_cast<SemMesh<DIM>*>(&(this->mrMesh));
}

template<unsigned DIM>
SemBasedCellPopulation<DIM>::~SemBasedCellPopulation()
{
    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
SemMesh<DIM>& SemBasedCellPopulation<DIM>::rGetMesh()
{
    return *mpSemMesh;
}

template<unsigned DIM>
const SemMesh<DIM>& SemBasedCellPopulation<DIM>::rGetMesh() const
{
    return *mpSemMesh;
}

template<unsigned DIM>
SemElement<DIM>* SemBasedCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mpSemMesh->GetElement(elementIndex);
}

template<unsigned DIM>
unsigned SemBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
unsigned SemBasedCellPopulation<DIM>::GetNumElements()
{
    return mpSemMesh->GetNumElements();
}


template<unsigned DIM>
c_vector<double, DIM> SemBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mpSemMesh->GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

template<unsigned DIM>
Node<DIM>* SemBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
std::set<unsigned> SemBasedCellPopulation<DIM>::GetNeighbouringLocationIndices(CellPtr pCell)
{
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);
    std::set<unsigned> neighbouring_element_indices;

    for (const auto& r_node_pair : this->mNodePairs)
    {
        const std::set<unsigned>& r_first_element_indices = r_node_pair.first->rGetContainingElementIndices();
        const std::set<unsigned>& r_second_element_indices = r_node_pair.second->rGetContainingElementIndices();

        if (r_first_element_indices.find(elem_index) != r_first_element_indices.end())
        {
            for (unsigned candidate_index : r_second_element_indices)
            {
                if (candidate_index != elem_index && !this->GetElement(candidate_index)->IsDeleted())
                {
                    neighbouring_element_indices.insert(candidate_index);
                }
            }
        }

        if (r_second_element_indices.find(elem_index) != r_second_element_indices.end())
        {
            for (unsigned candidate_index : r_first_element_indices)
            {
                if (candidate_index != elem_index && !this->GetElement(candidate_index)->IsDeleted())
                {
                    neighbouring_element_indices.insert(candidate_index);
                }
            }
        }
    }

    return neighbouring_element_indices;
}

template<unsigned DIM>
unsigned SemBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mpSemMesh->AddNode(pNewNode);
}

template<unsigned DIM>
void SemBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    this->GetNode(nodeIndex)->SetPoint(rNewLocation);
}

template<unsigned DIM>
SemElement<DIM>* SemBasedCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mpSemMesh->GetElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Get the SemElement index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get the cell's volume from the SemMesh
    double cell_volume = mpSemMesh->GetVolumeOfElement(elem_index);

    return cell_volume;
}

template<unsigned DIM>
void SemBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void SemBasedCellPopulation<DIM>::Validate()
{
    // Check each SemElement has only one cell attached
    std::vector<unsigned> validated_element = std::vector<unsigned>(this->GetNumElements(), 0);
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);
        if (elem_index >= this->GetNumElements())
        {
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() <<", Cell is associated with element index " << elem_index << ", which is outside the SemMesh");
        }
        if (this->GetElement(elem_index)->IsDeleted())
        {
            EXCEPTION("At time " << SimulationTime::Instance()->GetTime() <<", Cell is associated with deleted element " << elem_index);
        }
        validated_element[elem_index]++;
    }

    for (unsigned i=0; i<validated_element.size(); i++)
    {
        if (this->GetElement(i)->IsDeleted())
        {
            continue;
        }

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
void SemBasedCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    SemMeshWriter<DIM> mesh_writer(rDirectory, "results", false);

    for (const auto& p_writer : this->mNodePointDataWriters)
    {
        mesh_writer.AddPointData(p_writer->rGetFieldName(), p_writer->GetPointData(this));
    }

    // Write per-cell quantities as per-element VTK cell data. One cell maps to one SemElement, so
    // each array has one value per element, indexed by element (location) index; SemMeshWriter
    // expands these to the element's point-cloud and surface VTK cells.
    unsigned num_elements = mpSemMesh->GetNumElements();
    if (num_elements > 0)
    {
        // Any registered cell writers (e.g. CellVolumesWriter, which uses GetVolumeOfCell())
        for (const auto& p_cell_writer : this->mCellWriters)
        {
            std::vector<double> vtk_element_data(num_elements);
            for (typename SemMesh<DIM>::SemElementIterator elem_iter = mpSemMesh->GetElementIteratorBegin();
                 elem_iter != mpSemMesh->GetElementIteratorEnd();
                 ++elem_iter)
            {
                unsigned elem_index = elem_iter->GetIndex();
                CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
                vtk_element_data[elem_index] = p_cell_writer->GetCellDataForVtkOutput(p_cell, this);
            }
            mesh_writer.AddElementData(p_cell_writer->GetVtkCellDataName(), vtk_element_data);
        }

        // Any CellData items (assuming the first cell is representative of all cells)
        unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
        std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

        std::vector<std::vector<double> > cell_data(num_cell_data_items, std::vector<double>(num_elements));
        for (typename SemMesh<DIM>::SemElementIterator elem_iter = mpSemMesh->GetElementIteratorBegin();
             elem_iter != mpSemMesh->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
            for (unsigned var = 0; var < num_cell_data_items; ++var)
            {
                cell_data[var][elem_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
            }
        }
        for (unsigned var = 0; var < num_cell_data_items; ++var)
        {
            mesh_writer.AddElementData(cell_data_names[var], cell_data[var]);
        }
    }

    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;
    mesh_writer.WriteVtkUsingMesh(*mpSemMesh, time.str());

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#else
    (void)rDirectory;
#endif //CHASTE_VTK
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* SemBasedCellPopulation<DIM>::GetTetrahedralMeshForPdeModifier()
{
    return nullptr;
}

template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetCellDataItemAtPdeNode(unsigned pdeNodeIndex,std::string& item, bool, double)
{
    return 0.0;
}

template<unsigned DIM>
bool SemBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    unsigned location_index = this->GetLocationIndexUsingCell(pCell);
    return this->GetElement(location_index)->IsDeleted();
}
template<unsigned DIM>
CellPtr SemBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    (void)pNewCell;
    (void)pParentCell;
    EXCEPTION("SemBasedCellPopulation does not support AddCell() because SEM element division is not implemented");
    return CellPtr();
}
template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetDefaultTimeStep()
{
    return 0.002;
}
template<unsigned DIM>
unsigned SemBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         )
    {
        if ((*cell_iter)->IsDead())
        {
            unsigned location_index = this->GetLocationIndexUsingCell(*cell_iter);
            this->GetElement(location_index)->MarkAsDeleted();
            this->RemoveCellUsingLocationIndex(location_index, *cell_iter);
            cell_iter = this->mCells.erase(cell_iter);
            num_removed++;
        }
        else
        {
            ++cell_iter;
        }
    }

    return num_removed;
}

template<unsigned DIM>
void SemBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    this->mNodePairs.clear();
    mpSemMesh->UpdateBoxCollection();
    mpSemMesh->CalculateNodePairs(this->mNodePairs);
}

template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    return this->mrMesh.GetWidth(rDimension);
}
template<unsigned DIM>
std::set<unsigned> SemBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    std::set<unsigned> neighbouring_node_indices;
    for (const auto& r_node_pair : this->mNodePairs)
    {
        if (r_node_pair.first->GetIndex() == index)
        {
            neighbouring_node_indices.insert(r_node_pair.second->GetIndex());
        }
        if (r_node_pair.second->GetIndex() == index)
        {
            neighbouring_node_indices.insert(r_node_pair.first->GetIndex());
        }
    }
    return neighbouring_node_indices;
}
template<unsigned DIM>
void SemBasedCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
}
template<unsigned DIM>
void SemBasedCellPopulation<DIM>::AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter)
{
    pPopulationCountWriter->Visit(this);
}
template<unsigned DIM>
void SemBasedCellPopulation<DIM>::AcceptPopulationEventWriter(boost::shared_ptr<AbstractCellPopulationEventWriter<DIM, DIM> > pPopulationEventWriter)
{
    pPopulationEventWriter->Visit(this);
}
template<unsigned DIM>
void SemBasedCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}
template<unsigned DIM>
void SemBasedCellPopulation<DIM>::CheckForStepSizeException(unsigned nodeIndex, c_vector<double,DIM>& rDisplacement, double dt)
{
    // A node moving further than the box-collection interaction distance in a single step
    // invalidates the neighbour search (potential interactions would be missed), so treat
    // this as a step-size error rather than allowing a silent node explosion.
    double length = norm_2(rDisplacement);
    double max_interaction_distance = mpSemMesh->GetMaximumInteractionDistance();

    if (length > max_interaction_distance)
    {
        std::ostringstream message;
        message << "Node " << nodeIndex << " is moving by " << length
                << ", which is more than the SEM maximum interaction distance (" << max_interaction_distance
                << "): use a smaller timestep or diffusion constant to avoid this exception.";

        // Suggest a net time step that will give a movement smaller than the interaction distance
        double new_step = 0.95*dt*(max_interaction_distance/length);

        throw StepSizeException(new_step, message.str(), true); // terminate
    }
}

template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    std::set<unsigned> containing_elements = GetNode(nodeIndex)->rGetContainingElementIndices();
    double damping_constant = 0.0;
    unsigned num_live_containing_elements = 0u;

    for (unsigned elem_index : containing_elements)
    {
        if (elem_index < this->GetNumElements() && !this->GetElement(elem_index)->IsDeleted())
        {
            CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
            bool cell_is_wild_type = p_cell->GetMutationState()->IsType<WildTypeCellMutationState>();
            damping_constant += cell_is_wild_type ? this->GetDampingConstantNormal() : this->GetDampingConstantMutant();
            num_live_containing_elements++;
        }
    }

    if (num_live_containing_elements == 0u)
    {
        EXCEPTION("At time " << SimulationTime::Instance()->GetTime() << ", Node " << nodeIndex << " is not contained in any live SEM elements, so GetDampingConstant() returns zero");
    }

    return damping_constant/static_cast<double>(num_live_containing_elements);
}

// Explicit instantiation
template class SemBasedCellPopulation<1>;
template class SemBasedCellPopulation<2>;
template class SemBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemBasedCellPopulation)
