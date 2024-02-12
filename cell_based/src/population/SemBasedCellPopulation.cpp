/*

Copyright (c) 2005-2023, University of Oxford.
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
#include "VtkMeshWriter.hpp"

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
    
    unsigned index = 0;
    for (auto cell_ptr : this->mCells) {
        this->AddCellUsingLocationIndex(index, cell_ptr);
        index++;
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
    //unsigned elem_index = this->GetLocationIndexUsingCell(pCell);
    return {}; //this->rGetMesh().GetNeighbouringElementIndices(elem_index);
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
void SemBasedCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
    SemMeshWriter<DIM> mesh_writer(rDirectory, "results", false);
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;
    mesh_writer.WriteVtkUsingMesh(*mpSemMesh, time.str());
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
    return false;
}
template<unsigned DIM>
CellPtr SemBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{
    auto node_index = this->mCells.size();

    // Associate the new cell with the neighbouring node
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->AddCellUsingLocationIndex(node_index,p_created_cell);

    return p_created_cell;
}
template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetDefaultTimeStep()
{
    return 0.0;
}
template<unsigned DIM>
unsigned SemBasedCellPopulation<DIM>::RemoveDeadCells()
{
    return 0;
}
template<unsigned DIM>
void SemBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
}
template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    return 0.0;
}
template<unsigned DIM>
std::set<unsigned> SemBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    return {};
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
    
}

template<unsigned DIM>
double SemBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    return 1.0;
}

// Explicit instantiation
template class SemBasedCellPopulation<1>;
template class SemBasedCellPopulation<2>;
template class SemBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemBasedCellPopulation)