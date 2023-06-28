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
    return this->rGetMesh().GetNeighbouringElementIndices(elem_index);
}

template<unsigned DIM>
unsigned SemBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mpSemMesh->AddNode(pNewNode);
}

template<unsigned DIM>
void SemBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mpSemMesh->SetNode(nodeIndex, rNewLocation);
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

// Explicit instantiation
template class SemBasedCellPopulation<1>;
template class SemBasedCellPopulation<2>;
template class SemBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemBasedCellPopulation)