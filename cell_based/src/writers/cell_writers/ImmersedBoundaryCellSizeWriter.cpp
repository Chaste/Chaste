//
// Created by bartmanski on 02/05/17.
//

#include "ImmersedBoundaryCellSizeWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryCellSizeWriter<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryCellSizeWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("ib_cell_size.dat")
{
    this->mVtkCellDataName = "Neighbour number";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryCellSizeWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return static_cast<double>(pCellPopulation->GetVolumeOfCell(pCell));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryCellSizeWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);

    double cell_size = this->GetCellDataForVtkOutput(pCell, pCellPopulation);

    *this->mpOutStream << location_index << " " << cell_id << " ";
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << centre_location[i] << " ";
    }

    *this->mpOutStream << cell_size << " ";
}

// Explicit instantiation
template class ImmersedBoundaryCellSizeWriter<1,1>;
template class ImmersedBoundaryCellSizeWriter<1,2>;
template class ImmersedBoundaryCellSizeWriter<2,2>;
template class ImmersedBoundaryCellSizeWriter<1,3>;
template class ImmersedBoundaryCellSizeWriter<2,3>;
template class ImmersedBoundaryCellSizeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ImmersedBoundaryCellSizeWriter)
