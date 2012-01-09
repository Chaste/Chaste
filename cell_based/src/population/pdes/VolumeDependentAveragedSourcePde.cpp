/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "VolumeDependentAveragedSourcePde.hpp"
#include "ApoptoticCellProperty.hpp"

template<unsigned DIM>
VolumeDependentAveragedSourcePde<DIM>::VolumeDependentAveragedSourcePde(AbstractCellPopulation<DIM>& rCellPopulation, double coefficient)
    : AveragedSourcePde<DIM>(rCellPopulation, coefficient)
{
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)));
    mpStaticCastCellPopulation = static_cast<NodeBasedCellPopulation<DIM>*>(&(this->mrCellPopulation));
}

template<unsigned DIM>
void VolumeDependentAveragedSourcePde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh,  std::map< CellPtr, unsigned >* pCellPdeElementMap) // must be called before solve
{
    // Allocate memory
    this->mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<this->mCellDensityOnCoarseElements.size(); elem_index++)
    {
        this->mCellDensityOnCoarseElements[elem_index] = 0.0;
    }

    // Loop over cells, find which coarse element it is in, and add volume to the mSourceTermOnCoarseElements[elem_index];
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        unsigned elem_index = 0;
        const ChastePoint<DIM>& r_position_of_cell = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

        if (pCellPdeElementMap != NULL)
        {
            elem_index = (*pCellPdeElementMap)[*cell_iter];
        }
        else
        {
            elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
        }

        unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Uptake normalised to 1 for unit cell
        double radius = mpStaticCastCellPopulation->rGetMesh().GetCellRadius(node_index);
        double cell_weight = radius*radius;

        bool cell_is_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();
        if (!cell_is_apoptotic)
        {
            this->mCellDensityOnCoarseElements[elem_index] += cell_weight;
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<this->mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        this->mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VolumeDependentAveragedSourcePde<1>;
template class VolumeDependentAveragedSourcePde<2>;
template class VolumeDependentAveragedSourcePde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeDependentAveragedSourcePde)
