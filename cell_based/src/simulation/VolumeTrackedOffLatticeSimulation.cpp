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

#include "VolumeTrackedOffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CellwiseData.hpp"

template<unsigned DIM>
VolumeTrackedOffLatticeSimulation<DIM>::VolumeTrackedOffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                      bool deleteCellPopulationInDestructor,
                                                      bool initialiseCells)
    : OffLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{
}

template<unsigned DIM>
VolumeTrackedOffLatticeSimulation<DIM>::~VolumeTrackedOffLatticeSimulation()
{
}

template<unsigned DIM>
void VolumeTrackedOffLatticeSimulation<DIM>::SetupSolve()
{
    /*
     * We must update CellwiseData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellwiseData();
}

template<unsigned DIM>
void VolumeTrackedOffLatticeSimulation<DIM>::PostSolve()
{
    UpdateCellwiseData();
}

template<unsigned DIM>
void VolumeTrackedOffLatticeSimulation<DIM>::UpdateCellwiseData()
{
    // Make sure the cell population is updated
    this->mrCellPopulation.Update();

    // Prepare CellwiseData by reallocating memory according to the number of cells
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    /**
     * This hack is needed because in the case of a MeshBasedCellPopulation in which
     * multiple cell divisions have occurred over one time step, the Voronoi tessellation
     * (while existing) is out-of-date. Thus, if we did not regenerate the Voronoi
     * tessellation here, an assertion may trip as we try to access a Voronoi element
     * whose index exceeds the number of elements in the out-of-date tessellation.
     * 
     * \todo work out how to properly fix this (#1986)
     */
    if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->CreateVoronoiTessellation();
    }

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        // Get the volume of this cell
        double cell_volume = this->mrCellPopulation.GetVolumeOfCell(*cell_iter);

        // Get the location index corresponding to this cell
        unsigned location_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Store the cell's volume in CellwiseData
        CellwiseData<DIM>::Instance()->SetValue(cell_volume, location_index, 0);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VolumeTrackedOffLatticeSimulation<1>;
template class VolumeTrackedOffLatticeSimulation<2>;
template class VolumeTrackedOffLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeTrackedOffLatticeSimulation)
