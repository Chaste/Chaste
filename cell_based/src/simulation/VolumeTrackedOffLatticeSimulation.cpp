/*

Copyright (C) University of Oxford, 2005-2011

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
#include "VertexBasedCellPopulation.hpp"
#include "CellwiseData.hpp"

template<unsigned DIM>
VolumeTrackedOffLatticeSimulation<DIM>::VolumeTrackedOffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                      bool deleteCellPopulationInDestructor,
                                                      bool initialiseCells)
    : OffLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{
    if (!dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation) && !dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) )
    {
        EXCEPTION("VolumeTrackedOffLatticeSimulation require a subclass of MeshBasedCellPopulation or VertexBasedSimulation.");
    }
}

template<unsigned DIM>
VolumeTrackedOffLatticeSimulation<DIM>::~VolumeTrackedOffLatticeSimulation()
{
}

template<unsigned DIM>
void VolumeTrackedOffLatticeSimulation<DIM>::PostSolve()
{
    // Make sure the cell population is updated
    this->mrCellPopulation.Update();
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        // Static cast on the cell population
        MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation));

        // Create Voronoi tessellation for volumes
        p_static_cast_cell_population->CreateVoronoiTessellation();

        // Loop over cells and set volume value in CellWiseData
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
             cell_iter != this->mrCellPopulation.End();
             ++cell_iter)
        {
            // Get the index of the node corresponding to this cell
            unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // Store in CellwiseData
            CellwiseData<DIM>::Instance()->SetValue(p_static_cast_cell_population->GetVolumeOfVoronoiElement(node_index), node_index, 0);
        }
    }
    else if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        // Static cast on the cell population
        VertexBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&(this->mrCellPopulation));

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
             cell_iter != this->mrCellPopulation.End();
             ++cell_iter)
        {

            // Get the index of vertex element corresponding to this cell
            unsigned element_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // Store in CellwiseData
            CellwiseData<DIM>::Instance()->SetValue(p_static_cast_cell_population->rGetMesh().GetVolumeOfElement(element_index), element_index, 0);
        }
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
