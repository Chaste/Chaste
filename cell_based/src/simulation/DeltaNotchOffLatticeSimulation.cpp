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

#include "DeltaNotchOffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "DeltaNotchCellCycleModel.hpp"

template<unsigned DIM>
DeltaNotchOffLatticeSimulation<DIM>::DeltaNotchOffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                                  bool deleteCellPopulationInDestructor,
                                                                  bool initialiseCells)
    : OffLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{
}

template<unsigned DIM>
DeltaNotchOffLatticeSimulation<DIM>::~DeltaNotchOffLatticeSimulation()
{
}

template<unsigned DIM>
void DeltaNotchOffLatticeSimulation<DIM>::SetupSolve()
{
    UpdateCellwiseData();
}

template<unsigned DIM>
void DeltaNotchOffLatticeSimulation<DIM>::PostSolve()
{
    UpdateCellwiseData();
}

template<unsigned DIM>
void DeltaNotchOffLatticeSimulation<DIM>::UpdateCellwiseData()
{
    // Make sure the cell population is updated
    this->mrCellPopulation.Update();

    // Prepare CellwiseData by reallocating memory according to the number of cells
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // First store each cell's Notch and Delta concentrations in CellwiseData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        unsigned index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        DeltaNotchCellCycleModel* p_model = static_cast<DeltaNotchCellCycleModel*>(cell_iter->GetCellCycleModel());
        double this_delta = p_model->GetDelta();
        double this_notch = p_model->GetNotch();

        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        CellwiseData<DIM>::Instance()->SetValue(this_notch, index, 0);
        CellwiseData<DIM>::Instance()->SetValue(this_delta, index, 1);
    }

    // Next iterate over the population to compute and store each cell's neighbouring Delta concentration in CellwiseData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        // Get the location index corresponding to this cell
        unsigned index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Get the set of neighbouring location indices            
        std::set<unsigned> neighbour_indices;
        if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
        {
            neighbour_indices = this->mrCellPopulation.GetNeighbouringNodeIndices(index);
        }
        else
        {
            neighbour_indices = static_cast<VertexBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetMesh().GetNeighbouringElementIndices(index);
        }

        // Compute this cell's average neighbouring Delta concentration and store in CellwiseData
        if (!neighbour_indices.empty())
        {
            double mean_delta = 0.0;
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {
                CellPtr p_cell = this->mrCellPopulation.GetCellUsingLocationIndex(*iter);
                double this_delta = CellwiseData<DIM>::Instance()->GetValue(p_cell, 1);
                mean_delta += this_delta/neighbour_indices.size();
            }
            CellwiseData<DIM>::Instance()->SetValue(mean_delta, index, 2);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DeltaNotchOffLatticeSimulation<1>;
template class DeltaNotchOffLatticeSimulation<2>;
template class DeltaNotchOffLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchOffLatticeSimulation)
