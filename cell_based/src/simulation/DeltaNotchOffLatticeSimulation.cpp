/*

Copyright (c) 2005-2012, University of Oxford.
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
void DeltaNotchOffLatticeSimulation<DIM>::UpdateAtEndOfTimeStep()
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
