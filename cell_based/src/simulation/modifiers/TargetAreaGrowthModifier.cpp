/*

Copyright (c) 2005-2014, University of Oxford.
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

#include "TargetAreaGrowthModifier.hpp"

template<unsigned DIM>
TargetAreaGrowthModifier<DIM>::TargetAreaGrowthModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mMatureCellTargetArea(1.0)
{
}

template<unsigned DIM>
TargetAreaGrowthModifier<DIM>::~TargetAreaGrowthModifier()
{
}

template<unsigned DIM>
void TargetAreaGrowthModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateTargetAreas(rCellPopulation);
}

template<unsigned DIM>
void TargetAreaGrowthModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateTargetAreas(rCellPopulation);
}

template<unsigned DIM>
void TargetAreaGrowthModifier<DIM>::UpdateTargetAreas(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    ///\todo #2488: double check that this update call doesn't break anything (i.e. counting of swaps etc.)
    rCellPopulation.Update();

    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("TargetAreaGrowthModifier is to be used with a VertexBasedCellPopulation only");
    }

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        UpdateTargetAreaOfCell( p_cell_population->GetCellUsingLocationIndex(elem_index) );
    }
}

template<unsigned DIM>
void TargetAreaGrowthModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = mMatureCellTargetArea;

    double cell_age = pCell->GetAge();
    double g1_duration = pCell->GetCellCycleModel()->GetG1Duration();

    // If the cell is differentiated then its G1 duration is infinite
    if (g1_duration == DBL_MAX) // don't use magic number, compare to DBL_MAX
    {
        // This is just for fixed cell-cycle models, need to work out how to find the g1 duration
        g1_duration = pCell->GetCellCycleModel()->GetTransitCellG1Duration();
    }

    if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        // Age of cell when apoptosis begins
        if (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime() < g1_duration)
        {
            cell_target_area *= 0.5*(1 + (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime())/g1_duration);
        }

        // The target area of an apoptotic cell decreases linearly to zero (and past it negative)
        cell_target_area = cell_target_area - 0.5*cell_target_area/(pCell->GetApoptosisTime())*(SimulationTime::Instance()->GetTime()-pCell->GetStartOfApoptosisTime());

        // Don't allow a negative target area
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else
    {
        // The target area of a proliferating cell increases linearly from A/2 to A over the course of the G1 phase
        if (cell_age < g1_duration)
        {
            cell_target_area *= 0.5*(1 + cell_age/g1_duration);
        }
        else
        {
            /**
             * At division, daughter cells inherit the cell data array from the mother cell.
             * Here, we assign the target area that we want daughter cells to have to cells
             * that we know to divide in this time step.
             *
             * \todo This is a little hack that we might want to clean up in the future.
             */
            if (pCell->ReadyToDivide())
            {
                cell_target_area = 0.5*mMatureCellTargetArea;
            }
        }
    }

    // Set cell data
    pCell->GetCellData()->SetItem("target area", cell_target_area);
}

template<unsigned DIM>
double TargetAreaGrowthModifier<DIM>::GetMatureCellTargetArea()
{
    return mMatureCellTargetArea;
}

template<unsigned DIM>
void TargetAreaGrowthModifier<DIM>::SetMatureCellTargetArea(double matureCellTargetArea)
{
    assert(matureCellTargetArea >= 0.0);
    mMatureCellTargetArea = matureCellTargetArea;
}

// Explicit instantiation
template class TargetAreaGrowthModifier<1>;
template class TargetAreaGrowthModifier<2>;
template class TargetAreaGrowthModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TargetAreaGrowthModifier)
