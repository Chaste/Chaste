/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "FarhadifarTypeModifier.hpp"

template<unsigned DIM>
FarhadifarTypeModifier<DIM>::FarhadifarTypeModifier()
    : AbstractTargetAreaModifier<DIM>()
{
}

template<unsigned DIM>
FarhadifarTypeModifier<DIM>::~FarhadifarTypeModifier()
{
}

template<unsigned DIM>
void FarhadifarTypeModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = this->mReferenceTargetArea;

    if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        // The target area of an apoptotic cell decreases linearly to zero (and past it negative)
        ///\todo: which cells are apoptotic? if they get apoptotic during G2-phase then this line has to be changed
        cell_target_area = cell_target_area - 0.5*cell_target_area/(pCell->GetApoptosisTime())*(SimulationTime::Instance()->GetTime()-pCell->GetStartOfApoptosisTime());

        // Don't allow a negative target area
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else
    {
        double cell_age = pCell->GetAge();

        // Get the combined duration of the cell's M, G1 and S phases
        AbstractCellCycleModel* p_model = pCell->GetCellCycleModel();
        double growth_start_time = p_model->GetMDuration() + p_model->GetG1Duration() + p_model->GetSDuration();

        // The target area of a proliferating cell increases linearly from A to 2A over the course of the G2 phase
        if (cell_age > growth_start_time)
        {
            double g2_duration = p_model->GetG2Duration();
            cell_target_area *= (1 + (cell_age-growth_start_time)/g2_duration);
        }

        /**
         * At division, daughter cells inherit the cell data array from the mother cell.
         * Here, we assign the target area that we want daughter cells to have to cells
         * that we know to divide in this time step.
         *
         * \todo This is a little hack that we might want to clean up in the future.
         */
        if (pCell->ReadyToDivide())
        {
            cell_target_area = this->mReferenceTargetArea;
        }
    }

    // Set cell data
    pCell->GetCellData()->SetItem("target area", cell_target_area);
}

template<unsigned DIM>
void FarhadifarTypeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractTargetAreaModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class FarhadifarTypeModifier<1>;
template class FarhadifarTypeModifier<2>;
template class FarhadifarTypeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FarhadifarTypeModifier)
