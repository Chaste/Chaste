/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "ChemotacticForce.hpp"

#include "CellwiseDataGradient.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
ChemotacticForce<DIM>::ChemotacticForce()
    : AbstractForce<DIM>()
{
}

template<unsigned DIM>
ChemotacticForce<DIM>::~ChemotacticForce()
{
}

template<unsigned DIM>
double ChemotacticForce<DIM>::GetChemotacticForceMagnitude(const double concentration, const double concentrationGradientMagnitude)
{
    return concentration; // temporary force law - can be changed to something realistic
                          // without tests failing
}

template<unsigned DIM>
void ChemotacticForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    CellwiseDataGradient<DIM> gradients;
    gradients.SetupGradients(rCellPopulation, "nutrient");

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Only labelled cells move chemotactically
        if (cell_iter->template HasCellProperty<CellLabel>())
        {
            unsigned node_global_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            c_vector<double,DIM>& r_gradient = gradients.rGetGradient(node_global_index);
            double nutrient_concentration = cell_iter->GetCellData()->GetItem("nutrient");
            double magnitude_of_gradient = norm_2(r_gradient);

            double force_magnitude = GetChemotacticForceMagnitude(nutrient_concentration, magnitude_of_gradient);

            // force += chi * gradC/|gradC|
            if (magnitude_of_gradient > 0)
            {
                c_vector<double,DIM> force = (force_magnitude/magnitude_of_gradient)*r_gradient;
                rCellPopulation.GetNode(node_global_index)->AddAppliedForceContribution(force);
            }
            // else Fc=0
        }
    }
}

template<unsigned DIM>
void ChemotacticForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ChemotacticForce<1>;
template class ChemotacticForce<2>;
template class ChemotacticForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemotacticForce)
