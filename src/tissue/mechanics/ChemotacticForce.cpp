/*

Copyright (C) University of Oxford, 2005-2010

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


#include "ChemotacticForce.hpp"
#include "CellwiseDataGradient.hpp"

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
void ChemotacticForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                 AbstractTissue<DIM>& rTissue)
{
    CellwiseDataGradient<DIM> gradients;
    gradients.SetupGradients();

    for (typename AbstractTissue<DIM>::Iterator cell_iter = rTissue.Begin();
         cell_iter != rTissue.End();
         ++cell_iter)
    {
        // Only labelled cells move chemotactically
        boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();
        if (p_state->IsType<LabelledCellMutationState>())
        {
            unsigned node_global_index = rTissue.GetLocationIndexUsingCell(*cell_iter);

            c_vector<double,DIM>& r_gradient = gradients.rGetGradient(node_global_index);
            double nutrient_concentration = CellwiseData<DIM>::Instance()->GetValue(*cell_iter, 0);
            double magnitude_of_gradient = norm_2(r_gradient);

            double force_magnitude = GetChemotacticForceMagnitude(nutrient_concentration, magnitude_of_gradient);

            // force += chi * gradC/|gradC|
            if (magnitude_of_gradient > 0)
            {
                rForces[node_global_index] += (force_magnitude/magnitude_of_gradient)*r_gradient;
            }
            // else Fc=0
        }
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ChemotacticForce<1>;
template class ChemotacticForce<2>;
template class ChemotacticForce<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemotacticForce)
