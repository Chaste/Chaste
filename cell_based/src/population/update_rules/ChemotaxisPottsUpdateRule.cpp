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

#include "ChemotaxisPottsUpdateRule.hpp"

template<unsigned DIM>
ChemotaxisPottsUpdateRule<DIM>::ChemotaxisPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>()
{
}

template<unsigned DIM>
ChemotaxisPottsUpdateRule<DIM>::~ChemotaxisPottsUpdateRule()
{
}

template<unsigned DIM>
double ChemotaxisPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    double delta_H = 0.0;
    c_vector<double, DIM> current_location = rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
    c_vector<double, DIM> target_location = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation();

    for (unsigned dimension = 0; dimension < DIM; dimension++)
    {
        if(target_location[dimension] > current_location[dimension])
        {
            delta_H -= 0.2;
        }
        else if(target_location[dimension] < current_location[dimension])
        {
            delta_H += 0.2;
        }
    }

    return delta_H;
}

template<unsigned DIM>
void ChemotaxisPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ChemotaxisPottsUpdateRule<1>;
template class ChemotaxisPottsUpdateRule<2>;
template class ChemotaxisPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChemotaxisPottsUpdateRule)
