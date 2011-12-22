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

#include "AbstractCellPopulationBoundaryCondition.hpp"

template<unsigned DIM>
AbstractCellPopulationBoundaryCondition<DIM>::AbstractCellPopulationBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation)
        : mpCellPopulation(pCellPopulation)
{
}

template<unsigned DIM>
AbstractCellPopulationBoundaryCondition<DIM>::~AbstractCellPopulationBoundaryCondition()
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM>* AbstractCellPopulationBoundaryCondition<DIM>::GetCellPopulation() const
{
    return mpCellPopulation;
}

template<unsigned DIM>
void AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionInfo(out_stream& rParamsFile)
{
    std::string cell_killer_type = GetIdentifier();

    *rParamsFile << "\t\t<" << cell_killer_type << ">\n";
    OutputCellPopulationBoundaryConditionParameters(rParamsFile);
    *rParamsFile << "\t\t</" << cell_killer_type << ">\n";
}

template<unsigned DIM>
void AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    // No parameters to output
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AbstractCellPopulationBoundaryCondition<1>;
template class AbstractCellPopulationBoundaryCondition<2>;
template class AbstractCellPopulationBoundaryCondition<3>;
