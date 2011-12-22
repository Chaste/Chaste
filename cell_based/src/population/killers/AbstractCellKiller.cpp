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

#include "AbstractCellKiller.hpp"

template<unsigned SPACE_DIM>
AbstractCellKiller<SPACE_DIM>::AbstractCellKiller(AbstractCellPopulation<SPACE_DIM>* pCellPopulation)
        : mpCellPopulation(pCellPopulation)
{
}

template<unsigned SPACE_DIM>
AbstractCellKiller<SPACE_DIM>::~AbstractCellKiller()
{
}

template<unsigned SPACE_DIM>
const AbstractCellPopulation<SPACE_DIM>* AbstractCellKiller<SPACE_DIM>::GetCellPopulation() const
{
    return mpCellPopulation;
}

template<unsigned DIM>
void AbstractCellKiller<DIM>::OutputCellKillerInfo(out_stream& rParamsFile)
{
    std::string cell_killer_type = GetIdentifier();

    *rParamsFile << "\t\t<" << cell_killer_type << ">\n";
    OutputCellKillerParameters(rParamsFile);
    *rParamsFile << "\t\t</" << cell_killer_type << ">\n";
}

template<unsigned DIM>
void AbstractCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // No parameters to output
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AbstractCellKiller<1>;
template class AbstractCellKiller<2>;
template class AbstractCellKiller<3>;
