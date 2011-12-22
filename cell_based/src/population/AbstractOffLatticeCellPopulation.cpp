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

#include "AbstractOffLatticeCellPopulation.hpp"

template<unsigned DIM>
AbstractOffLatticeCellPopulation<DIM>::AbstractOffLatticeCellPopulation(std::vector<CellPtr>& rCells,
                                                                  const std::vector<unsigned> locationIndices)
    : AbstractCellPopulation<DIM>(rCells, locationIndices),
      mDampingConstantNormal(1.0),
      mDampingConstantMutant(1.0)
{
}

template<unsigned DIM>
AbstractOffLatticeCellPopulation<DIM>::AbstractOffLatticeCellPopulation()
    : AbstractCellPopulation<DIM>()
{
}

template<unsigned DIM>
void AbstractOffLatticeCellPopulation<DIM>::SetDampingConstantNormal(double dampingConstantNormal)
{
    assert(dampingConstantNormal > 0.0);
    mDampingConstantNormal = dampingConstantNormal;
}

template<unsigned DIM>
void AbstractOffLatticeCellPopulation<DIM>::SetDampingConstantMutant(double dampingConstantMutant)
{
    assert(dampingConstantMutant > 0.0);
    mDampingConstantMutant = dampingConstantMutant;
}

template<unsigned DIM>
double AbstractOffLatticeCellPopulation<DIM>::GetDampingConstantNormal()
{
    return mDampingConstantNormal;
}

template<unsigned DIM>
double AbstractOffLatticeCellPopulation<DIM>::GetDampingConstantMutant()
{
    return mDampingConstantMutant;
}

template<unsigned DIM>
void AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<DampingConstantNormal>" << mDampingConstantNormal << "</DampingConstantNormal>\n";
    *rParamsFile << "\t\t<DampingConstantMutant>" << mDampingConstantMutant << "</DampingConstantMutant>\n";

    // Call method on direct parent class
    AbstractCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractOffLatticeCellPopulation<1>;
template class AbstractOffLatticeCellPopulation<2>;
template class AbstractOffLatticeCellPopulation<3>;
