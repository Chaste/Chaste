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

#include "AbstractOffLatticeCellPopulation.hpp"

template<unsigned DIM>
AbstractOffLatticeCellPopulation<DIM>::AbstractOffLatticeCellPopulation( AbstractMesh<DIM, DIM>& rMesh,
																	std::vector<CellPtr>& rCells,
                                                                  const std::vector<unsigned> locationIndices)
    : AbstractCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDampingConstantNormal(1.0),
      mDampingConstantMutant(1.0),
      mAbsoluteMovementThreshold(0.5)
{
}

template<unsigned DIM>
AbstractOffLatticeCellPopulation<DIM>::AbstractOffLatticeCellPopulation(AbstractMesh<DIM, DIM>& rMesh)
	: AbstractCellPopulation<DIM>(rMesh)
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
void AbstractOffLatticeCellPopulation<DIM>::SetAbsoluteMovementThreshold(double absoluteMovementThreshold)
{
    mAbsoluteMovementThreshold = absoluteMovementThreshold;
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
double AbstractOffLatticeCellPopulation<DIM>::GetAbsoluteMovementThreshold()
{
    return mAbsoluteMovementThreshold;
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
