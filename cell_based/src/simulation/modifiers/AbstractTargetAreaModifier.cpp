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

#include "AbstractTargetAreaModifier.hpp"

template<unsigned DIM>
AbstractTargetAreaModifier<DIM>::AbstractTargetAreaModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mReferenceTargetArea(1.0)
{
}

template<unsigned DIM>
AbstractTargetAreaModifier<DIM>::~AbstractTargetAreaModifier()
{
}

template<unsigned DIM>
void AbstractTargetAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateTargetAreas(rCellPopulation);
}

template<unsigned DIM>
void AbstractTargetAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateTargetAreas(rCellPopulation);
}

template<unsigned DIM>
void AbstractTargetAreaModifier<DIM>::UpdateTargetAreas(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Loop over the list of cells, rather than using the population iterator, so as to include dead cells
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
         cell_iter != rCellPopulation.rGetCells().end();
         ++cell_iter)
    {
        UpdateTargetAreaOfCell(*cell_iter);
    }
}

template<unsigned DIM>
double AbstractTargetAreaModifier<DIM>::GetReferenceTargetArea()
{
    return mReferenceTargetArea;
}

template<unsigned DIM>
void AbstractTargetAreaModifier<DIM>::SetReferenceTargetArea(double referenceTargetArea)
{
    assert(referenceTargetArea >= 0.0);
    mReferenceTargetArea = referenceTargetArea;
}

template<unsigned DIM>
void AbstractTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ReferenceTargetArea>" << mReferenceTargetArea << "</ReferenceTargetArea>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractTargetAreaModifier<1>;
template class AbstractTargetAreaModifier<2>;
template class AbstractTargetAreaModifier<3>;
