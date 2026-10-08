/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "CellPopulationExtentWriter.hpp"

#include "AbstractCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "SemBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::CellPopulationExtentWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("cellpopulationextent.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    WriteExtent(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    WriteExtent(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    WriteExtent(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    WriteExtent(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    WriteExtent(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    ImmersedBoundaryCellPopulation<SPACE_DIM>* pCellPopulation)
{
    WriteExtent(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationExtentWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    SemBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    WriteExtent(pCellPopulation);
}

// Explicit instantiation
template class CellPopulationExtentWriter<1, 1>;
template class CellPopulationExtentWriter<1, 2>;
template class CellPopulationExtentWriter<2, 2>;
template class CellPopulationExtentWriter<1, 3>;
template class CellPopulationExtentWriter<2, 3>;
template class CellPopulationExtentWriter<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPopulationExtentWriter)
