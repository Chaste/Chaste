/*

Copyright (c) 2005-2013, University of Oxford.
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
#include "CellProliferativeTypesCountWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::CellProliferativeTypesCountWriter(std::string directory)
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>(directory)
{
    this->mFileName = "celltypes.dat";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM>* pCellPopulation)
{
        this->WriteTimeStamp();

    std::vector<unsigned> proliferative_type_count = pCellPopulation->GetCellProliferativeTypeCount();

    for (unsigned i=0; i<proliferative_type_count.size(); i++)
    {
        *this->mpOutStream << proliferative_type_count[i] << "\t";
    }
    *this->mpOutStream << "\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>* pCellPopulation)
{
#define COVERAGE_IGNORE    //\ todo remove this when integrated with cell population.    #2183
    VisitAnyPopulation(pCellPopulation);
#undef COVERAGE_IGNORE
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MultipleCaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
#define COVERAGE_IGNORE    //\ todo remove this when integrated with cell population.    #2183
    VisitAnyPopulation(pCellPopulation);
#undef COVERAGE_IGNORE
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
#define COVERAGE_IGNORE    //\ todo remove this when integrated with cell population.    #2183
    VisitAnyPopulation(pCellPopulation);
#undef COVERAGE_IGNORE
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulationWithBuskeUpdate<SPACE_DIM>* pCellPopulation)
{
#define COVERAGE_IGNORE    //\ todo remove this when integrated with cell population.    #2183
    VisitAnyPopulation(pCellPopulation);
#undef COVERAGE_IGNORE
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulationWithParticles<SPACE_DIM>* pCellPopulation)
{
#define COVERAGE_IGNORE    //\ todo remove this when integrated with cell population.    #2183
    VisitAnyPopulation(pCellPopulation);
#undef COVERAGE_IGNORE
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
#define COVERAGE_IGNORE    //\ todo remove this when integrated with cell population.    #2183
    VisitAnyPopulation(pCellPopulation);
#undef COVERAGE_IGNORE
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
#define COVERAGE_IGNORE    //\ todo remove this when integrated with cell population.    #2183
    VisitAnyPopulation(pCellPopulation);
#undef COVERAGE_IGNORE
}

// Explicit instantiation
template class CellProliferativeTypesCountWriter<1,1>;
template class CellProliferativeTypesCountWriter<2,2>;
template class CellProliferativeTypesCountWriter<3,3>;
