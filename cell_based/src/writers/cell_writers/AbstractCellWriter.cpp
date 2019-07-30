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
#include "AbstractCellWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "UblasVectorInclude.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::AbstractCellWriter(const std::string& rFileName)
    : AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>(rFileName),
      mOutputScalarData(true),
      mOutputVectorData(false),
      mVtkCellDataName("DefaultVtkCellDataName"),
      mVtkVectorCellDataName("DefaultVtkVectorCellDataName")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputScalarData()
{
    return mOutputScalarData;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputVectorData()
{
    return mOutputVectorData;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::SetVtkCellDataName(std::string vtkCellDataName)
{
    mVtkCellDataName = vtkCellDataName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::SetVtkVectorCellDataName(std::string vtkCellDataName)
{
    mVtkVectorCellDataName = vtkCellDataName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(
        CellPtr pCell,
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return DOUBLE_UNSET;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorCellDataForVtkOutput(
        CellPtr pCell,
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return scalar_vector<double>(SPACE_DIM, DOUBLE_UNSET);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::GetVtkCellDataName()
{
    return mVtkCellDataName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::GetVtkVectorCellDataName()
{
    return mVtkVectorCellDataName;
}

// Explicit instantiation
template class AbstractCellWriter<1,1>;
template class AbstractCellWriter<1,2>;
template class AbstractCellWriter<2,2>;
template class AbstractCellWriter<1,3>;
template class AbstractCellWriter<2,3>;
template class AbstractCellWriter<3,3>;
