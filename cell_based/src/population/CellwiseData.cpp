/*

Copyright (C) University of Oxford, 2005-2012

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

#include "CellwiseData.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

template<unsigned DIM>
CellwiseData<DIM>* CellwiseData<DIM>::mpInstance = NULL;

template<unsigned DIM>
CellwiseData<DIM>* CellwiseData<DIM>::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CellwiseData<DIM>;
    }
    return mpInstance;
}

template<unsigned DIM>
CellwiseData<DIM>::CellwiseData()
    : mpCellPopulation(NULL),
      mAllocatedMemory(false),
      mNumberOfVariables(UNSIGNED_UNSET),
      mUseConstantDataForTesting(false)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}

template<unsigned DIM>
CellwiseData<DIM>::~CellwiseData()
{
}

template<unsigned DIM>
void CellwiseData<DIM>::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

template<unsigned DIM>
double CellwiseData<DIM>::GetValue(CellPtr pCell, unsigned variableNumber)
{
    if (variableNumber >= mNumberOfVariables)
    {
        EXCEPTION("Request for variable above mNumberOfVariables. Call SetNumCellsAndVars() to increase it.");
    }

    // To test a cell and cell-cycle models without a cell population
    if (mUseConstantDataForTesting)
    {
        return mConstantDataForTesting[variableNumber];
    }

    assert(IsSetUp());
    assert(mpCellPopulation != NULL);
    assert(mAllocatedMemory);

    unsigned location_index = mpCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned vector_index = location_index*mNumberOfVariables + variableNumber;
    return mData[vector_index];
}

template<unsigned DIM>
void CellwiseData<DIM>::SetValue(double value, unsigned locationIndex, unsigned variableNumber)
{
    assert(IsSetUp());
    if (variableNumber >= mNumberOfVariables)
    {
        EXCEPTION("Request for variable above mNumberOfVariables. Call SetNumCellsAndVars() to increase it.");
    }

    unsigned vector_index = locationIndex*mNumberOfVariables + variableNumber;
    assert(vector_index < mData.size());
    mData[vector_index] = value;
}

template<unsigned DIM>
void CellwiseData<DIM>::SetCellPopulation(AbstractCellPopulation<DIM>* pCellPopulation)
{
    if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(pCellPopulation))
    {
        EXCEPTION("CellwiseData does not work with ghost nodes.");
    }

    if (mAllocatedMemory == false)
    {
        EXCEPTION("SetCellPopulation must be called after SetNumCellsAndVars()");
    }

    mpCellPopulation = pCellPopulation;
}

template<unsigned DIM>
AbstractCellPopulation<DIM>& CellwiseData<DIM>::rGetCellPopulation()
{
    return *mpCellPopulation;
}

template<unsigned DIM>
void CellwiseData<DIM>::SetNumCellsAndVars(unsigned numCells, unsigned numberOfVariables)
{
    if (mpCellPopulation!=NULL)
    {
        EXCEPTION("SetNumCellsAndVars() must be called before setting the CellPopulation (and after a Destroy)");
    }

    assert(numberOfVariables > 0);
    assert(mAllocatedMemory == false);

    mNumberOfVariables = numberOfVariables;
    mData.clear();
    mData.resize(numCells * mNumberOfVariables, 0.0);

    mAllocatedMemory = true;
}

template<unsigned DIM>
bool CellwiseData<DIM>::IsSetUp()
{
    return ((mAllocatedMemory) && (mpInstance!=NULL) && (mpCellPopulation!=NULL));
}

template<unsigned DIM>
void CellwiseData<DIM>::ReallocateMemory()
{
    assert(mAllocatedMemory==true);
    assert(mpCellPopulation!=NULL);

    unsigned num_cells = mpCellPopulation->GetNumRealCells();

    if (mData.size() != num_cells*mNumberOfVariables)
    {
        mData.clear();
        mData.resize(num_cells*mNumberOfVariables, 0.0);
    }
}

template<unsigned DIM>
void CellwiseData<DIM>::SetConstantDataForTesting(std::vector<double>& rValues)
{
    mConstantDataForTesting = rValues;
    mUseConstantDataForTesting = true;
    mNumberOfVariables = 1;
}

template<unsigned DIM>
unsigned CellwiseData<DIM>::GetNumVariables()
{
    return mNumberOfVariables;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellwiseData<1>;
template class CellwiseData<2>;
template class CellwiseData<3>;
