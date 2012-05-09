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

#include "CellwiseData.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellPropertyCollection.hpp"

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
      mAllocatedMemory(true),
      mNumberOfVariables(0u), ///\todo This is a big positive number - It should be zero so that we can comparison test against
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
void CellwiseData<DIM>::SetValue(double value, unsigned locationIndex, unsigned variableNumber)
{
    assert(IsSetUp());

    if (variableNumber >= mNumberOfVariables)
    {
        EXCEPTION("Request for variable above mNumberOfVariables.");
    }

    // Get the cell associated with locationIndex
    CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(locationIndex);

	p_cell->GetCellData()->SetItem(variableNumber, value);
}

template<unsigned DIM>
AbstractCellPopulation<DIM>& CellwiseData<DIM>::rGetCellPopulation()
{
    return *mpCellPopulation;
}

template<unsigned DIM>
void CellwiseData<DIM>::SetPopulationAndNumVars(AbstractCellPopulation<DIM>* pCellPopulation, unsigned numberOfVariables)
{
    assert(pCellPopulation != NULL);
    assert(numberOfVariables > 0);

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(pCellPopulation))
	{
		EXCEPTION("CellwiseData does not work with ghost nodes.");
	}

	if (IsSetUp())
	{
		EXCEPTION("Can't call SetPopulationAndNumVars() once CellwiseData is setup.");
	}

	mpCellPopulation = pCellPopulation;
    mNumberOfVariables = numberOfVariables;

    for( typename AbstractCellPopulation<DIM>::Iterator cell_iter = pCellPopulation->Begin();
    		cell_iter != pCellPopulation->End();
    		++cell_iter)
    {
    	CellPropertyCollection& collection = cell_iter->rGetCellPropertyCollection();
		if (collection.HasProperty<CellData>())
		{
			collection.RemoveProperty<CellData>();
		}
		assert(collection.GetPropertiesType<CellData>().GetSize() == 0);

    	MAKE_PTR_ARGS(CellData, p_cell_data, (numberOfVariables));
    	cell_iter->AddCellProperty(p_cell_data);
    }

}

template<unsigned DIM>
bool CellwiseData<DIM>::IsSetUp()
{
    return ((mpInstance!=NULL) && (mpCellPopulation!=NULL));
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
