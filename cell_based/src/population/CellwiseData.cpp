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
//	if(mpCellPopulation != NULL)
//	{
//
//		for( typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
//				cell_iter != mpCellPopulation->End();
//				++cell_iter)
//		{
//			CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection();
//			if (collection.HasProperty<CellData>())
//			{
//				collection.RemoveProperty<CellData>();
//			}
//		}
//		MARK;
//	}
//	MARK;
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

    CellPropertyCollection cell_data_collection = pCell->rGetCellPropertyCollection().GetPropertiesType<CellData>();
	boost::shared_ptr<CellData> p_cell_data = boost::static_pointer_cast<CellData>(cell_data_collection.GetProperty());
	return p_cell_data->GetCellData(variableNumber);
}

template<unsigned DIM>
void CellwiseData<DIM>::SetValue(double value, unsigned locationIndex, unsigned variableNumber)
{
    assert(IsSetUp());
    if (variableNumber >= mNumberOfVariables)
    {
        EXCEPTION("Request for variable above mNumberOfVariables. Call SetNumCellsAndVars() to increase it.");
    }

    // Get the cell associated with locationIndex
    CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(locationIndex);

    CellPropertyCollection cell_data_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellData>();
    boost::shared_ptr<CellData> p_cell_data = boost::static_pointer_cast<CellData>(cell_data_collection.GetProperty());
    p_cell_data->SetCellData(variableNumber, value);

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
void CellwiseData<DIM>::SetNumCellsAndVars(unsigned numCells, unsigned numberOfVariables, AbstractCellPopulation<DIM>* pCellPopulation)
{
    if (mpCellPopulation!=NULL)
    {
        EXCEPTION("SetNumCellsAndVars() must be called before setting the CellPopulation (and after a Destroy)");
    }

    assert(pCellPopulation != NULL);
    assert(numCells == pCellPopulation->GetNumRealCells());
    assert(numberOfVariables > 0);
    assert(mAllocatedMemory == false);

    mNumberOfVariables = numberOfVariables;

    for( typename AbstractCellPopulation<DIM>::Iterator cell_iter = pCellPopulation->Begin();
    		cell_iter != pCellPopulation->End();
    		++cell_iter)
    {
//    	CellPropertyCollection cell_property_collection = (*cell_iter)->rGetCellPropertyCollection();
//    	CellPropertyCollection cell_data_collection = cell_property_collection.GetPropertiesType<CellData>();
//    	assert(cell_data_collection.GetSize() == 0);


    	CellPropertyCollection& collection = cell_iter->rGetCellPropertyCollection();
		if (collection.HasProperty<CellData>())
		{
			collection.RemoveProperty<CellData>();
		}
		assert(collection.GetPropertiesType<CellData>().GetSize() == 0);

    	MAKE_PTR_ARGS(CellData, p_cell_data, (numberOfVariables));
    	cell_iter->AddCellProperty(p_cell_data);
    }

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
