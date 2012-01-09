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

#include "PdeAndBoundaryConditions.hpp"

template<unsigned DIM>
PdeAndBoundaryConditions<DIM>::PdeAndBoundaryConditions(AbstractLinearEllipticPde<DIM,DIM>* pPde,
                                                        AbstractBoundaryCondition<DIM>* pBoundaryCondition,
                                                        bool isNeumannBoundaryCondition,
                                                        Vec solution,
                                                        bool deleteMemberPointersInDestructor)
    : mpPde(pPde),
      mpBoundaryCondition(pBoundaryCondition),
      mIsNeumannBoundaryCondition(isNeumannBoundaryCondition),
      mSolution(NULL),
      mDeleteMemberPointersInDestructor(deleteMemberPointersInDestructor)
{
    if (solution)
    {
        mSolution = solution;
    }
}

template<unsigned DIM>
PdeAndBoundaryConditions<DIM>::~PdeAndBoundaryConditions()
{
    // Avoid memory leaks if the object was loaded from an archive
    if (mDeleteMemberPointersInDestructor)
    {
        delete mpPde;
        delete mpBoundaryCondition;
    }

    DestroySolution();
}

template<unsigned DIM>
AbstractLinearEllipticPde<DIM,DIM>* PdeAndBoundaryConditions<DIM>::GetPde()
{
    return mpPde;
}

template<unsigned DIM>
AbstractBoundaryCondition<DIM>* PdeAndBoundaryConditions<DIM>::GetBoundaryCondition() const
{
    return mpBoundaryCondition;
}

template<unsigned DIM>
Vec PdeAndBoundaryConditions<DIM>::GetSolution()
{
    return mSolution;
}

template<unsigned DIM>
Vec PdeAndBoundaryConditions<DIM>::GetSolution() const
{
    return mSolution;
}

template<unsigned DIM>
void PdeAndBoundaryConditions<DIM>::SetSolution(Vec solution)
{
    mSolution = solution;
}

template<unsigned DIM>
bool PdeAndBoundaryConditions<DIM>::IsNeumannBoundaryCondition()
{
    return mIsNeumannBoundaryCondition;
}

template<unsigned DIM>
bool PdeAndBoundaryConditions<DIM>::HasAveragedSourcePde()
{
    return (dynamic_cast<AveragedSourcePde<DIM>*>(mpPde) != NULL);
}

template<unsigned DIM>
void PdeAndBoundaryConditions<DIM>::DestroySolution()
{
    if (mSolution)
    {
        VecDestroy(mSolution);
    }
}

template<unsigned DIM>
void PdeAndBoundaryConditions<DIM>::SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM,DIM>* pMesh, std::map< CellPtr, unsigned >* pCellPdeElementMap)
{
    assert(HasAveragedSourcePde());
    static_cast<AveragedSourcePde<DIM>*>(mpPde)->SetupSourceTerms(*pMesh, pCellPdeElementMap);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class PdeAndBoundaryConditions<1>;
template class PdeAndBoundaryConditions<2>;
template class PdeAndBoundaryConditions<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PdeAndBoundaryConditions)
