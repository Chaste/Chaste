/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "ParabolicPdeAndBoundaryConditions.hpp"
#include "AveragedSourceParabolicPde.hpp"

template<unsigned DIM>
ParabolicPdeAndBoundaryConditions<DIM>::ParabolicPdeAndBoundaryConditions(AbstractLinearParabolicPde<DIM,DIM>* pPde,
                                                        AbstractBoundaryCondition<DIM>* pBoundaryCondition,
                                                        bool isNeumannBoundaryCondition,
                                                        Vec solution,
                                                        bool deleteMemberPointersInDestructor)
    : mpPde(pPde),
      mpBoundaryCondition(pBoundaryCondition),
      mIsNeumannBoundaryCondition(isNeumannBoundaryCondition),
      mSolution(NULL),
      mDeleteMemberPointersInDestructor(deleteMemberPointersInDestructor),
      mDependentVariableName("")
{
    if (solution)
    {
        mSolution = solution;
    }
}

template<unsigned DIM>
ParabolicPdeAndBoundaryConditions<DIM>::~ParabolicPdeAndBoundaryConditions()
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
AbstractLinearParabolicPde<DIM,DIM>* ParabolicPdeAndBoundaryConditions<DIM>::GetPde()
{
    return mpPde;
}

template<unsigned DIM>
AbstractBoundaryCondition<DIM>* ParabolicPdeAndBoundaryConditions<DIM>::GetBoundaryCondition() const
{
    return mpBoundaryCondition;
}

template<unsigned DIM>
Vec ParabolicPdeAndBoundaryConditions<DIM>::GetSolution()
{
    return mSolution;
}

template<unsigned DIM>
Vec ParabolicPdeAndBoundaryConditions<DIM>::GetSolution() const
{
    return mSolution;
}

template<unsigned DIM>
void ParabolicPdeAndBoundaryConditions<DIM>::SetSolution(Vec solution)
{
    mSolution = solution;
}

template<unsigned DIM>
bool ParabolicPdeAndBoundaryConditions<DIM>::IsNeumannBoundaryCondition()
{
    return mIsNeumannBoundaryCondition;
}

template<unsigned DIM>
bool ParabolicPdeAndBoundaryConditions<DIM>::HasAveragedSourcePde()
{
    return (dynamic_cast<AveragedSourceParabolicPde<DIM>*>(mpPde) != NULL);
}

template<unsigned DIM>
void ParabolicPdeAndBoundaryConditions<DIM>::DestroySolution()
{
    if (mSolution)
    {
        PetscTools::Destroy(mSolution);
    }
}

template<unsigned DIM>
void ParabolicPdeAndBoundaryConditions<DIM>::SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM,DIM>* pMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap)
{
    assert(HasAveragedSourcePde());
    static_cast<AveragedSourceParabolicPde<DIM>*>(mpPde)->SetupSourceTerms(*pMesh, pCellPdeElementMap);
}

template<unsigned DIM>
void ParabolicPdeAndBoundaryConditions<DIM>::SetDependentVariableName(const std::string& rName)
{
    mDependentVariableName = rName;
}

template<unsigned DIM>
std::string& ParabolicPdeAndBoundaryConditions<DIM>::rGetDependentVariableName()
{
    return mDependentVariableName;
}

// Explicit instantiation
template class ParabolicPdeAndBoundaryConditions<1>;
template class ParabolicPdeAndBoundaryConditions<2>;
template class ParabolicPdeAndBoundaryConditions<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParabolicPdeAndBoundaryConditions)
