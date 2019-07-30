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

#ifndef ABSTRACTBOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_
#define ABSTRACTBOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_

#include "AbstractBoundaryConditionsContainer.hpp"
#include "PetscTools.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractBoundaryConditionsContainer(bool deleteConditions)
    : mHasDirichletBCs(false),
      mCheckedAndCommunicatedIfDirichletBcs(false),
      mDeleteConditions(deleteConditions)
{
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        mpDirichletMap[index_of_unknown] = new std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~AbstractBoundaryConditionsContainer()
{
    DeleteDirichletBoundaryConditions();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::HasDirichletBoundaryConditions()
{
    if (!mCheckedAndCommunicatedIfDirichletBcs)
    {
        bool i_have_dirichlet=false;
        for (unsigned i=0; i<PROBLEM_DIM; i++)
        {
            if (!mpDirichletMap[i]->empty())
            {
                i_have_dirichlet=true;
                break;
            }
        }
        mHasDirichletBCs = PetscTools::ReplicateBool(i_have_dirichlet);
        mCheckedAndCommunicatedIfDirichletBcs = true;
    }
    return mHasDirichletBCs;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DeleteDirichletBoundaryConditions(std::set<const AbstractBoundaryCondition<SPACE_DIM>*> alreadyDeletedConditions)
{
    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        if (mpDirichletMap[i])
        {
            mDirichIterator = mpDirichletMap[i]->begin();
            while (mDirichIterator != mpDirichletMap[i]->end() )
            {
                if (alreadyDeletedConditions.count(mDirichIterator->second) == 0)
                {
                    alreadyDeletedConditions.insert(mDirichIterator->second);
                    if (mDeleteConditions)
                    {
                        delete mDirichIterator->second;
                    }
                }
                mDirichIterator++;
            }

            delete(mpDirichletMap[i]);
            mpDirichletMap[i] = nullptr;
        }
    }

    // Recommunicate that Dirichlet BCs have changed (next time we ask)
    ResetDirichletCommunication();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode, unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);
    //assert(pBoundaryNode->IsBoundaryNode());

    mDirichIterator = mpDirichletMap[indexOfUnknown]->find(pBoundaryNode);
    assert(mDirichIterator != mpDirichletMap[indexOfUnknown]->end());

    return mDirichIterator->second->GetValue(pBoundaryNode->GetPoint());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode, unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    this->mDirichIterator = this->mpDirichletMap[indexOfUnknown]->find(pNode);

    return (this->mDirichIterator != this->mpDirichletMap[indexOfUnknown]->end());
}

#endif /*ABSTRACTBOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_*/
