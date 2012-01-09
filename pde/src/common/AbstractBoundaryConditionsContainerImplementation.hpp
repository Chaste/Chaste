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
            mpDirichletMap[i] = NULL;
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
