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

#ifndef ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_
#define ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <set>
#include <boost/utility.hpp>
#include "AbstractBoundaryCondition.hpp"
#include "Node.hpp"
//#include "ConstBoundaryCondition.hpp"
//#include "TetrahedralMesh.hpp"
//#include "LinearSystem.hpp"
//#include "PetscException.hpp"

/**
 * Helper struct storing an operator for computing whether one node
 * has a lower index than another.
 */
template<unsigned SPACE_DIM>
struct LessThanNode
{
    /**
     * Less-then node index comparison operator.
     * @return true if n1<n2
     * @param n1 pointer to a node
     * @param n2 pointer to a node
     */
    bool operator()(const Node<SPACE_DIM> * const &n1, const Node<SPACE_DIM> * const &n2) const
    {
        return (n1->GetIndex() < n2->GetIndex() );
    }
};

/**
 * Abstract boundary conditions container.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractBoundaryConditionsContainer : boost::noncopyable
{
protected:

    /** To save typing */
    typedef typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >
        DirichletMapType;
    DirichletMapType* mpDirichletMap[PROBLEM_DIM]; /**< List (map) of Dirichlet boundary conditions */

    /** To save typing */
    typedef typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >::const_iterator
        DirichletIteratorType;
    DirichletIteratorType mDirichIterator; /**< Internal iterator over Dirichlet boundary conditions */

    /** @return true if there are any Dirichlet BCs anywhere on the mesh*/
    bool mHasDirichletBCs;

    /** @return true if we calculated mHasDirichletBCs. */
    bool mCheckedAndCommunicatedIfDirichletBcs;

    /** @return true if we need to delete BCs in destructor. */
    bool mDeleteConditions;

    /**
     * Delete the list of Dirichlet boundary conditions.
     *
     * @note This should stay as a protected method to avoid it being called with default arguments and causing seg faults
     *  (requires careful bookkeeping when calling this method).
     * @param alreadyDeletedConditions  This is a set of BCs that have already been deleted that we should avoid trying
     *  to delete inside this method. (defaults to empty = delete everything)
     */
    void DeleteDirichletBoundaryConditions(std::set<const AbstractBoundaryCondition<SPACE_DIM>*> alreadyDeletedConditions
                                            = std::set<const AbstractBoundaryCondition<SPACE_DIM>*>());

public:

    /**
     * Constructor allocates memory for the Dirichlet boundary conditions lists.
     *
     * @param deleteConditions whether to delete BCs in destructor (defaults to true)
     */
    AbstractBoundaryConditionsContainer(bool deleteConditions=true);

    /**
     * Destructor.
     */
    ~AbstractBoundaryConditionsContainer();

    /**
     * @return whether any Dirichlet conditions are defined (for ANY of the unknowns, on ANY of the processes).
     * Must be called collectively. The first time this is called, the result is communicated to all processes
     * and then cached locally (the bool mHasDirichletBCs). If this needs recalculating
     * mCheckedAndCommunicatedIfDirichletBcs must be reset to zero.
     */
    bool HasDirichletBoundaryConditions();

    /**
     * @return value of Dirichlet boundary condition at specified node.
     *
     * This is unlikely to be needed by the user, the methods ApplyDirichletToLinearProblem or
     * ApplyDirichletToNonlinearProblem can be called instead to apply all Dirichlet boundary conditions
     * at the same time.
     *
     * @param pBoundaryNode pointer to a boundary node
     * @param indexOfUnknown index of the unknown for which to obtain the value of the boundary condition (defaults to 0)
     */
    double GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode, unsigned indexOfUnknown = 0);

    /**
     * @return true if there is a Dirichlet boundary condition defined on the given node.
     *
     * @param pNode pointer to a node
     * @param indexOfUnknown index of the unknown for which to obtain the value of the boundary condition (defaults to 0)
     */
    bool HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode, unsigned indexOfUnknown = 0);

    /**
     * When Dirichlet boundary conditions are likely to be added on one or more processes then we should call this
     * method collectively in order to ensure that all processes do a collective communication on the next call
     * to HasDirichletBoundaryConditions()
     */
    void ResetDirichletCommunication()
    {
        mCheckedAndCommunicatedIfDirichletBcs = false;
    }
};

#endif /*ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_*/
