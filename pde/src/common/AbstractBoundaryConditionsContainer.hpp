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

#ifndef ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_
#define ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <set>
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
     *
     * @param n1 pointer to a node
     * @param n2 pointer to a node
     */
    bool operator()(const Node<SPACE_DIM> * const &n1, const Node<SPACE_DIM> * const &n2)
    {
        return (n1->GetIndex() < n2->GetIndex() );
    }
};

/**
 * Abstract boundary conditions container.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractBoundaryConditionsContainer
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

    /** Whether there are any Dirichlet BCs anywhere on the mesh*/
    bool mHasDirichletBCs;

    /** Have we calculated mHasDirichletBCs. */
    bool mCheckedAndCommunicatedIfDirichletBcs;

    /** Whether to delete BCs in destructor. */
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
     * Return whether any Dirichlet conditions are defined (for ANY of the unknowns, on ANY of the processes).
     * Must be called collectively. The first time this is called, the result is communicated to all processes
     * and then cached locally (the bool mHasDirichletBCs). If this needs recalculating
     * mCheckedAndCommunicatedIfDirichletBcs must be reset to zero.
     */
    bool HasDirichletBoundaryConditions();

    /**
     * Obtain value of Dirichlet boundary condition at specified node.
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
     * Test if there is a Dirichlet boundary condition defined on the given node.
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
