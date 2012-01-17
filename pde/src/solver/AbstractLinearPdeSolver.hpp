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

#ifndef ABSTRACTLINEARPDESOLVER_HPP_
#define ABSTRACTLINEARPDESOLVER_HPP_

#include "LinearSystem.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "HeartEventHandler.hpp"

/**
 * Simple abstract class containing some common functionality between
 * AbstractStaticLinearPdeSolver and AbstractDynamicLinearPdeSolver.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractLinearPdeSolver : boost::noncopyable
{
protected:

    /** The linear system that will be set up and solved as part of the PDE solve. */
    LinearSystem* mpLinearSystem;

    /** Pointer to the mesh. */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

public:

    /**
     * Constructor.
     *
     * @param pMesh the mesh
     */
    AbstractLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        : mpLinearSystem(NULL),
          mpMesh(pMesh)
    {
        assert(pMesh!=NULL);
    }

    /**
     * Destructor.
     */
    virtual ~AbstractLinearPdeSolver()
    {
        if (mpLinearSystem)
        {
            delete mpLinearSystem;
        }
    }

    /**
     * Initialise method: sets up the linear system (using the mesh to
     * determine the number of unknowns per row to preallocate) if it
     * is not already set up. Can use an initial solution as PETSc template,
     * or base it on the mesh size.
     *
     * @param initialSolution Initial solution (defaults to NULL) for PETSc
     *  to use as a template.
     */
    virtual void InitialiseForSolve(Vec initialSolution = NULL);

    /**
     * The static and dynamic Solve() implementations both call this
     * before after SetupLinearSystem(). It can be overloaded if needed.
     *
     * @param currentSolution The current solution
     */
    virtual void PrepareForSetupLinearSystem(Vec currentSolution)
    {
    }

    /**
     * The static and dynamic Solve() implementations both call this
     * immediately after SetupLinearSystem(). It can be overloaded if
     * further work needs to be done.
     *
     * @param currentSolution The current solution
     */
    virtual void FinaliseLinearSystem(Vec currentSolution)
    {
    }

    /**
     * The static and dynamic Solve() implementations both call this immediately after
     * the linear solve is carried out (but before the timestep counter is incremented.
     * This can be overloaded if further work on the solution vector needs to be done
     * (for example, in operator splitting of the diffusion and reaction terms in the
     * OperatorSplittingMonodomainSolver.
     *
     * @param currentSolution The current solution (solution of the linear system solve)
     */
    virtual void FollowingSolveLinearSystem(Vec currentSolution)
    {
    }

    /**
     * The main Solve() methods in the child classes use this method. The concrete
     * solver classes must implement it, depending on the the choice of numerical
     * approach. The method should completely set up the linear system that has to
     * be solved (that timestep, if dynamic PDEs).
     *
     * @param currentSolution The current solution which can be used in setting up
     *  the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *  (mainly for dynamic solves).
     */
    virtual void SetupLinearSystem(Vec currentSolution, bool computeMatrix)=0;

    /**
     * Get a pointer to the linear system.
     */
    LinearSystem* GetLinearSystem()
    {
        return mpLinearSystem;
    }
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem == NULL)
    {
        unsigned preallocation=(mpMesh->CalculateMaximumContainingElementsPerProcess() + ELEMENT_DIM);
        if (ELEMENT_DIM > 1)
        {
            // Highest connectivity is closed
            preallocation--;
        }
        preallocation *= PROBLEM_DIM;

        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        if (initialSolution == NULL)
        {
            /*
             * Static problem, create linear system. The following ensures
             * all the unknowns for a particular node are on the same processor.
             */
            Vec template_vec = mpMesh->GetDistributedVectorFactory()->CreateVec(PROBLEM_DIM);

            this->mpLinearSystem = new LinearSystem(template_vec, preallocation);

            PetscTools::Destroy(template_vec);
        }
        else
        {
            /*
             * Use the current solution (ie the initial solution)
             * as the template in the alternative constructor of
             * LinearSystem. This is to avoid problems with VecScatter.
             */
            this->mpLinearSystem = new LinearSystem(initialSolution, preallocation);
        }

        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    }
}

#endif /*ABSTRACTLINEARPDESOLVER_HPP_*/
