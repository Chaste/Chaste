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
class AbstractLinearPdeSolver : private boost::noncopyable
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
        : mpLinearSystem(nullptr),
          mpMesh(pMesh)
    {
        assert(pMesh!=nullptr);
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
    virtual void InitialiseForSolve(Vec initialSolution = nullptr);

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
     * @return a pointer to the linear system.
     */
    LinearSystem* GetLinearSystem()
    {
        return mpLinearSystem;
    }
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem == nullptr)
    {
        unsigned preallocation = PROBLEM_DIM * mpMesh->CalculateMaximumNodeConnectivityPerProcess();

        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        if (initialSolution == nullptr)
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
