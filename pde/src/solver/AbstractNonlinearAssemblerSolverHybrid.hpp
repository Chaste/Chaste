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

#ifndef ABSTRACTNONLINEARASSEMBLERSOLVERHYBRID_HPP_
#define ABSTRACTNONLINEARASSEMBLERSOLVERHYBRID_HPP_

#include "BoundaryConditionsContainer.hpp"
#include "AbstractNonlinearEllipticPde.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "PetscTools.hpp"
#include "PetscMatTools.hpp"
#include "PetscVecTools.hpp"
#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "LinearSystem.hpp"
#include "SimplePetscNonlinearSolver.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"

/*
 * Definitions of GLOBAL functions used by PETSc nonlinear solver
 * (implementations are at the bottom of this file).
 */

/**
 * Function called by PETSc to compute the residual vector given the current solution guess.
 * Calls a method on a AbstractNonlinearAssemblerSolverHybrid object to do the calculation.
 *
 * @param snes This is not used by us, but required by PETSc.
 * @param currentGuess The solution guess for the current iteration.
 * @param residualVector We fill this with the residual vector.
 * @param pContext Pointer to a AbstractNonlinearAssemblerSolverHybrid object.
 * @return PETSc expects this function to return a PetscErrorCode. We always return 0;
 *   exceptions are thrown if there is an error.
 *
 * Note: this is a global function, hence the need for a long name to avoid
 * potential conflicting names
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeResidual(SNES snes,
                                                                      Vec currentGuess,
                                                                      Vec residualVector,
                                                                      void* pContext);
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
    /**
     * Function called by PETSc to compute the Jacobian matrix given the current solution guess.
     * Calls a method on a AbstractNonlinearAssemblerSolverHybrid object to do the calculation.
     *
     * Note this method interface was updated in PETSc 3.5, as was the jacobian matrix interface
     * itself, so you'll need to update your Jacobian too.
     *
     * @param snes This is not used by us, but required by PETSc.
     * @param currentGuess The solution guess for the current iteration.
     * @param globalJacobian Pointer to object to fill with the Jacobian matrix.
     * @param preconditioner This is not used by us, but required by PETSc.
     * @param pContext Pointer to a AbstractNonlinearAssemblerSolverHybrid object.
     * @return PETSc expects this function to return a PetscErrorCode. We always return 0;
     *   exceptions are thrown if there is an error.
     *
     * Note: this is a global function, hence the need a long name to avoid
     * potential conflicting names
     */
    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>

    PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian(SNES snes,
                                                                          Vec currentGuess,
                                                                          Mat globalJacobian,
                                                                          Mat preconditioner,
                                                                          void* pContext);
#else
    /**
     * Function called by PETSc to compute the Jacobian matrix given the current solution guess.
     * Calls a method on a AbstractNonlinearAssemblerSolverHybrid object to do the calculation.
     *
     * @param snes This is not used by us, but required by PETSc.
     * @param currentGuess The solution guess for the current iteration.
     * @param pGlobalJacobian Pointer to object to fill with the Jacobian matrix.
     * @param pPreconditioner This is not used by us, but required by PETSc.
     * @param pMatStructure This is not used by us, but required by PETSc.
     * @param pContext Pointer to a AbstractNonlinearAssemblerSolverHybrid object.
     * @return PETSc expects this function to return a PetscErrorCode. We always return 0;
     *   exceptions are thrown if there is an error.
     *
     * Note: this is a global function, hence the need a long name to avoid
     * potential conflicting names
     */
    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
    PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian(SNES snes,
                                                                          Vec currentGuess,
                                                                          Mat* pGlobalJacobian,
                                                                          Mat* pPreconditioner,
                                                                          MatStructure* pMatStructure,
                                                                          void* pContext);
#endif

/**
 * The ASSEMBLER-SOLVER class.
 *
 * An abstract solver for solving nonlinear elliptic PDEs.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractNonlinearAssemblerSolverHybrid : public AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,true,NONLINEAR>
{
protected:

    /** The mesh. */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /** Boundary conditions container. */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* mpBoundaryConditions;

    /** The nonlinear solver. */
    AbstractNonlinearSolver* mpSolver;

    /** Whether memory has been allocated for the solver. */
    bool mWeAllocatedSolverMemory;

    /** Whether to use an analytical expression for the Jacobian. */
    bool mUseAnalyticalJacobian;

    /** An assembler for the surface integral terms in the residual */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* mpNeumannSurfaceTermsAssembler;

    /**
     * Apply Dirichlet boundary conditions to either the residual or Jacobian.
     *
     * @param currentGuess  the solution guess for the current nonlinear solver iteration
     * @param residual the residual to apply boundary conditions to (can be NULL if bcs to
     *   be applied to Jacobian only)
     * @param pMat the Jacobian to apply boundary conditions to (can be NULL if bcs to
     *   be applied to residual only)
     */
    void ApplyDirichletConditions(Vec currentGuess, Vec residual, Mat* pMat);

    /**
     * Computes the Jacobian numerically i.e. an approximation, using numerical
     * partial derivatives.
     *
     * @param currentGuess Independent variable, u in f(u), for example
     * @param pJacobian A pointer to the Jacobian matrix
     */
    void ComputeJacobianNumerically(const Vec currentGuess, Mat* pJacobian);

public:

    /**
     * Compute the residual vector given the current solution guess.
     *
     * @param currentGuess The solution guess for the current nonlinear solve iteration.
     * @param residualVector We fill this with the residual vector.
     *
     * NOTE: this method is called indirectly by the PETSc iterative
     * solvers, so must be public.
     */
    void ComputeResidual(const Vec currentGuess, Vec residualVector);

    /**
     * Compute the Jacobian matrix given a current guess at the solution.
     * Choose whether to use a numerical or analytical method based on a flag
     * provided by the user (in Solve()).
     *
     * @param currentGuess The solution guess for the current iteration.
     * @param pJacobian Pointer to object to fill with the Jacobian matrix.
     *
     * NOTE: this method is called indirectly by the PETSc iterative
     * solvers, so must be public.
     */
    void ComputeJacobian(const Vec currentGuess, Mat* pJacobian);

    /**
     * Constructor.
     *
     * @param pMesh The mesh
     * @param pBoundaryConditions The boundary conditions to apply
     */
    AbstractNonlinearAssemblerSolverHybrid(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                           BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions);

    /**
     * Destructor.
     */
    virtual ~AbstractNonlinearAssemblerSolverHybrid();

    /**
      * Assemble and solve the system for a nonlinear elliptic PDE.
      *
      * @param initialGuess An initial guess for the iterative solver
      * @param useAnalyticalJacobian Set to true to use an analytically calculated
      *     Jacobian matrix rather than a numerically approximated one.
      *
      * @return A PETSc vector giving the solution at each mesh node.
      */
    virtual Vec Solve(Vec initialGuess, bool useAnalyticalJacobian=false);

    /**
      * SetNonlinearSolver - by default a SimplePetscNonlinearSolver is created
      * and used in this class, this method can be called to use a different
      * AbstractNonlinearSolver.
      *
      * @param pNonlinearSolver a nonlinear solver
      */
    void SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver);

    /**
     * A helper method for use when writing concrete assemblers. Once the user has calculated
     * (on paper) the weak form and the form of the ComputeMatrixTerm method, they can check
     * whether the analytic Jacobian matches the numerical Jacobian to verify the
     * correctness of the code.
     *
     * @param tol A tolerance which defaults to 1e-5
     * @return true if the componentwise difference between the matrices is less than
     *   the tolerance, false otherwise.
     *
     * This method should NOT be run in simulations - it is only to verify the correctness
     * of the concrete assembler code. Allocates dense matrices.
     */
    bool VerifyJacobian(double tol);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractNonlinearAssemblerSolverHybrid(
            AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
            BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pBoundaryConditions)
    :  AbstractFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,true,NONLINEAR>(pMesh),
       mpMesh(pMesh),
       mpBoundaryConditions(pBoundaryConditions)
{
    /*
     * Note: if this is run with SPACE_DIM != ELEMENT_DIM the class has to be checked:
     * there may be lots of places where we should be using SPACE_DIM not ELEMENT_DIM.
     */
    assert(SPACE_DIM==ELEMENT_DIM);
    assert(pMesh!=nullptr);
    assert(pBoundaryConditions!=nullptr);

    mpSolver = new SimplePetscNonlinearSolver;
    mWeAllocatedSolverMemory = true;

    assert(mpMesh->GetNumNodes() == mpMesh->GetDistributedVectorFactory()->GetProblemSize());

    mpNeumannSurfaceTermsAssembler = new NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(pMesh,pBoundaryConditions);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~AbstractNonlinearAssemblerSolverHybrid()
{
    if (mWeAllocatedSolverMemory)
    {
        delete mpSolver;
    }
    delete mpNeumannSurfaceTermsAssembler;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ApplyDirichletConditions(Vec currentGuess, Vec residual, Mat* pJacobian)
{
    if (residual)
    {
        mpBoundaryConditions->ApplyDirichletToNonlinearResidual(currentGuess,
                                                                residual,
                                                                *(mpMesh->GetDistributedVectorFactory()));
    }
    if (pJacobian)
    {
        mpBoundaryConditions->ApplyDirichletToNonlinearJacobian(*pJacobian);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeResidual(const Vec currentGuess, Vec residualVector)
{
    // Add the volume integral part of the residual. This will call ComputeVectorTerm, which needs to be
    // implemented in the concrete class.
    this->SetVectorToAssemble(residualVector,true);
    this->SetCurrentSolution(currentGuess);
    this->AssembleVector();

    // Add the surface integral contribution - note the negative sign (scale factor).
    // This contribution is, for example for a 1 unknown problem
    // -integral(g\phi_id dS), where g is the specfied neumann boundary condition
    mpNeumannSurfaceTermsAssembler->SetVectorToAssemble(residualVector, false);
    mpNeumannSurfaceTermsAssembler->SetScaleFactor(-1.0);
    mpNeumannSurfaceTermsAssembler->Assemble();

    PetscVecTools::Finalise(residualVector);

    ApplyDirichletConditions(currentGuess, residualVector, nullptr);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeJacobian(const Vec currentGuess, Mat* pJacobian)
{
    if (mUseAnalyticalJacobian)
    {
        this->SetMatrixToAssemble(*pJacobian);
        this->SetCurrentSolution(currentGuess);
        this->AssembleMatrix();

        PetscMatTools::SwitchWriteMode(*pJacobian);

        ApplyDirichletConditions(currentGuess, nullptr, pJacobian);

        PetscMatTools::Finalise(*pJacobian);
    }
    else
    {
        ComputeJacobianNumerically(currentGuess, pJacobian);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeJacobianNumerically(const Vec currentGuess, Mat* pJacobian)
{
    unsigned num_unknowns = PROBLEM_DIM * this->mpMesh->GetNumNodes();

    // Set up working vectors
    Vec residual;
    Vec perturbed_residual;
    Vec result;
    residual = PetscTools::CreateVec(num_unknowns);
    result = PetscTools::CreateVec(num_unknowns);
    perturbed_residual = PetscTools::CreateVec(num_unknowns);

    // Copy the currentGuess vector; we perturb the copy
    Vec current_guess_copy;
    PETSCEXCEPT( VecDuplicate(currentGuess, &current_guess_copy) );
    PETSCEXCEPT( VecCopy(currentGuess, current_guess_copy) );

    // Compute the current residual
    ComputeResidual(currentGuess, residual);

    // Amount to perturb each input element by
    double h = 1e-5;
    double near_hsquared = 1e-9;

    PetscInt ilo, ihi;
    VecGetOwnershipRange(current_guess_copy, &ilo, &ihi);
    unsigned lo = ilo;
    unsigned hi = ihi;

    // Iterate over entries in the input vector
    for (unsigned global_index_outer = 0; global_index_outer < num_unknowns; global_index_outer++)
    {
        // Only perturb if we own it
        PetscVecTools::AddToElement(current_guess_copy, global_index_outer, h);
        ComputeResidual(current_guess_copy, perturbed_residual);

        // result = (perturbed_residual - residual) / h
        PetscVecTools::WAXPY(result, -1.0, residual, perturbed_residual);
        PetscVecTools::Scale(result, 1.0/h);

        double* p_result;
        ///\todo This loop is setting the column "global_index_outer" of
        /// the Jacobian matrix to the result vector in a non-sparse way.
        PETSCEXCEPT( VecGetArray(result, &p_result) );
        for (unsigned global_index=lo; global_index < hi; global_index++)
        {
            double result_entry = p_result[ global_index - lo];
            if (!CompareDoubles::IsNearZero(result_entry, near_hsquared))
            {
                PetscMatTools::SetElement(*pJacobian, global_index, global_index_outer, result_entry);
            }
        }
        PETSCEXCEPT( VecRestoreArray(result, &p_result) );

        PetscVecTools::AddToElement(current_guess_copy, global_index_outer, -h);
    }

    PetscTools::Destroy(residual);
    PetscTools::Destroy(perturbed_residual);
    PetscTools::Destroy(result);
    PetscTools::Destroy(current_guess_copy);

    PetscMatTools::Finalise(*pJacobian);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve(Vec initialGuess,
                                                                                       bool useAnalyticalJacobian)
{
    assert(initialGuess != nullptr);
    mUseAnalyticalJacobian = useAnalyticalJacobian;

    PetscInt size_of_init_guess;
    VecGetSize(initialGuess, &size_of_init_guess);
    PetscInt problem_size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
    if (size_of_init_guess != problem_size)
    {
        EXCEPTION("Size of initial guess vector, " << size_of_init_guess
                      << ", does not match size of problem, " << problem_size);
    }

    // Run the solver, telling it which global functions to call in order to assemble the residual or jacobian
    return mpSolver->Solve(&AbstractNonlinearAssemblerSolverHybrid_ComputeResidual<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>,
                           &AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>,
                           initialGuess,
                           PROBLEM_DIM * this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                           this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetNonlinearSolver(AbstractNonlinearSolver* pNonlinearSolver)
{
    if (mWeAllocatedSolverMemory)
    {
        delete mpSolver;
    }
    mpSolver = pNonlinearSolver;
    mWeAllocatedSolverMemory = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::VerifyJacobian(double tol)
{
    unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();

    Vec initial_guess = PetscTools::CreateAndSetVec(size, 0.0);

    Mat analytic_jacobian;
    Mat numerical_jacobian;

    PetscTools::SetupMat(analytic_jacobian, size, size, size);
    PetscTools::SetupMat(numerical_jacobian, size, size, size);

    mUseAnalyticalJacobian = true;
    ComputeJacobian(initial_guess, &analytic_jacobian);

    mUseAnalyticalJacobian = false;
    ComputeJacobian(initial_guess, &numerical_jacobian);

    bool ok = PetscMatTools::CheckEquality(numerical_jacobian, analytic_jacobian, tol);
    PetscTools::Destroy(numerical_jacobian);
    PetscTools::Destroy(analytic_jacobian);
    PetscTools::Destroy(initial_guess);

    return ok;
}

// Implementations of GLOBAL functions called by PETSc nonlinear solver

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeResidual(SNES snes,
                                                                      Vec currentGuess,
                                                                      Vec residualVector,
                                                                      void* pContext)
{
    // Extract the solver from the void*
    AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_solver =
        (AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;

    p_solver->ComputeResidual(currentGuess, residualVector);

    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian(SNES snes,
                                                                      Vec currentGuess,
                                                                      Mat globalJacobian,
                                                                      Mat preconditioner,
                                                                      void* pContext)
{
    // Extract the solver from the void*
    AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_solver =
        (AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;

    p_solver->ComputeJacobian(currentGuess, &globalJacobian);

    return 0;
}
#else
PetscErrorCode AbstractNonlinearAssemblerSolverHybrid_ComputeJacobian(SNES snes,
                                                                      Vec currentGuess,
                                                                      Mat* pGlobalJacobian,
                                                                      Mat* pPreconditioner,
                                                                      MatStructure* pMatStructure,
                                                                      void* pContext)
{
    // Extract the solver from the void*
    AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* p_solver =
        (AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>*) pContext;

    p_solver->ComputeJacobian(currentGuess, pGlobalJacobian);

    return 0;
}
#endif

#endif /*ABSTRACTNONLINEARASSEMBLERSOLVERHYBRID_HPP_*/
