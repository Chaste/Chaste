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

#ifndef STOKESFLOWSOLVER_HPP_
#define STOKESFLOWSOLVER_HPP_

#include "AbstractContinuumMechanicsSolver.hpp"
#include "LinearSystem.hpp"
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include "Timer.hpp"
#include "PetscMatTools.hpp"
#include "PetscVecTools.hpp"
#include "StokesFlowProblemDefinition.hpp"
#include "StokesFlowAssembler.hpp"
#include "ContinuumMechanicsNeumannBcsAssembler.hpp"

#define STOKES_VERBOSE

/**
 * Finite element solver for Stokes flow problems
 */
template<unsigned DIM>
class StokesFlowSolver : public AbstractContinuumMechanicsSolver<DIM>
{
    friend class TestStokesFlowSolver;

private:
    /** Object containing all the information about the problem to solve */
    StokesFlowProblemDefinition<DIM>& mrProblemDefinition;

    /** Assembler for computing volume integral part of matrix and RHS vector */
    StokesFlowAssembler<DIM>* mpStokesFlowAssembler;

    /**
     *  Assembler for adding the surface integral arising from natural
     *  Neumman boundary conditions to the RHS vector
     */
    ContinuumMechanicsNeumannBcsAssembler<DIM>* mpNeumannBcsAssembler;

    /**
     * Absolute tolerance for linear systems. Can be set by calling
     * SetKspAbsoluteTolerances(), but default to -1, in which case
     * a relative tolerance is used.
     */
    double mKspAbsoluteTol;

    /**
     * Assemble the linear system and preconditioner matrix.
     */
    void AssembleSystem();

    /**
     * Apply the Dirichlet boundary conditions to the linear system and preconditioner matrix.
     */
    void ApplyBoundaryConditions();

public:

    /**
     * Constructor.
     *
     * @param rQuadMesh Quadratic mesh
     * @param rProblemDefinition Problem definition
     * @param outputDirectory the output directory to use
     */
    StokesFlowSolver(QuadraticMesh<DIM>& rQuadMesh,
                     StokesFlowProblemDefinition<DIM>& rProblemDefinition,
                     std::string outputDirectory);

    /**
     * Destructor.
     */
    virtual ~StokesFlowSolver();

    /**
     * Solve the system.
     */
    void Solve();

    /**
     * Set the absolute tolerance to be used when solving the linear system.
     * If this is not called a relative tolerance is used.
     *
     * @param kspAbsoluteTolerance the tolerance
     */
    void SetKspAbsoluteTolerance(double kspAbsoluteTolerance);


    /**
     * Implemented method, returns the flow.
     * Note: return_value[i](j) = u_j for node i.
     */
    std::vector<c_vector<double,DIM> >& rGetSpatialSolution();

    /**
     * Get the flow. Note: return_value[i](j) = u_j for node i. Just
     * calls rGetSpatialSolution().
     */
    std::vector<c_vector<double,DIM> >& rGetVelocities();
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void StokesFlowSolver<DIM>::ApplyBoundaryConditions()
{
    std::vector<unsigned> rows;
    rows.resize(DIM*mrProblemDefinition.rGetDirichletNodes().size());

    for (unsigned i=0; i<mrProblemDefinition.rGetDirichletNodes().size(); i++)
    {
        unsigned node_index = mrProblemDefinition.rGetDirichletNodes()[i];
        for (unsigned j=0; j<DIM; j++)
        {
            unsigned dof_index = DIM*node_index + j;
            rows[DIM*i + j] = dof_index;

            double value = mrProblemDefinition.rGetDirichletNodeValues()[i](j);
            PetscVecTools::SetElement(this->mLinearSystemRhsVector, dof_index, value);
        }
    }

    PetscMatTools::ZeroRowsWithValueOnDiagonal(this->mSystemLhsMatrix, rows, 1.0);
    PetscMatTools::ZeroRowsWithValueOnDiagonal(this->mPreconditionMatrix, rows, 1.0);
}

template<unsigned DIM>
StokesFlowSolver<DIM>::StokesFlowSolver(QuadraticMesh<DIM>& rQuadMesh,
                                        StokesFlowProblemDefinition<DIM>& rProblemDefinition,
                                        std::string outputDirectory)
    : AbstractContinuumMechanicsSolver<DIM>(rQuadMesh, rProblemDefinition, outputDirectory, INCOMPRESSIBLE),
      mrProblemDefinition(rProblemDefinition),
      mKspAbsoluteTol(-1)
{
    assert(DIM==2 || DIM==3);
    assert(!mrProblemDefinition.rGetDirichletNodes().empty());

    mpStokesFlowAssembler = new StokesFlowAssembler<DIM>(&this->mrQuadMesh, &mrProblemDefinition);
    mpNeumannBcsAssembler = new ContinuumMechanicsNeumannBcsAssembler<DIM>(&this->mrQuadMesh, &mrProblemDefinition);
}

template<unsigned DIM>
StokesFlowSolver<DIM>::~StokesFlowSolver()
{
    delete mpStokesFlowAssembler;
    delete mpNeumannBcsAssembler;
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::Solve()
{
    mrProblemDefinition.Validate();

    #ifdef STOKES_VERBOSE
    Timer::Reset();
    #endif

    // Assemble Jacobian (and preconditioner)
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ASSEMBLE);
    AssembleSystem();
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ASSEMBLE);
    #ifdef STOKES_VERBOSE
    Timer::PrintAndReset("AssembleSystem");
    #endif

    /*
     * Solve the linear system using PETSc GMRES and an LU factorisation
     * of the preconditioner. Note we don't call Solve on the linear_system
     * as we want to set PETSc options.
     */
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::SOLVE);

    Vec solution;
    VecDuplicate(this->mLinearSystemRhsVector,&solution);


    KSP solver;
    KSPCreate(PETSC_COMM_WORLD,&solver);

    KSPSetOperators(solver, this->mSystemLhsMatrix, this->mPreconditionMatrix, DIFFERENT_NONZERO_PATTERN /*in precond between successive solves*/);

    KSPSetType(solver, KSPGMRES);

    if (mKspAbsoluteTol < 0)
    {
        double ksp_rel_tol = 1e-6;
        KSPSetTolerances(solver, ksp_rel_tol, PETSC_DEFAULT, PETSC_DEFAULT, 10000 /*max iter*/); //hopefully with the preconditioner this max is way too high
    }
    else
    {
        KSPSetTolerances(solver, 1e-16, mKspAbsoluteTol, PETSC_DEFAULT, 10000 /*max iter*/); //hopefully with the preconditioner this max is way too high
    }

    unsigned num_restarts = 100;
    KSPGMRESSetRestart(solver,num_restarts); // gmres num restarts

    KSPSetFromOptions(solver);
    KSPSetUp(solver);
    #ifdef STOKES_VERBOSE
    Timer::PrintAndReset("KSP Setup");
    #endif

    PC pc;
    KSPGetPC(solver, &pc);

/////// What was going on before, when hypre was being used...
//    #ifndef *****
    PCSetType(pc, PCBJACOBI); // BJACOBI = ILU on each block (block = part of matrix on each process)
//    #else
//    /////////////////////////////////////////////////////////////////////////////////////////////////////
//    // Speed up linear solve time massively for larger simulations (in fact GMRES may stagnate without
//    // this for larger problems), by using a AMG preconditioner -- needs HYPRE installed
//    /////////////////////////////////////////////////////////////////////////////////////////////////////
//    PetscOptionsSetValue("-pc_hypre_type", "boomeramg");
//    // PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
//    // PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");
//
//    PCSetType(pc, PCHYPRE);
//
//    //PCBlockDiagonalMechanics* p_custom_pc = new PCBlockDiagonalMechanics(solver, r_precond_jac, mBlock1Size, mBlock2Size);
//    //PCLDUFactorisationMechanics* p_custom_pc = new PCLDUFactorisationMechanics(solver, r_precond_jac, mBlock1Size, mBlock2Size);
//    //remember to delete memory..
//    //KSPSetPreconditionerSide(solver, PC_RIGHT);
//    #endif

    KSPSetFromOptions(solver);

    KSPSolve(solver,this->mLinearSystemRhsVector,solution);

//    std::cout << "RHS\n";
//    PetscVecTools::Display(this->mLinearSystemRhsVector);
//
//    std::cout << "Matrix\n";
//    for (unsigned i=0; i<22; i++)
//    {
//        for (unsigned j=0; j<22; j++)
//        {
//            double val = PetscMatTools::GetElement(this->mSystemLhsMatrix, i, j);
//            if (fabs(val)<1e-9)
//            {
//                val = 0.0;
//            }
//            std::cout << val << " ";
//        }
//        std::cout << "\n";
//    }
//
//    PetscMatTools::Display(r_jac);
//    std::cout << "Solution\n";
//     PetscVecTools::Display(solution);

    KSPConvergedReason reason;
    KSPGetConvergedReason(solver,&reason);
    KSPEXCEPT(reason);
    std::cout << "Converged Reason = " << reason << "\n" << std::flush;

    #ifdef STOKES_VERBOSE
    Timer::PrintAndReset("KSP Solve");
    int num_iters;
    KSPGetIterationNumber(solver, &num_iters);
    std::cout << "[" << PetscTools::GetMyRank() << "]: Num iterations = " << num_iters << "\n" << std::flush;
    #endif

    MechanicsEventHandler::EndEvent(MechanicsEventHandler::SOLVE);

///\todo: three copies?!
    // Copy solution into the std::vector
    ReplicatableVector solution_repl(solution);
    for (unsigned i=0; i<this->mNumDofs; i++)
    {
        this->mCurrentSolution[i] = solution_repl[i];
    }

    PetscTools::Destroy(solution);
    KSPDestroy(PETSC_DESTROY_PARAM(solver));

    this->WriteCurrentSpatialSolution("flow_solution", "nodes");
    this->WriteCurrentPressureSolution();
}



template<unsigned DIM>
void StokesFlowSolver<DIM>::AssembleSystem()
{
    PetscVecTools::Zero(this->mLinearSystemRhsVector);
    PetscMatTools::Zero(this->mSystemLhsMatrix);
    PetscMatTools::Zero(this->mPreconditionMatrix);


    // Use assembler to assemble volume integral part....
    mpStokesFlowAssembler->SetMatrixToAssemble(this->mSystemLhsMatrix, true);
    mpStokesFlowAssembler->SetVectorToAssemble(this->mLinearSystemRhsVector, true);
    mpStokesFlowAssembler->Assemble();

///\todo! don't use the same assembler for this, use one that puts in C=M...
    mpStokesFlowAssembler->SetMatrixToAssemble(this->mPreconditionMatrix, true);
    mpStokesFlowAssembler->AssembleMatrix();

    mpNeumannBcsAssembler->SetVectorToAssemble(this->mLinearSystemRhsVector, false /*don't zero!*/);
    mpNeumannBcsAssembler->AssembleVector();

    PetscVecTools::Finalise(this->mLinearSystemRhsVector);
    PetscMatTools::SwitchWriteMode(this->mSystemLhsMatrix);
    PetscMatTools::SwitchWriteMode(this->mPreconditionMatrix);

///\todo! Do we really want a symmetric matrix (the true below) - can't use CG after all..
    // Apply Dirichlet boundary conditions
    this->ApplyDirichletBoundaryConditions(LINEAR_PROBLEM, true);

    PetscVecTools::Finalise(this->mLinearSystemRhsVector);
    PetscMatTools::Finalise(this->mSystemLhsMatrix);
    PetscMatTools::Finalise(this->mPreconditionMatrix);
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::SetKspAbsoluteTolerance(double kspAbsoluteTolerance)
{
    assert(kspAbsoluteTolerance > 0);
    mKspAbsoluteTol = kspAbsoluteTolerance;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& StokesFlowSolver<DIM>::rGetSpatialSolution()
{
    this->mSpatialSolution.resize(this->mrQuadMesh.GetNumNodes(), zero_vector<double>(DIM));
    for (unsigned i=0; i<this->mrQuadMesh.GetNumNodes(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            this->mSpatialSolution[i](j) = this->mCurrentSolution[DIM*i+j];
        }
    }
    return this->mSpatialSolution;
}

template<unsigned DIM>
std::vector<c_vector<double,DIM> >& StokesFlowSolver<DIM>::rGetVelocities()
{
    return rGetSpatialSolution();
}



#endif /* STOKESFLOWSOLVER_HPP_ */
