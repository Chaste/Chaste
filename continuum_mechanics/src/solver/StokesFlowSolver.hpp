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
#include "StokesFlowPreconditionerAssembler.hpp"
#include "ContinuumMechanicsNeumannBcsAssembler.hpp"

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

    /** Assembler for computing volume integral part of preconditioner matrix (which
     *  is the same as the system matrix except has a mass matrix in the pressure-pressure
     *  block
     */
    StokesFlowPreconditionerAssembler<DIM>* mpStokesFlowPreconditionerAssembler;

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

public:

    /**
     * Constructor.
     *
     * @param rQuadMesh Quadratic mesh
     * @param rProblemDefinition Problem definition
     * @param outputDirectory the output directory to use
     */
    StokesFlowSolver(AbstractTetrahedralMesh<DIM,DIM>& rQuadMesh,
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
     * @return the flow.
     * Note: return_value[i](j) = u_j for node i.
     */
    std::vector<c_vector<double,DIM> >& rGetSpatialSolution();

    /**
     * @return the flow. Note: return_value[i](j) = u_j for node i. Just
     * calls rGetSpatialSolution().
     */
    std::vector<c_vector<double,DIM> >& rGetVelocities();
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
StokesFlowSolver<DIM>::StokesFlowSolver(AbstractTetrahedralMesh<DIM,DIM>& rQuadMesh,
                                        StokesFlowProblemDefinition<DIM>& rProblemDefinition,
                                        std::string outputDirectory)
    : AbstractContinuumMechanicsSolver<DIM>(rQuadMesh, rProblemDefinition, outputDirectory, INCOMPRESSIBLE),
      mrProblemDefinition(rProblemDefinition),
      mKspAbsoluteTol(-1)
{
    assert(DIM==2 || DIM==3);
    assert(!mrProblemDefinition.rGetDirichletNodes().empty());

    mpStokesFlowAssembler = new StokesFlowAssembler<DIM>(&this->mrQuadMesh, &mrProblemDefinition);
    mpStokesFlowPreconditionerAssembler = new StokesFlowPreconditionerAssembler<DIM>(&this->mrQuadMesh, &mrProblemDefinition);
    mpNeumannBcsAssembler = new ContinuumMechanicsNeumannBcsAssembler<DIM>(&this->mrQuadMesh, &mrProblemDefinition);
}

template<unsigned DIM>
StokesFlowSolver<DIM>::~StokesFlowSolver()
{
    delete mpStokesFlowAssembler;
    delete mpStokesFlowPreconditionerAssembler;
    delete mpNeumannBcsAssembler;
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::Solve()
{
    mrProblemDefinition.Validate();

    if (this->mVerbose)
    {
        Timer::Reset();
    }

    // Assemble Jacobian (and preconditioner)
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ASSEMBLE);
    AssembleSystem();
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ASSEMBLE);

    if (this->mVerbose)
    {
        Timer::PrintAndReset("AssembleSystem");
    }

    /*
     * Solve the linear system using PETSc GMRES and an LU factorisation
     * of the preconditioner. Note we don't call Solve on the linear_system
     * as we want to set PETSc options.
     */
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::SOLVE);

    Vec solution;
    VecDuplicate(this->mLinearSystemRhsVector,&solution);
    PetscVecTools::Zero(solution);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD,&solver);
#if ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR>=5))
    KSPSetOperators(solver, this->mSystemLhsMatrix, this->mPreconditionMatrix);
#else
    KSPSetOperators(solver, this->mSystemLhsMatrix, this->mPreconditionMatrix, DIFFERENT_NONZERO_PATTERN /*in precond between successive solves*/);
#endif
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

    PC pc;
    KSPGetPC(solver, &pc);
    PCSetType(pc, PCJACOBI);

    KSPSetUp(solver);

    KSPSetFromOptions(solver);

//    ///// For printing matrix when debugging
//    OutputFileHandler handler("TEMP");
//    out_stream p_file = handler.OpenOutputFile("matrix.txt");
//    for (unsigned i=0; i<this->mNumDofs; i++)
//    {
//        for (unsigned j=0; j<this->mNumDofs; j++)
//        {
//            *p_file << PetscMatTools::GetElement(this->mSystemLhsMatrix, i, j) << " ";
//        }
//        *p_file << "\n";
//    }
//    p_file->close();
//
//    out_stream p_file2 = handler.OpenOutputFile("rhs.txt");
//    for (unsigned i=0; i<this->mNumDofs; i++)
//    {
//        *p_file2 << PetscVecTools::GetElement(this->mLinearSystemRhsVector, i) << "\n";
//    }
//    p_file2->close();

    if (this->mVerbose)
    {
        Timer::PrintAndReset("KSP Setup");
    }

    KSPSolve(solver,this->mLinearSystemRhsVector,solution);

    KSPConvergedReason reason;
    KSPGetConvergedReason(solver,&reason);
    KSPEXCEPT(reason);

    if (this->mVerbose)
    {
        Timer::PrintAndReset("KSP Solve");
        int num_iters;
        KSPGetIterationNumber(solver, &num_iters);
        std::cout << "[" << PetscTools::GetMyRank() << "]: Num iterations = " << num_iters << "\n" << std::flush;
    }

    MechanicsEventHandler::EndEvent(MechanicsEventHandler::SOLVE);

    // Copy solution into the std::vector
    ReplicatableVector solution_repl(solution);
    for (unsigned i=0; i<this->mNumDofs; i++)
    {
        this->mCurrentSolution[i] = solution_repl[i];
    }

    // Remove pressure dummy values (P=0 at internal nodes, which should have been
    // been the result of the solve above), by linear interpolating from vertices of
    // edges to the internal node
    this->RemovePressureDummyValuesThroughLinearInterpolation();

    PetscTools::Destroy(solution);
    KSPDestroy(PETSC_DESTROY_PARAM(solver));

    this->WriteCurrentSpatialSolution("flow_solution", "nodes");
    this->WriteCurrentPressureSolution();
}

template<unsigned DIM>
void StokesFlowSolver<DIM>::AssembleSystem()
{
    // Use assembler to assemble volume integral part....
    mpStokesFlowAssembler->SetMatrixToAssemble(this->mSystemLhsMatrix, true);
    mpStokesFlowAssembler->SetVectorToAssemble(this->mLinearSystemRhsVector, true);
    mpStokesFlowAssembler->Assemble();

    mpStokesFlowPreconditionerAssembler->SetMatrixToAssemble(this->mPreconditionMatrix, true);
    mpStokesFlowPreconditionerAssembler->AssembleMatrix();

    mpNeumannBcsAssembler->SetVectorToAssemble(this->mLinearSystemRhsVector, false /*don't zero!*/);
    mpNeumannBcsAssembler->AssembleVector();

    PetscVecTools::Finalise(this->mLinearSystemRhsVector);
    PetscMatTools::SwitchWriteMode(this->mSystemLhsMatrix);
    PetscMatTools::SwitchWriteMode(this->mPreconditionMatrix);

    // Note: maintaining symmetry for Dirichlet BCs is possible at the moment (the second parameter)
    // but doing so for the identity block is not yet implemented (the second parameter must be false)
    // Not sure if maintaining symmetry is worth it - may allow CG to work, but matrix is indefinite
    // so GC not guaranteed to work..
    //
    // Note: the identity block needs to be added before the BCs - see comments in
    // PetscMatTools::ZeroRowsWithValueOnDiagonal()
    this->AddIdentityBlockForDummyPressureVariables(LINEAR_PROBLEM);
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
            // DIM+1 is the problem dimension
            this->mSpatialSolution[i](j) = this->mCurrentSolution[(DIM+1)*i+j];
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
