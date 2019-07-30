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

#include "MonodomainSolver.hpp"
#include "MassMatrixAssembler.hpp"
#include "PetscMatTools.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainSolver<ELEMENT_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);


    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix)
    /////////////////////////////////////////
    if (computeMatrix)
    {
        mpMonodomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        mpMonodomainAssembler->AssembleMatrix();

        MassMatrixAssembler<ELEMENT_DIM,SPACE_DIM> mass_matrix_assembler(this->mpMesh, HeartConfig::Instance()->GetUseMassLumping());
        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.Assemble();

        this->mpLinearSystem->FinaliseLhsMatrix();
        PetscMatTools::Finalise(mMassMatrix);

        if (HeartConfig::Instance()->GetUseMassLumpingForPrecond() && !HeartConfig::Instance()->GetUseMassLumping())
        {
            this->mpLinearSystem->SetPrecondMatrixIsDifferentFromLhs();

            MonodomainAssembler<ELEMENT_DIM,SPACE_DIM> lumped_mass_assembler(this->mpMesh,this->mpMonodomainTissue);
            lumped_mass_assembler.SetMatrixToAssemble(this->mpLinearSystem->rGetPrecondMatrix());

            HeartConfig::Instance()->SetUseMassLumping(true);
            lumped_mass_assembler.AssembleMatrix();
            HeartConfig::Instance()->SetUseMassLumping(false);

            this->mpLinearSystem->FinalisePrecondMatrix();
        }
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

    //////////////////////////////////////////
    // Set up z in b=Mz
    //////////////////////////////////////////
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    // dist stripe for z (return value)
    DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(mVecForConstructingRhs);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm = HeartConfig::Instance()->GetCapacitance();

    for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
         index!= dist_vec_matrix_based.End();
         ++index)
    {
        double V = distributed_current_solution[index];
        double F = - Am*this->mpMonodomainTissue->rGetIionicCacheReplicated()[index.Global]
                   - this->mpMonodomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];

        dist_vec_matrix_based[index] = Am*Cm*V*PdeSimulationTime::GetPdeTimeStepInverse() + F;
    }
    dist_vec_matrix_based.Restore();

    //////////////////////////////////////////
    // b = Mz
    //////////////////////////////////////////
    MatMult(mMassMatrix, mVecForConstructingRhs, this->mpLinearSystem->rGetRhsVector());

    // assembling RHS is not finished yet, as Neumann bcs are added below, but
    // the event will be begun again inside mpMonodomainAssembler->AssembleVector();
    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);

    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////
    mpNeumannSurfaceTermsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
    mpNeumannSurfaceTermsAssembler->AssembleVector();

    /////////////////////////////////////////
    // apply correction term
    /////////////////////////////////////////
    if (mpMonodomainCorrectionTermAssembler)
    {
        mpMonodomainCorrectionTermAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
        // don't need to set current solution
        mpMonodomainCorrectionTermAssembler->AssembleVector();
    }

    // finalise
    this->mpLinearSystem->FinaliseRhsVector();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }

    // call base class version...
    AbstractLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>::InitialiseForSolve(initialSolution);

    //..then do a bit extra
    if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
        this->mpLinearSystem->SetAbsoluteTolerance(HeartConfig::Instance()->GetAbsoluteTolerance());
    }
    else
    {
        this->mpLinearSystem->SetRelativeTolerance(HeartConfig::Instance()->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetUseFixedNumberIterations(HeartConfig::Instance()->GetUseFixedNumberIterationsLinearSolver(), HeartConfig::Instance()->GetEvaluateNumItsEveryNSolves());

    // initialise matrix-based RHS vector and matrix, and use the linear
    // system rhs as a template
    Vec& r_template = this->mpLinearSystem->rGetRhsVector();
    VecDuplicate(r_template, &mVecForConstructingRhs);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(r_template, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    PetscTools::SetupMat(mMassMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(),
                         this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
    // solve cell models
    mpMonodomainTissue->SolveCellSystems(currentSolution, PdeSimulationTime::GetTime(), PdeSimulationTime::GetNextTime());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainSolver<ELEMENT_DIM,SPACE_DIM>::MonodomainSolver(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
            BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>(pMesh),
      mpMonodomainTissue(pTissue),
      mpBoundaryConditions(pBoundaryConditions)
{
    assert(pTissue);
    assert(pBoundaryConditions);
    this->mMatrixIsConstant = true;

    mpMonodomainAssembler = new MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpMonodomainTissue);
    mpNeumannSurfaceTermsAssembler = new NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,1>(pMesh,pBoundaryConditions);


    // Tell tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs = NULL;

    if (HeartConfig::Instance()->GetUseStateVariableInterpolation())
    {
        mpMonodomainCorrectionTermAssembler
            = new MonodomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpMonodomainTissue);
        //We are going to need those caches after all
        pTissue->SetCacheReplication(true);
    }
    else
    {
        mpMonodomainCorrectionTermAssembler = NULL;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainSolver<ELEMENT_DIM,SPACE_DIM>::~MonodomainSolver()
{
    delete mpMonodomainAssembler;
    delete mpNeumannSurfaceTermsAssembler;

    if (mVecForConstructingRhs)
    {
        PetscTools::Destroy(mVecForConstructingRhs);
        PetscTools::Destroy(mMassMatrix);
    }

    if (mpMonodomainCorrectionTermAssembler)
    {
        delete mpMonodomainCorrectionTermAssembler;
    }
}

// Explicit instantiation
template class MonodomainSolver<1,1>;
template class MonodomainSolver<1,2>;
template class MonodomainSolver<1,3>;
template class MonodomainSolver<2,2>;
template class MonodomainSolver<3,3>;
