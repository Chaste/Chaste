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

#include "MonodomainPurkinjeSolver.hpp"
#include "MonodomainPurkinjeVolumeMassMatrixAssembler.hpp"
#include "MonodomainPurkinjeCableMassMatrixAssembler.hpp"
#include "PetscMatTools.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::SetupLinearSystem(Vec currentSolution, bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);

    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix)
    /////////////////////////////////////////
    if(computeMatrix)
    {
        mpVolumeAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix(),false);
        mpVolumeAssembler->AssembleMatrix();

        mpCableAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix(),false);
        // False here is to say to the cable assembler don't zero the matrix before assembling,
        // just add the terms to the previous value of the matrix.
        mpCableAssembler->AssembleMatrix();

        SetIdentityBlockToLhsMatrix();
        this->mpLinearSystem->FinaliseLhsMatrix();


        MonodomainPurkinjeVolumeMassMatrixAssembler<ELEMENT_DIM,SPACE_DIM> volume_mass_matrix_assembler(mpMixedMesh, HeartConfig::Instance()->GetUseMassLumping());
        volume_mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        volume_mass_matrix_assembler.Assemble();

        MonodomainPurkinjeCableMassMatrixAssembler<ELEMENT_DIM,SPACE_DIM> cable_mass_matrix_assembler(mpMixedMesh, HeartConfig::Instance()->GetUseMassLumping());
        cable_mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix,false /* don't zero the matrix*/);
        cable_mass_matrix_assembler.Assemble();

        PetscMatTools::Finalise(mMassMatrix);
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

    //////////////////////////////////////////
    // Set up z in b=Mz
    //////////////////////////////////////////
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();

    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    DistributedVector::Stripe distributed_current_solution_volume(distributed_current_solution, 0);
    DistributedVector::Stripe distributed_current_solution_cable(distributed_current_solution, 1);
    // dist stripe for z
    DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(mVecForConstructingRhs);
    DistributedVector::Stripe dist_vec_matrix_based_volume(dist_vec_matrix_based, 0);
    DistributedVector::Stripe dist_vec_matrix_based_cable(dist_vec_matrix_based, 1);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm  = HeartConfig::Instance()->GetCapacitance();

    double Am_purkinje = HeartConfig::Instance()->GetPurkinjeSurfaceAreaToVolumeRatio();
    double Cm_purkinje  = HeartConfig::Instance()->GetPurkinjeCapacitance();


    for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
         index!= dist_vec_matrix_based.End();
         ++index)
    {
        double V_volume = distributed_current_solution_volume[index];
        double F_volume = - Am*this->mpMonodomainTissue->rGetIionicCacheReplicated()[index.Global]
                          - this->mpMonodomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];
        dist_vec_matrix_based_volume[index] = Am*Cm*V_volume*PdeSimulationTime::GetPdeTimeStepInverse() + F_volume;

        double V_cable = distributed_current_solution_cable[index];
        double F_cable = - Am*this->mpMonodomainTissue->rGetPurkinjeIionicCacheReplicated()[index.Global] //Purkinje intra-cell stimulus not defined yet
                         - this->mpMonodomainTissue->rGetPurkinjeIntracellularStimulusCacheReplicated()[index.Global];

        dist_vec_matrix_based_cable[index] = Am_purkinje*Cm_purkinje*V_cable*PdeSimulationTime::GetPdeTimeStepInverse() + F_cable;
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

   // finalise
   this->mpLinearSystem->FinaliseRhsVector();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::SetIdentityBlockToLhsMatrix()
{
    this->mpLinearSystem->FinaliseLhsMatrix();

    Vec diagonal;
    VecDuplicate(mVecForConstructingRhs, &diagonal);
    MatGetDiagonal(this->mpLinearSystem->rGetLhsMatrix(), diagonal);

    // if A(i,i)=0, i must be within the block to be altered, so set A(i,i)=1.
    PetscMatTools::SwitchWriteMode(this->mpLinearSystem->rGetLhsMatrix());
    PetscInt lo, hi;
    this->mpLinearSystem->GetOwnershipRange(lo, hi);
    for (int row=lo; row<hi; row++)
    {
        if ( fabs( PetscVecTools::GetElement(diagonal, row)) < 1e-8 )
        {
            PetscMatTools::SetElement(this->mpLinearSystem->rGetLhsMatrix(),row, row, 1.0);
        }
    }

    PetscTools::Destroy(diagonal);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }

    // call base class version...
    AbstractLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>::InitialiseForSolve(initialSolution);

    //..then do a bit extra
    if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
        this->mpLinearSystem->SetAbsoluteTolerance(HeartConfig::Instance()->GetAbsoluteTolerance());
    }
    else
    {
        NEVER_REACHED;
//        this->mpLinearSystem->SetRelativeTolerance(HeartConfig::Instance()->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetUseFixedNumberIterations(HeartConfig::Instance()->GetUseFixedNumberIterationsLinearSolver(), HeartConfig::Instance()->GetEvaluateNumItsEveryNSolves());

    // Initialise sizes/partitioning of mass matrix & vector, using the initial condition as a template
    VecDuplicate(initialSolution, &mVecForConstructingRhs);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(initialSolution, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    PetscTools::SetupMat(mMassMatrix, 2*this->mpMesh->GetNumNodes(), 2*this->mpMesh->GetNumNodes(),
                         2*this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
    // solve cell models
    mpMonodomainTissue->SolveCellSystems(currentSolution, PdeSimulationTime::GetTime(), PdeSimulationTime::GetNextTime());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::MonodomainPurkinjeSolver(
            MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
            BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>(pMesh),
      mpMixedMesh(pMesh),
      mpMonodomainTissue(pTissue),
      mpBoundaryConditions(pBoundaryConditions)
{
    assert(pTissue);
    assert(pBoundaryConditions);
    if(HeartConfig::Instance()->GetUseStateVariableInterpolation())
    {
        EXCEPTION("State-variable interpolation is not yet supported with Purkinje");
    }
    this->mMatrixIsConstant = true;

    mpVolumeAssembler = new MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>(mpMixedMesh,this->mpMonodomainTissue);
    mpCableAssembler = new MonodomainPurkinjeCableAssembler<ELEMENT_DIM,SPACE_DIM>(mpMixedMesh);
    mpNeumannSurfaceTermsAssembler = new NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,2>(pMesh,pBoundaryConditions);

    // Tell tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs = NULL;

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>::~MonodomainPurkinjeSolver()
{
    delete mpVolumeAssembler;
    delete mpCableAssembler;
    delete mpNeumannSurfaceTermsAssembler;

    if (mVecForConstructingRhs)
    {
        PetscTools::Destroy(mVecForConstructingRhs);
        PetscTools::Destroy(mMassMatrix);
    }
}


///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class MonodomainPurkinjeSolver<2,2>;
template class MonodomainPurkinjeSolver<3,3>;
