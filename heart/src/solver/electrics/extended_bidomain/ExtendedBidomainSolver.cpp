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


#include "ExtendedBidomainSolver.hpp"
#include "ExtendedBidomainMassMatrixAssembler.hpp"
#include "ExtendedBidomainAssembler.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }
    AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(initialSolution);

    // initialise matrix-based RHS vector and matrix, and use the linear
    // system rhs as a template
    Vec& r_template = this->mpLinearSystem->rGetRhsVector();
    VecDuplicate(r_template, &mVecForConstructingRhs);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(r_template, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    PetscTools::SetupMat(mMassMatrix, 3*this->mpMesh->GetNumNodes(), 3*this->mpMesh->GetNumNodes(),
                         3*this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::SetupLinearSystem(
        Vec currentSolution,
        bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);
    assert(currentSolution != NULL);

    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix)
    /////////////////////////////////////////
    if (computeMatrix)
    {
        mpExtendedBidomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        mpExtendedBidomainAssembler->AssembleMatrix();

        // the ExtendedBidomainMassMatrixAssembler deals with the mass matrix
        // for both bath and nonbath problems
        assert(SPACE_DIM==ELEMENT_DIM);
        ExtendedBidomainMassMatrixAssembler<SPACE_DIM> mass_matrix_assembler(this->mpMesh);
        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.Assemble();

        this->mpLinearSystem->SwitchWriteModeLhsMatrix();
        PetscMatTools::Finalise(mMassMatrix);
    }


    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

    //////////////////////////////////////////
    // Set up z in b=Mz
    //////////////////////////////////////////
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();

    // get bidomain parameters
    double Am1 = this->mpExtendedBidomainTissue->GetAmFirstCell();
    double Am2 = this->mpExtendedBidomainTissue->GetAmSecondCell();
    double AmGap = this->mpExtendedBidomainTissue->GetAmGap();
    double Cm1 = this->mpExtendedBidomainTissue->GetCmFirstCell();
    double Cm2 = this->mpExtendedBidomainTissue->GetCmSecondCell();

    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    DistributedVector::Stripe distributed_current_solution_v_first_cell(distributed_current_solution, 0);
    DistributedVector::Stripe distributed_current_solution_v_second_cell(distributed_current_solution, 1);
    DistributedVector::Stripe distributed_current_solution_phi_e(distributed_current_solution, 2);

    // dist stripe for z
    DistributedVector dist_vec_matrix_based = p_factory->CreateDistributedVector(mVecForConstructingRhs);
    DistributedVector::Stripe dist_vec_matrix_based_v_first_cell(dist_vec_matrix_based, 0);
    DistributedVector::Stripe dist_vec_matrix_based_v_second_cell(dist_vec_matrix_based, 1);
    DistributedVector::Stripe dist_vec_matrix_based_phi_e(dist_vec_matrix_based, 2);

    for (DistributedVector::Iterator index = dist_vec_matrix_based.Begin();
         index!= dist_vec_matrix_based.End();
         ++index)
    {
        double V_first_cell = distributed_current_solution_v_first_cell[index];
        double V_second_Cell = distributed_current_solution_v_second_cell[index];

        double i_ionic_first_cell = this->mpExtendedBidomainTissue->rGetIionicCacheReplicated()[index.Global];
        double i_ionic_second_cell = this->mpExtendedBidomainTissue->rGetIionicCacheReplicatedSecondCell()[index.Global];
        double intracellular_stimulus_first_cell = this->mpExtendedBidomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];
        double intracellular_stimulus_second_cell = this->mpExtendedBidomainTissue->rGetIntracellularStimulusCacheReplicatedSecondCell()[index.Global];
        double extracellular_stimulus =  this->mpExtendedBidomainTissue->rGetExtracellularStimulusCacheReplicated()[index.Global];
        double g_gap = this->mpExtendedBidomainTissue->rGetGgapCacheReplicated()[index.Global];
        double delta_t = PdeSimulationTime::GetPdeTimeStep();
        dist_vec_matrix_based_v_first_cell[index] = Am1*Cm1*V_first_cell/delta_t  - Am1*i_ionic_first_cell + AmGap*g_gap*(V_second_Cell - V_first_cell) -  intracellular_stimulus_first_cell;
        dist_vec_matrix_based_v_second_cell[index] = Am2*Cm2*V_second_Cell/delta_t  - Am2*i_ionic_second_cell + AmGap*g_gap*(V_first_cell - V_second_Cell) - intracellular_stimulus_second_cell;

        if (this->mpExtendedBidomainTissue->HasTheUserSuppliedExtracellularStimulus() )
        {
            assert((fabs(intracellular_stimulus_first_cell) < 1e-12)
            && (fabs(intracellular_stimulus_second_cell) < 1e-12));///\todo turn these into exceptions somewhere else
            /**
             * The following line should also have
             *  - intracellular_stimulus_first_cell - intracellular_stimulus_second_cell in the summation,
             *  but they are zero...
             */
            dist_vec_matrix_based_phi_e[index] = -extracellular_stimulus;
        }
        else
        {
            dist_vec_matrix_based_phi_e[index] = 0.0;
        }
    }


    dist_vec_matrix_based.Restore();

    //////////////////////////////////////////
    // b = Mz
    //////////////////////////////////////////

    MatMult(mMassMatrix, mVecForConstructingRhs, this->mpLinearSystem->rGetRhsVector());

    // assembling RHS is not finished yet, as Neumann bcs are added below, but
    // the event will be begun again inside mpExtendedBidomainAssembler->AssembleVector();
    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);

    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////
    mpExtendedBidomainNeumannSurfaceTermAssembler->ResetBoundaryConditionsContainer(this->mpBoundaryConditions); // as the BCC can change
    mpExtendedBidomainNeumannSurfaceTermAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
    mpExtendedBidomainNeumannSurfaceTermAssembler->AssembleVector();

    this->mpLinearSystem->FinaliseRhsVector();

    this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

    if (computeMatrix)
    {
        this->mpLinearSystem->FinaliseLhsMatrix();
    }
    this->mpLinearSystem->FinaliseRhsVector();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::ExtendedBidomainSolver(
        bool bathSimulation,
        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
        ExtendedBidomainTissue<SPACE_DIM>* pTissue,
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,3>* pBoundaryConditions)
    : AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>(bathSimulation,pMesh,pTissue,pBoundaryConditions)
{
    // Tell Tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs = NULL;

    // create assembler
    if (this->mBathSimulation)
    {
        //this->mpExtendedExtendedBidomainAssembler = new ExtendedExtendedBidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpExtendedExtendedBidomainTissue,this->mDt);
        EXCEPTION("Bath simulations are not yet supported for extended bidomain problems");
    }
    else
    {
        mpExtendedBidomainAssembler = new ExtendedBidomainAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpExtendedBidomainTissue);
    }

    mpExtendedBidomainNeumannSurfaceTermAssembler = new ExtendedBidomainNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh,pBoundaryConditions);

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::~ExtendedBidomainSolver()
{
    delete mpExtendedBidomainAssembler;
    delete mpExtendedBidomainNeumannSurfaceTermAssembler;

    if (mVecForConstructingRhs)
    {
        PetscTools::Destroy(mVecForConstructingRhs);
        PetscTools::Destroy(mMassMatrix);
    }
}

// Explicit instantiation
template class ExtendedBidomainSolver<1,1>;
template class ExtendedBidomainSolver<2,2>;
template class ExtendedBidomainSolver<3,3>;
