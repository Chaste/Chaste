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

#include <iostream>

#include "PetscVecTools.hpp" // Includes Ublas so must come first
#include "PCLDUFactorisation.hpp"
#include "Exception.hpp"
#include "Warnings.hpp"

#ifdef TRACE_KSP
#include "Timer.hpp"
#endif

PCLDUFactorisation::PCLDUFactorisation(KSP& rKspObject)
{
#ifdef TRACE_KSP
    mPCContext.mScatterTime = 0.0;
    mPCContext.mA1PreconditionerTime = 0.0;;
    mPCContext.mA2PreconditionerTime = 0.0;;
    mPCContext.mGatherTime = 0.0;;
#endif

    PCLDUFactorisationCreate(rKspObject);
    PCLDUFactorisationSetUp();
}

PCLDUFactorisation::~PCLDUFactorisation()
{
#ifdef TRACE_KSP
    if (PetscTools::AmMaster())
    {
        std::cout << " -- LDU factorisation preconditioner profile information: " << std::endl;
        std::cout << "\t mScatterTime: " << mPCContext.mScatterTime << std::endl;
        std::cout << "\t mA1PreconditionerTime: " << mPCContext.mA1PreconditionerTime << std::endl;
        std::cout << "\t mA2PreconditionerTime: " << mPCContext.mA2PreconditionerTime << std::endl;
        std::cout << "\t mExtraLAOperations: " << mPCContext.mExtraLAOperations << std::endl;
        std::cout << "\t mGatherTime: " << mPCContext.mGatherTime << std::endl;
    }
#endif

    PetscTools::Destroy(mPCContext.A11_matrix_subblock);
    PetscTools::Destroy(mPCContext.A22_matrix_subblock);
    PetscTools::Destroy(mPCContext.B_matrix_subblock);

    PCDestroy(PETSC_DESTROY_PARAM(mPCContext.PC_amg_A11));
    PCDestroy(PETSC_DESTROY_PARAM(mPCContext.PC_amg_A22));

    PetscTools::Destroy(mPCContext.x1_subvector);
    PetscTools::Destroy(mPCContext.y1_subvector);
    PetscTools::Destroy(mPCContext.x2_subvector);
    PetscTools::Destroy(mPCContext.y2_subvector);
    PetscTools::Destroy(mPCContext.z);
    PetscTools::Destroy(mPCContext.temp);

    VecScatterDestroy(PETSC_DESTROY_PARAM(mPCContext.A11_scatter_ctx));
    VecScatterDestroy(PETSC_DESTROY_PARAM(mPCContext.A22_scatter_ctx));
}

void PCLDUFactorisation::PCLDUFactorisationCreate(KSP& rKspObject)
{
    KSPGetPC(rKspObject, &mPetscPCObject);

    Mat system_matrix, dummy;
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
    KSPGetOperators(rKspObject, &system_matrix, &dummy);
#else
    MatStructure flag;
    KSPGetOperators(rKspObject, &system_matrix, &dummy, &flag);
#endif

    PetscInt num_rows, num_columns;
    MatGetSize(system_matrix, &num_rows, &num_columns);

    PetscInt num_local_rows, num_local_columns;
    MatGetLocalSize(system_matrix, &num_local_rows, &num_local_columns);

    // Odd number of rows: impossible in Bidomain.
    // Odd number of local rows: impossible if V_m and phi_e for each node are stored in the same processor.
    if ((num_rows%2 != 0) || (num_local_rows%2 != 0))
    {
        TERMINATE("Wrong matrix parallel layout detected in PCLDUFactorisation."); // LCOV_EXCL_LINE
    }

    // Allocate memory
    unsigned subvector_num_rows = num_rows/2;
    unsigned subvector_local_rows = num_local_rows/2;
    mPCContext.x1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.x2_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.y1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.y2_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.z = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.temp = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);

    // Create scatter contexts
    {
        // Needed by VecScatterCreate in order to find out parallel layout.
        Vec dummy_vec = PetscTools::CreateVec(num_rows, num_local_rows);

        PetscVecTools::SetupInterleavedVectorScatterGather(dummy_vec, mPCContext.A11_scatter_ctx, mPCContext.A22_scatter_ctx);

        PetscTools::Destroy(dummy_vec);
    }

    // Get matrix sublock A11
    {
        // Work out local row range for subblock A11 (same as x1 or y1)
        PetscInt low, high, global_size;
        VecGetOwnershipRange(mPCContext.x1_subvector, &low, &high);
        VecGetSize(mPCContext.x1_subvector, &global_size);
        assert(global_size == num_rows/2);

        IS A11_local_rows;
        IS A11_columns;
        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low, 2, &A11_local_rows);
        ISCreateStride(PETSC_COMM_WORLD, global_size, 0, 2, &A11_columns);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 8) //PETSc 3.8 or later
        MatCreateSubMatrix(system_matrix, A11_local_rows, A11_local_rows,
                           MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A11_local_rows, A11_local_rows,
                        MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#endif
        ISDestroy(PETSC_DESTROY_PARAM(A11_local_rows));
        ISDestroy(PETSC_DESTROY_PARAM(A11_columns));
    }

    // Get matrix sublock A22
    {
        // Work out local row range for subblock A22 (same as x2 or y2)
        PetscInt low, high, global_size;
        VecGetOwnershipRange(mPCContext.x2_subvector, &low, &high);
        VecGetSize(mPCContext.x2_subvector, &global_size);
        assert(global_size == num_rows/2);

        IS A22_local_rows;
        IS A22_columns;
        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low+1, 2, &A22_local_rows);
        ISCreateStride(PETSC_COMM_WORLD, global_size, 1, 2, &A22_columns);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 8) //PETSc 3.8 or later
        MatCreateSubMatrix(system_matrix, A22_local_rows, A22_local_rows,
                           MAT_INITIAL_MATRIX, &mPCContext.A22_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A22_local_rows, A22_local_rows,
                        MAT_INITIAL_MATRIX, &mPCContext.A22_matrix_subblock);
#endif

        ISDestroy(PETSC_DESTROY_PARAM(A22_local_rows));
        ISDestroy(PETSC_DESTROY_PARAM(A22_columns));
    }

    // Get matrix sublock B (the upper triangular one)
    {
        // Work out local row range for subblock B (same as A11, x1 or y1)
        PetscInt low, high, global_size;
        VecGetOwnershipRange(mPCContext.x1_subvector, &low, &high);
        VecGetSize(mPCContext.x1_subvector, &global_size);
        assert(global_size == num_rows/2);

        IS B_local_rows;
        IS B_columns;
        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low, 2, &B_local_rows);

        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low+1, 2, &B_columns);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 8) //PETSc 3.8 or later
        MatCreateSubMatrix(system_matrix, B_local_rows, B_columns,
                           MAT_INITIAL_MATRIX, &mPCContext.B_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, B_local_rows, B_columns,
                        MAT_INITIAL_MATRIX, &mPCContext.B_matrix_subblock);
#endif

        ISDestroy(PETSC_DESTROY_PARAM(B_local_rows));
        ISDestroy(PETSC_DESTROY_PARAM(B_columns));
    }

    /*
     * Experimental (#1082): in PP removing the mass matrix from the A22 block seems to work better.
     *                       This is equivalent to do A22 = A22 + B in this implementation.
     */
//#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
//     PetscScalar petsc_one = 1.0;
//     MatAXPY(&petsc_one, mPCContext.B_matrix_subblock, mPCContext.A22_matrix_subblock, DIFFERENT_NONZERO_PATTERN);
// #else
//     MatAXPY(mPCContext.A22_matrix_subblock, 1.0, mPCContext.B_matrix_subblock, DIFFERENT_NONZERO_PATTERN);
//#endif

//     // Shift the block
//#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
//     PetscScalar shift = -1e-8;
//     MatShift(&shift, mPCContext.A22_matrix_subblock);
// #else
//     PetscScalar shift = -1e-8;
//     MatShift(mPCContext.A22_matrix_subblock, shift);
//#endif

    PCSetType(mPetscPCObject, PCSHELL);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    // Register PC context and call-back function
    PCShellSetApply(mPetscPCObject, PCLDUFactorisationApply, (void*) &mPCContext);
#else
    // Register PC context so it gets passed to PCBlockDiagonalApply
    PCShellSetContext(mPetscPCObject, &mPCContext);
    // Register call-back function
    PCShellSetApply(mPetscPCObject, PCLDUFactorisationApply);
#endif
}

void PCLDUFactorisation::PCLDUFactorisationSetUp()
{
    // These options will get read by PCSetFromOptions
//     PetscTools::SetOption("-pc_hypre_boomeramg_max_iter", "1");
//     PetscTools::SetOption("-pc_hypre_boomeramg_strong_threshold", "0.0");
//     PetscTools::SetOption("-pc_hypre_type", "boomeramg");

    /*
     * Set up preconditioner for block A11
     */
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A11));
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
    // Attempt to emulate SAME_PRECONDITIONER below
    PCSetReusePreconditioner(mPCContext.PC_amg_A11, PETSC_TRUE);
    PCSetOperators(mPCContext.PC_amg_A11, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock);
#else
    PCSetOperators(mPCContext.PC_amg_A11, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock, SAME_PRECONDITIONER);
#endif

    // Choose between the two following blocks in order to approximate inv(A11) with one AMG cycle
    // or with an CG solve with high tolerance
////////
//    PCSetType(mPCContext.PC_amg_A11, PCBJACOBI);

    // We are expecting an error from PETSC on systems that don't have the hypre library, so suppress it
    // in case it aborts
    PetscPushErrorHandler(PetscIgnoreErrorHandler, nullptr);
    PetscErrorCode pc_set_error = PCSetType(mPCContext.PC_amg_A11, PCHYPRE);
    if (pc_set_error != 0)
    {
        WARNING("PETSc hypre preconditioning library is not installed");
    }
    // Stop supressing error
    PetscPopErrorHandler();

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=5)
    PetscTools::SetOption("-pc_hypre_type", "euclid");
    PetscTools::SetOption("-pc_hypre_euclid_levels", "0");
#else
    /* euclid was removed in 3.6. Use the previously-commented-out alternative */
    PetscTools::SetOption("-pc_hypre_type", "boomeramg");
    PetscTools::SetOption("-pc_hypre_boomeramg_max_iter", "1");
    PetscTools::SetOption("-pc_hypre_boomeramg_strong_threshold", "0.0");
    PetscTools::SetOption("-pc_hypre_boomeramg_coarsen_type", "HMIS");
#endif

////////
//    PCSetType(mPCContext.PC_amg_A11, PCKSP);
//    KSP ksp1;
//    PCKSPGetKSP(mPCContext.PC_amg_A11,&ksp1);
//    KSPSetType(ksp1, KSPCG);
//    KSPSetTolerances(ksp1, 0.1, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
//
//    PC prec1;
//    KSPGetPC(ksp1, &prec1);
//    PCSetType(prec1, PCBJACOBI);
//    PCSetFromOptions(prec1);
//    PCSetOperators(prec1, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
//    PCSetUp(prec1);
//
//    KSPSetFromOptions(ksp1);
//    KSPSetUp(ksp1);
////////

    PCSetFromOptions(mPCContext.PC_amg_A11);
    PCSetUp(mPCContext.PC_amg_A11);

    /*
     * Set up amg preconditioner for block A22
     */
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A22));
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
    // Attempt to emulate SAME_PRECONDITIONER below
    PCSetReusePreconditioner(mPCContext.PC_amg_A22, PETSC_TRUE);
    PCSetOperators(mPCContext.PC_amg_A22, mPCContext.A22_matrix_subblock, mPCContext.A22_matrix_subblock);
#else
    PCSetOperators(mPCContext.PC_amg_A22, mPCContext.A22_matrix_subblock, mPCContext.A22_matrix_subblock, SAME_PRECONDITIONER);
#endif

    // Choose between the two following blocks in order to approximate inv(A11) with one AMG cycle
    // or with an CG solve with high tolerance
////////
    // We are expecting an error from PETSC on systems that don't have the hypre library, so suppress it
    // in case it aborts
    PetscPushErrorHandler(PetscIgnoreErrorHandler, nullptr);
    PCSetType(mPCContext.PC_amg_A22, PCHYPRE);
    // Stop supressing error
    PetscPopErrorHandler();

    PetscTools::SetOption("-pc_hypre_type", "boomeramg");
    PetscTools::SetOption("-pc_hypre_boomeramg_max_iter", "1");
    PetscTools::SetOption("-pc_hypre_boomeramg_strong_threshold", "0.0");
    PetscTools::SetOption("-pc_hypre_boomeramg_coarsen_type", "HMIS");
    //    PetscTools::SetOption("-pc_hypre_boomeramg_interp_type","ext+i");

////////
//    PCSetType(mPCContext.PC_amg_A22, PCKSP);
//    KSP ksp2;
//    PCKSPGetKSP(mPCContext.PC_amg_A22,&ksp2);
//    KSPSetType(ksp2, KSPCG);
//    KSPSetTolerances(ksp2, 0.1, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
//
//    PC prec2;
//    KSPGetPC(ksp2, &prec2);
//    PCSetType(prec2, PCBJACOBI);
//    PCSetFromOptions(prec2);
//    PCSetOperators(prec2, mPCContext.A22_matrix_subblock, mPCContext.A22_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
//    PCSetUp(prec2);
//
//    KSPSetFromOptions(ksp2);
//    KSPSetUp(ksp2);
////////

    PCSetFromOptions(mPCContext.PC_amg_A22);
    PCSetUp(mPCContext.PC_amg_A22);
}
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 1) //PETSc 3.1 or later
PetscErrorCode PCLDUFactorisationApply(PC pc_object, Vec x, Vec y)
{
  void* pc_context;

  PCShellGetContext(pc_object, &pc_context);
#else
PetscErrorCode PCLDUFactorisationApply(void* pc_context, Vec x, Vec y)
{
#endif
    /// \todo refactoring: create a method for scattering and another for reversing

    // Cast the pointer to a PC context to our defined type
    PCLDUFactorisation::PCLDUFactorisationContext* block_diag_context = (PCLDUFactorisation::PCLDUFactorisationContext*) pc_context;
    assert(block_diag_context!=nullptr);

    /*
     * Split vector x into two. x = [x1 x2]'
     */
#ifdef TRACE_KSP
    Timer::Reset();
#endif

    PetscVecTools::DoInterleavedVecScatter(x, block_diag_context->A11_scatter_ctx, block_diag_context->x1_subvector, block_diag_context->A22_scatter_ctx, block_diag_context->x2_subvector);

#ifdef TRACE_KSP
    block_diag_context->mScatterTime += Timer::GetElapsedTime();
#endif

    /*
     * Apply preconditioner: [y1 y2]' = inv(P)[x1 x2]'
     *
     *    z  = inv(A11)*x1
     *    y2 = inv(A22)*(x2 - B*z)
     *    y1 = z - inv(A11)(B*y2)
     */
#ifdef TRACE_KSP
    Timer::Reset();
#endif
    //z  = inv(A11)*x1
    PCApply(block_diag_context->PC_amg_A11, block_diag_context->x1_subvector, block_diag_context->z);
#ifdef TRACE_KSP
    block_diag_context->mA1PreconditionerTime += Timer::GetElapsedTime();
#endif

    //y2 = inv(A22)*(x2 - B*z)
#ifdef TRACE_KSP
    Timer::Reset();
#endif
    MatMult(block_diag_context->B_matrix_subblock,block_diag_context->z,block_diag_context->temp); //temp = B*z
    double minus_one = -1.0;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    VecAYPX(&minus_one, block_diag_context->x2_subvector, block_diag_context->temp); // temp <-- x2 - temp
#else
    VecAYPX(block_diag_context->temp, minus_one, block_diag_context->x2_subvector); // temp <-- x2 - temp
#endif
#ifdef TRACE_KSP
    block_diag_context->mExtraLAOperations += Timer::GetElapsedTime();
#endif


#ifdef TRACE_KSP
    Timer::Reset();
#endif
    PCApply(block_diag_context->PC_amg_A22, block_diag_context->temp, block_diag_context->y2_subvector); // y2 = inv(A22)*temp
#ifdef TRACE_KSP
    block_diag_context->mA2PreconditionerTime += Timer::GetElapsedTime();
#endif

    // y1 = z - inv(A11)(B*y2)
#ifdef TRACE_KSP
    Timer::Reset();
#endif
    MatMult(block_diag_context->B_matrix_subblock,block_diag_context->y2_subvector,block_diag_context->temp); //temp = B*y2
#ifdef TRACE_KSP
    block_diag_context->mExtraLAOperations += Timer::GetElapsedTime();
#endif
#ifdef TRACE_KSP
    Timer::Reset();
#endif
    PCApply(block_diag_context->PC_amg_A11, block_diag_context->temp, block_diag_context->y1_subvector); // y1 = inv(A11)*temp
#ifdef TRACE_KSP
    block_diag_context->mA1PreconditionerTime += Timer::GetElapsedTime();
#endif

#ifdef TRACE_KSP
    Timer::Reset();
#endif
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    VecAYPX(&minus_one, block_diag_context->z, block_diag_context->y1_subvector); // y1 <-- z - y1
#else
    VecAYPX(block_diag_context->y1_subvector, minus_one, block_diag_context->z); // y1 <-- z - y1
#endif
#ifdef TRACE_KSP
    block_diag_context->mExtraLAOperations += Timer::GetElapsedTime();
#endif

    /*
     * Gather vectors y1 and y2. y = [y1 y2]'
     */
#ifdef TRACE_KSP
    Timer::Reset();
#endif

    PetscVecTools::DoInterleavedVecGather(y, block_diag_context->A11_scatter_ctx, block_diag_context->y1_subvector, block_diag_context->A22_scatter_ctx, block_diag_context->y2_subvector);

#ifdef TRACE_KSP
    block_diag_context->mGatherTime += Timer::GetElapsedTime();
#endif

    return 0;
}
