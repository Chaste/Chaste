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

#include "PCTwoLevelsBlockDiagonal.hpp"
#include "Exception.hpp"

#include <iostream>

PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonal(KSP& rKspObject, std::vector<PetscInt>& rBathNodes)
{
    PCTwoLevelsBlockDiagonalCreate(rKspObject, rBathNodes);
    PCTwoLevelsBlockDiagonalSetUp();
}

PCTwoLevelsBlockDiagonal::~PCTwoLevelsBlockDiagonal()
{
    PetscTools::Destroy(mPCContext.A11_matrix_subblock);
    PetscTools::Destroy(mPCContext.A22_B1_matrix_subblock);
    PetscTools::Destroy(mPCContext.A22_B2_matrix_subblock);

    PCDestroy(PETSC_DESTROY_PARAM(mPCContext.PC_amg_A11));
    PCDestroy(PETSC_DESTROY_PARAM(mPCContext.PC_amg_A22_B1));
    PCDestroy(PETSC_DESTROY_PARAM(mPCContext.PC_amg_A22_B2));

    PetscTools::Destroy(mPCContext.x1_subvector);
    PetscTools::Destroy(mPCContext.y1_subvector);

    PetscTools::Destroy(mPCContext.x21_subvector);
    PetscTools::Destroy(mPCContext.y21_subvector);

    PetscTools::Destroy(mPCContext.x22_subvector);
    PetscTools::Destroy(mPCContext.y22_subvector);

    VecScatterDestroy(PETSC_DESTROY_PARAM(mPCContext.A11_scatter_ctx));
    VecScatterDestroy(PETSC_DESTROY_PARAM(mPCContext.A22_B1_scatter_ctx));
    VecScatterDestroy(PETSC_DESTROY_PARAM(mPCContext.A22_B2_scatter_ctx));
}

void PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalCreate(KSP& rKspObject, std::vector<PetscInt>& rBathNodes)
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
    assert(num_rows==num_columns);

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

    unsigned subvector_num_rows_tissue = subvector_num_rows - rBathNodes.size();
    unsigned subvector_local_rows_tissue = subvector_num_rows_tissue; /// \todo: #1082 won't work in parallel

    unsigned subvector_num_rows_bath = rBathNodes.size();
    unsigned subvector_local_rows_bath = subvector_num_rows_bath; /// \todo: #1082 won't work in parallel

    assert(PetscTools::IsSequential());

    mPCContext.x1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.x21_subvector = PetscTools::CreateVec(subvector_num_rows_tissue, subvector_local_rows_tissue);
    mPCContext.x22_subvector = PetscTools::CreateVec(subvector_num_rows_bath, subvector_local_rows_bath);
    mPCContext.y1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.y21_subvector = PetscTools::CreateVec(subvector_num_rows_tissue, subvector_local_rows_tissue);
    mPCContext.y22_subvector = PetscTools::CreateVec(subvector_num_rows_bath, subvector_local_rows_bath);

    // Define IS objects that will be used throughout the method
    IS A11_all_rows;
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 2, &A11_all_rows);

    IS A22_all_rows;
    PetscInt A22_size;
    VecGetSize(mPCContext.x1_subvector, &A22_size);
    /// \todo: #1082 assert size(x1) = size(x21) + size(x22)
    ISCreateStride(PETSC_COMM_WORLD, A22_size, 1, 2, &A22_all_rows);

    IS A22_bath_rows;
    PetscInt* phi_e_bath_rows = new PetscInt[rBathNodes.size()];
    for (unsigned index=0; index<rBathNodes.size(); index++)
    {
        phi_e_bath_rows[index] = 2*rBathNodes[index] + 1;
    }
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
    ISCreateGeneral(PETSC_COMM_WORLD, rBathNodes.size(), phi_e_bath_rows, PETSC_USE_POINTER, &A22_bath_rows);
 #else
    ISCreateGeneralWithArray(PETSC_COMM_WORLD, rBathNodes.size(), phi_e_bath_rows, &A22_bath_rows);
#endif

    IS A22_tissue_rows;
    ISDifference(A22_all_rows, A22_bath_rows, &A22_tissue_rows);

    // Create scatter contexts
    {
        // Needed by VecScatterCreate in order to find out parallel layout.
        Vec dummy_vec = PetscTools::CreateVec(num_rows, num_local_rows);

        /// \todo: #1082 legacy, no need to use the references
        IS& A11_rows=A11_all_rows;
        IS& A22_B1_rows=A22_tissue_rows;
        IS& A22_B2_rows=A22_bath_rows;

        IS all_vector;
        ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 1, &all_vector);

        IS tissue_vector;
        ISCreateStride(PETSC_COMM_WORLD, (num_rows/2)-rBathNodes.size(), 0, 1, &tissue_vector);

        IS bath_vector;
        ISCreateStride(PETSC_COMM_WORLD, rBathNodes.size(), 0, 1, &bath_vector);

        VecScatterCreate(dummy_vec, A11_rows, mPCContext.x1_subvector, all_vector, &mPCContext.A11_scatter_ctx);
        VecScatterCreate(dummy_vec, A22_B1_rows, mPCContext.x21_subvector, tissue_vector, &mPCContext.A22_B1_scatter_ctx);
        VecScatterCreate(dummy_vec, A22_B2_rows, mPCContext.x22_subvector, bath_vector, &mPCContext.A22_B2_scatter_ctx);

        ISDestroy(PETSC_DESTROY_PARAM(all_vector));
        ISDestroy(PETSC_DESTROY_PARAM(tissue_vector));
        ISDestroy(PETSC_DESTROY_PARAM(bath_vector));

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
        IS& A11_columns=A11_all_rows;
        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low, 2, &A11_local_rows); /// \todo: #1082 OK in parallel. Use as an example for the other two blocks

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 8) //PETSc 3.8 or later
        MatCreateSubMatrix(system_matrix, A11_local_rows, A11_columns,
                           MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A11_local_rows, A11_columns,
                        MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#endif

        ISDestroy(PETSC_DESTROY_PARAM(A11_local_rows));
    }

    // Get matrix sublock A22_B1
    {
//        // Work out local row range for subblock A22 (same as x2 or y2)
//        PetscInt low, high, global_size;
//        VecGetOwnershipRange(mPCContext.x21_subvector, &low, &high);
//        VecGetSize(mPCContext.x21_subvector, &global_size);
//        assert(global_size == (num_rows/2) - (PetscInt) rBathNodes.size());

        assert(PetscTools::IsSequential());
        IS& A22_B1_local_rows = A22_tissue_rows; // wrong in parallel, need to give local rows
        IS& A22_B1_columns = A22_tissue_rows;

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 8) //PETSc 3.8 or later
        MatCreateSubMatrix(system_matrix, A22_B1_local_rows, A22_B1_columns,
                           MAT_INITIAL_MATRIX, &mPCContext.A22_B1_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A22_B1_local_rows, A22_B1_columns,
                        MAT_INITIAL_MATRIX, &mPCContext.A22_B1_matrix_subblock);
#endif

    }

    // Get matrix sublock A22_B2
    {
//        // Work out local row range for subblock A22 (same as x2 or y2)
//        PetscInt low, high, global_size;
//        VecGetOwnershipRange(mPCContext.x21_subvector, &low, &high);
//        VecGetSize(mPCContext.x21_subvector, &global_size);
//        assert(global_size == (num_rows/2) - (PetscInt) rBathNodes.size());

        assert(PetscTools::IsSequential());
        IS& A22_B2_local_rows = A22_bath_rows; // wrong in parallel, need to give local rows
        IS& A22_B2_columns = A22_bath_rows;

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 8) //PETSc 3.8 or later
        MatCreateSubMatrix(system_matrix, A22_B2_local_rows, A22_B2_columns,
                           MAT_INITIAL_MATRIX, &mPCContext.A22_B2_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A22_B2_local_rows, A22_B2_columns,
                        MAT_INITIAL_MATRIX, &mPCContext.A22_B2_matrix_subblock);
#endif
    }

    ISDestroy(PETSC_DESTROY_PARAM(A11_all_rows));
    ISDestroy(PETSC_DESTROY_PARAM(A22_all_rows));
    ISDestroy(PETSC_DESTROY_PARAM(A22_bath_rows));
    delete[] phi_e_bath_rows;
    ISDestroy(PETSC_DESTROY_PARAM(A22_tissue_rows));

    // Register call-back function and its context
    PCSetType(mPetscPCObject, PCSHELL);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PCShellSetApply(mPetscPCObject, PCTwoLevelsBlockDiagonalApply, (void*) &mPCContext);
#else
    // Register PC context so it gets passed to PCTwoLevelsBlockDiagonalApply
    PCShellSetContext(mPetscPCObject, &mPCContext);

    // Register call-back function
    PCShellSetApply(mPetscPCObject, PCTwoLevelsBlockDiagonalApply);
#endif
}

void PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalSetUp()
{
    // These options will get read by PCSetFromOptions
    PetscTools::SetOption("-pc_hypre_boomeramg_max_iter", "1");
    PetscTools::SetOption("-pc_hypre_boomeramg_strong_threshold", "0.0");
    PetscTools::SetOption("-pc_hypre_type", "boomeramg");

    // Set up amg preconditioner for block A11
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A11));
    PCSetType(mPCContext.PC_amg_A11, PCBJACOBI);
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
    PCSetOperators(mPCContext.PC_amg_A11, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock);
#else
    PCSetOperators(mPCContext.PC_amg_A11, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
#endif
    PCSetFromOptions(mPCContext.PC_amg_A11);
    PCSetUp(mPCContext.PC_amg_A11);

    // Set up amg preconditioner for block A22_B1
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A22_B1));
    PCSetType(mPCContext.PC_amg_A22_B1, PCBJACOBI);
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
    PCSetOperators(mPCContext.PC_amg_A22_B1, mPCContext.A22_B1_matrix_subblock, mPCContext.A22_B1_matrix_subblock);
#else
    PCSetOperators(mPCContext.PC_amg_A22_B1, mPCContext.A22_B1_matrix_subblock, mPCContext.A22_B1_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
#endif
    PCSetFromOptions(mPCContext.PC_amg_A22_B1);
    PCSetUp(mPCContext.PC_amg_A22_B1);

    // Set up amg preconditioner for block A22_B2
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A22_B2));
    PCSetType(mPCContext.PC_amg_A22_B2, PCHYPRE);
    //PCHYPRESetType(mPCContext.PC_amg_A22_B2, "boomeramg");
    PetscTools::SetOption("-pc_hypre_type", "boomeramg");

    PetscTools::SetOption("-pc_hypre_boomeramg_max_iter", "1");
    PetscTools::SetOption("-pc_hypre_boomeramg_strong_threshold", "0.0");
    PetscTools::SetOption("-pc_hypre_boomeramg_coarsen_type", "HMIS");

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5)
    PCSetOperators(mPCContext.PC_amg_A22_B2, mPCContext.A22_B2_matrix_subblock, mPCContext.A22_B2_matrix_subblock);
#else
    PCSetOperators(mPCContext.PC_amg_A22_B2, mPCContext.A22_B2_matrix_subblock, mPCContext.A22_B2_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
#endif
    PCSetFromOptions(mPCContext.PC_amg_A22_B2);
    PCSetUp(mPCContext.PC_amg_A22_B2);
}

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 1) //PETSc 3.1 or later
PetscErrorCode PCTwoLevelsBlockDiagonalApply(PC pc_object, Vec x, Vec y)
{
  void* pc_context;

  PCShellGetContext(pc_object, &pc_context);
#else
PetscErrorCode PCTwoLevelsBlockDiagonalApply(void* pc_context, Vec x, Vec y)
{
#endif

    // Cast the context pointer to PCTwoLevelsBlockDiagonalContext
    PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalContext* block_diag_context = (PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalContext*) pc_context;
    assert(block_diag_context!=nullptr);

    /*
     * Scatter x = [x1 x21 x22]'
     */
//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A11_scatter_ctx, x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A11_scatter_ctx, x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A11_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A11_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B1_scatter_ctx, x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A22_B1_scatter_ctx, x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B1_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B1_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B2_scatter_ctx, x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A22_B2_scatter_ctx, x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B2_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B2_scatter_ctx);
#endif

    /*
     * y1  = ILU(A11)*x1
     * y21 = ILU(A22)*x21
     * y22 = AMG(A22)*x22
     */
    PCApply(block_diag_context->PC_amg_A11, block_diag_context->x1_subvector, block_diag_context->y1_subvector);
    PCApply(block_diag_context->PC_amg_A22_B1, block_diag_context->x21_subvector, block_diag_context->y21_subvector);
    PCApply(block_diag_context->PC_amg_A22_B2, block_diag_context->x22_subvector, block_diag_context->y22_subvector);

    /*
     * Gather y = [y1 y21 y22]'
     */
//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A11_scatter_ctx, block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A11_scatter_ctx, block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A11_scatter_ctx);
    VecScatterEnd(block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A11_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B1_scatter_ctx, block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A22_B1_scatter_ctx, block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B1_scatter_ctx);
    VecScatterEnd(block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B1_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B2_scatter_ctx, block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A22_B2_scatter_ctx, block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B2_scatter_ctx);
    VecScatterEnd(block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B2_scatter_ctx);
#endif

    return 0;
}
