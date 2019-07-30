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

#ifndef PCTWOLEVELSBLOCKDIAGONAL_HPP_
#define PCTWOLEVELSBLOCKDIAGONAL_HPP_

#include <cassert>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include "PetscTools.hpp"

/**
 * PETSc will return the control to this function everytime it needs to precondition a vector (i.e. y = inv(M)*x)
 *
 * This function needs to be declared global, since I (migb) haven't found a way of defining it inside a class and
 * be able of passing it by reference.
 *
 * @param pc_context preconditioner context struct. Stores preconditioner state (i.e. PC, Mat, and Vec objects used)
 * @param x unpreconditioned residual.
 * @param y preconditioned residual. y = inv(M)*x
 */
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 1) //PETSc 3.1 or later
PetscErrorCode PCTwoLevelsBlockDiagonalApply(PC pc_context, Vec x, Vec y);
#else
PetscErrorCode PCTwoLevelsBlockDiagonalApply(void* pc_context, Vec x, Vec y);
#endif

/// \todo: #1082 update this description
/**
 *  This class defines a PETSc-compliant purpose-build preconditioner.
 *
 *  Let A be a matrix arising in the FEM discretisation of the bidomain equations with the following block structure:
 *
 *                 A = (A11  B')
 *                     (B   A22)
 *
 *  By creating an instance of this class, one will define the following preconditioner:
 *
 *                 inv(M) = inv( (A11   0)   = (inv(A11)        0)
 *                               (0   A22) )   (0        inv(A22))
 *
 *  The inverses are approximate with one cycle of AMG.
 *
 *  Note: This class requires PETSc to be build including HYPRE library. If it's not available, it will throw the following error:
 *
 *     [0]PETSC ERROR: --------------------- Error Message ------------------------------------
 *     [0]PETSC ERROR: Unknown type. Check for miss-spelling or missing external package needed for type!
 *     [0]PETSC ERROR: Unable to find requested PC type hypre!
 *
 *        and will approximate the inverse of the subblocks with PETSc's default preconditioner (bjacobi at the time of writing this).
 */
class PCTwoLevelsBlockDiagonal
{
public:

    /**
     * This struct defines the state of the preconditioner (initialised data and objects to be reused).
     */
    typedef struct{
        Mat A11_matrix_subblock; /**< Mat object that stores the A11 subblock*/
        Mat A22_B1_matrix_subblock; /**< Mat object that stores the tissue part of the A22 subblock*/
        Mat A22_B2_matrix_subblock; /**< Mat object that stores the bath part of the A22 subblock*/
        PC  PC_amg_A11; /**<  inv(A11) is approximated with ILU*/
        PC  PC_amg_A22_B1; /**<  inv(A22_B1) is approximated with ILU*/
        PC  PC_amg_A22_B2; /**<  inv(A22_B2) is approximated by an AMG cycle. We compute it with HYPRE via a PC object*/
        Vec x1_subvector;/**<  Used to store the first half of the vector to be preconditioned*/
        Vec x21_subvector;/**<  Used to store the tissue part of the second half of the vector to be preconditioned*/
        Vec x22_subvector;/**<  Used to store the bath part of the second half of the vector to be preconditioned*/
        Vec y1_subvector;/**<  Used to store the first half of the preconditioned vector*/
        Vec y21_subvector;/**<  Used to store the tissue part of the second half of the preconditioned vector*/
        Vec y22_subvector;/**<  Used to store the bath part second half of the preconditioned vector*/
        VecScatter A11_scatter_ctx;/**< Scattering context: gather x1 from x and scatter y1 back into y*/
        VecScatter A22_B1_scatter_ctx;/**< Scattering context: gather x21 from x and scatter y21 back into y*/
        VecScatter A22_B2_scatter_ctx;/**< Scattering context: gather x22 from x and scatter y22 back into y*/

    } PCTwoLevelsBlockDiagonalContext;

    PCTwoLevelsBlockDiagonalContext mPCContext; /**< PC context, this will be passed to PCTwoLevelsBlockDiagonalApply when PETSc returns control to our preconditioner subroutine.  See PCShellSetContext().*/
    PC mPetscPCObject;/**< Generic PETSc preconditioner object */

public:

    /**
     * Constructor.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     * @param rBathNodes a list of nodes defining the bath
     */
    PCTwoLevelsBlockDiagonal(KSP& rKspObject, std::vector<PetscInt>& rBathNodes);

    /**
     * Destructor.
     */
    ~PCTwoLevelsBlockDiagonal();

private:

    /**
     * Creates all the state data required by the preconditioner.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     * @param rBathNodes a list of nodes defining the bath
     */
    void PCTwoLevelsBlockDiagonalCreate(KSP& rKspObject, std::vector<PetscInt>& rBathNodes);

    /**
     * Setups preconditioner.
     */
    void PCTwoLevelsBlockDiagonalSetUp();
};

#endif /*PCTWOLEVELSBLOCKDIAGONAL_HPP_*/
