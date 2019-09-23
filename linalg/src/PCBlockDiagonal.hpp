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

#ifndef PCBLOCKDIAGONAL_HPP_
#define PCBLOCKDIAGONAL_HPP_

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
PetscErrorCode PCBlockDiagonalApply(PC pc_context, Vec x, Vec y);
#else
PetscErrorCode PCBlockDiagonalApply(void* pc_context, Vec x, Vec y);
#endif

/**
 * This class defines a PETSc-compliant purpouse-build preconditioner.
 *
 * Let A be a matrix arising in the FEM discretisation of the bidomain
 * equations with the following block structure:
 *
 *                 A = (A11  B')
 *                     (B   A22)
 *
 * By creating an instance of this class, one will define the following
 * preconditioner:
 *
 *                 inv(M) = inv( (A11   0)   = (inv(A11)        0)
 *                               (0   A22) )   (0        inv(A22))
 *
 * The inverses are approximate with one cycle of AMG.
 *
 * Note: This class requires PETSc to be build including HYPRE library.
 * If it's not available, it will show the following warning:
   Chaste warning: in file linalg/src/PCLDUFactorisation.cpp at line ???: PETSc HYPRE preconditioning library is not installed
 * and will approximate the inverse of the subblocks with PETSc's default
 * preconditioner (bjacobi at the time of writing this).
 */
class PCBlockDiagonal
{
public:

    /**
     * This struct defines the state of the preconditioner (initialised data and objects to be reused)
     */
    typedef struct{
        Mat A11_matrix_subblock; /**< Mat object that stores the A11 subblock*/
        Mat A22_matrix_subblock; /**< Mat object that stores the A22 subblock*/
        PC  PC_amg_A11; /**<  inv(A11) is approximated by an AMG cycle. We compute it with HYPRE via a PC object*/
        PC  PC_amg_A22; /**<  inv(A22) is approximated by an AMG cycle. We compute it with HYPRE via a PC object*/
        Vec x1_subvector;/**<  Used to store the first half of the vector to be preconditioned*/
        Vec x2_subvector;/**<  Used to store the second half of the vector to be preconditioned*/
        Vec y1_subvector;/**<  Used to store the first half of the preconditioned vector*/
        Vec y2_subvector;/**<  Used to store the second half of the preconditioned vector*/
        VecScatter A11_scatter_ctx;/**< Scattering context: gather x1 from x and scatter y1 back into y*/
        VecScatter A22_scatter_ctx;/**< Scattering context: gather x2 from x and scatter y2 back into y*/
#ifdef TRACE_KSP
        double mScatterTime;/**< Time counter used for profiling scatter operations*/
        double mA1PreconditionerTime;/**< Time counter used for profiling the application of the preconditioner on the A11 block*/
        double mA2PreconditionerTime;/**< Time counter used for profiling the application of the preconditioner on the A22 block*/
        double mGatherTime;/**< Time counter used for profiling gather operations*/
#endif

    } PCBlockDiagonalContext;

    PCBlockDiagonalContext mPCContext; /**< PC context, this will be passed to PCBlockDiagonalApply when PETSc returns control to our preconditioner subroutine.  See PCShellSetContext().*/
    PC mPetscPCObject;/**< Generic PETSc preconditioner object */

    /**
     * Constructor.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     */
    PCBlockDiagonal(KSP& rKspObject);

    ~PCBlockDiagonal();

private:

    /**
     * Creates all the state data required by the preconditioner.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     */
    void PCBlockDiagonalCreate(KSP& rKspObject);

    /**
     * Setups preconditioner.
     */
    void PCBlockDiagonalSetUp();
};

#endif /*PCBLOCKDIAGONAL_HPP_*/
