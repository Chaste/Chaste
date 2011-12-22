/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1) //PETSc 3.1
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
 * If it's not available, it will throw the following error:
 *
 *     [0]PETSC ERROR: --------------------- Error Message ------------------------------------
 *     [0]PETSC ERROR: Unknown type. Check for miss-spelling or missing external package needed for type!
 *     [0]PETSC ERROR: Unable to find requested PC type hypre!
 *
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
