/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef PCLDUFACTORISATION_HPP_
#define PCLDUFACTORISATION_HPP_

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
#if (PETSC_VERSION_MAJOR == 3 && (PETSC_VERSION_MINOR == 1 || PETSC_VERSION_MINOR ==  2) ) //PETSc 3.1
PetscErrorCode PCLDUFactorisationApply(PC pc_object, Vec x, Vec y);
#else
PetscErrorCode PCLDUFactorisationApply(void* pc_context, Vec x, Vec y);
#endif

/**
 * This class defines a PETSc-compliant purpose-built preconditioner.
 *
 * Let A be a matrix arising in the FEM discretisation of the bidomain
 * equations with the following block structure:
 *
 *                 A = (A11  B')
 *                     (B  A22)
 *
 * Let A=LDU be the following matrix factorisation
 *
 *                 LDU = (I      0)(A11 0)(I inv(A11)B')
 *                       (B*A11' I)(0  S)(0          I)
 *
 * with I the identity matrix and S=A22-B*inv(A11)*B'
 *
 * Let inv(A) be a preconditioner
 *
 *                 inv(A) = inv(U)inv(D)inv(L) = (I -inv(A11)B')(inv(A11)      0)(I           0)
 *                                               (0           I)(0       inv(S))(-B*inv(A11) I)
 *
 * This class implements an approximation of inv(A) where S=A22
 *
 *                 inv(P) ~ inv(A)
 *
 *                 inv(P) = (I   -inv(A11)B')(inv(A11)         0)(I           0)
 *                          (0             I)(0       inv(S=A22))(-B*inv(A11) I)
 *
 *                 inv(P) = (inv(A11)  -inv(A11)B'inv(S=A22))(I             0)
 *                          (0                    inv(S=A22))(-B*inv(A11)   I)
 *
 * In order to compute [y1 y2]' = inv(P)[x1 x2]' we do
 *
 *                 z  = inv(A11)*x1
 *                 y2 = inv(A22)*(x2 - B*z)
 *                 y1 = z - inv(A11)(B*y2)
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
class PCLDUFactorisation
{
public:

    /**
     * This struct defines the state of the preconditioner (initialised data and objects to be reused).
     */
    typedef struct{
        Mat A11_matrix_subblock; /**< Mat object that stores the A11 subblock.*/
        Mat A22_matrix_subblock; /**< Mat object that stores the A22 subblock.*/
        Mat B_matrix_subblock; /**< Mat object that stores the B subblock.*/
        PC  PC_amg_A11; /**<  inv(A11) is approximated by an AMG cycle. We compute it with HYPRE via a PC object. See \todo - don't create this every iteration but save it first time is needed. */
        PC  PC_amg_A22; /**<  inv(A22) is approximated by an AMG cycle. We compute it with HYPRE via a PC object. See \todo - don't create this every iteration but save it first time is needed. */
        Vec x1_subvector;/**<  Used to store the first half of the vector to be preconditioned*/
        Vec x2_subvector;/**<  Used to store the second half of the vector to be preconditioned*/
        Vec y1_subvector;/**<  Used to store the first half of the preconditioned vector*/
        Vec y2_subvector;/**<  Used to store the second half of the preconditioned vector*/
        Vec z;/**<    Used to store intermediate results*/
        Vec temp;/**< Used to store intermediate results*/
        VecScatter A11_scatter_ctx;/**< Scattering context: gather x1 from x and scatter y1 back into y*/
        VecScatter A22_scatter_ctx;/**< Scattering context: gather x2 from x and scatter y2 back into y*/
#ifdef TRACE_KSP
        double mScatterTime;/**< Time counter used for profiling scatter operations*/
        double mA1PreconditionerTime;/**< Time counter used for profiling the application of the preconditioner on the A11 block*/
        double mA2PreconditionerTime;/**< Time counter used for profiling the application of the preconditioner on the A22 block*/
        double mExtraLAOperations;/**< Time counter used for profiling extra matrix-vector and daxpy operations*/
        double mGatherTime;/**< Time counter used for profiling gather operations*/
#endif

    } PCLDUFactorisationContext;

    PCLDUFactorisationContext mPCContext; /**< PC context, this will be passed to PCBlockDiagonalApply when PETSc returns control to our preconditioner subroutine.  See PCShellSetContext().*/
    PC mPetscPCObject;/**< Generic PETSc preconditioner object */

public:

    /**
     * Constructor.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     */
    PCLDUFactorisation(KSP& rKspObject);

    ~PCLDUFactorisation();

private:

    /**
     * Creates all the state data required by the preconditioner.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     */
    void PCLDUFactorisationCreate(KSP& rKspObject);

    /**
     * Setups preconditioner.
     */
    void PCLDUFactorisationSetUp();
};

#endif /*PCLDUFACTORISATION_HPP_*/
